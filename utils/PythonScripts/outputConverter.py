#!/usr/bin/python3
import h5py  as h5      # for reading and writing h5 format
import numpy as np      # for handling arrays
import os, sys          # for directory walking
import subprocess as sp # for executing terminal command from python

"""
This script turns the different outputs from a COMPAS simulation into 
a single h5file. If you have a simulation in different subfolders due
to a large run, set the path to the data to the output root directory.
Note that in order to create the data we combine the data from different
folders in a single csv file first. Hence you need enough space to
store a duplicate of the data. The combined csv files will automatically
be removed afterwards.
"""

### User-defined parameters
def setDefaults():

    prefix         = 'BSE_'  			# Prefix of the data files                   # defaults to 'BSE_'  
    delimiter      = ','                # Delimeter used in the output csv files     # defaults to ','        
    extension      = 'csv'              # Extension of the data files                # defaults to 'csv'
    h5Name         = 'COMPAS_Output.h5' # Name of the output h5 file				 # defaults to 'COMPAS_Output.h5' 

    # Location of root directory of the data can be specified as a command line argument    
    if len(sys.argv) > 1:                   
        dataRootDir = str(sys.argv[1])      # data directory specified on command line
    else:
        dataRootDir    = '.'                # defaults to '.'            
    
    # To only combine a subset of the data files, specify them here    
    filesToCombine = None    # default None means to use all of them (apologies if that's counterintuitive...)
    #filesToCombine = [\
    #    'SystemParameters',\
    #    'CommonEnvelopes',\
    #    'DoubleCompactObjects',\
    #    'Supernovae',\
    #    'RLOF',\
    #    'errors',\            
    #    'output'\
    #]

    return filesToCombine, dataRootDir, prefix, delimiter, extension, h5Name




###############################################################
###############################################################
#
#     Changing code below this line is at own risk
# 
###############################################################
################################################################


### Global vars

# Number of lines in the data file headers.
#    Probably will never change, but
#    this avoids having "magic numbers" below
nLinesInHeader = 3            





##################################################################
###
### Main Function
###
##################################################################

def main(filesToCombine=None, dataRootDir=None, prefix=None,\
         delimiter=None, extension=None, h5Name='COMPAS_output.h5'):

    """ 
    Combines all of the different output files of the same type (e.g RLOF, Supernovae, etc.)
    into one CSV file per type, then synthesizes all of the combined CSV files into one
    single H5 file, before cleaning up the repository.
    """
    ### Step 0: create dictionary filesToCombine and Paths
    h5GroupDict = createDictionaryGroupPaths(filesToCombine=filesToCombine,\
                                             prefix=prefix, extension=extension)


    ### Step 1: Check that the rootDataDir exists and correct formatting if necessary
    print('Checking if directory and files to combine exist')
    dataRootDir = verifyPathsFiles(dataRootDir=dataRootDir, h5GroupDict=h5GroupDict)


    ### Step 2: Create the combined CSV file for each output type
    print('Combining %s files from subdirectories' %(extension))
    setOfUsedDatafiles = combineOutputsOfFile(dataRootDir=dataRootDir, h5GroupDict=h5GroupDict,\
                         delimiter=delimiter) 
    

    ### Step 3: Create a single H5 file for all the data
    print('Creating H5 file %s' %(h5Name))
    createH5file(dataRootDir=dataRootDir, h5Name=h5Name, h5GroupDict=h5GroupDict) 
    

    ### Step 4: Remove the temporary files 
    print('Cleaning up the combined %s files' %(extension))
    cleanUpInAisleNumber2Please(dataRootDir=dataRootDir, h5GroupDict=h5GroupDict) 
    print()
    print('-------------------------')
    print('--Overview of your data--')
    print('-------------------------')

    ### Step 5: Print columns in the h5 file
    printH5Columns(dataRootDir=dataRootDir, h5Name=h5Name)

    ### Step 6: Print out which data files were used
    printUsedDataFiles(setOfUsedDatafiles)
    print('Done, :smiling_imp:\n')



##################################################################
###
### Step 0: Create dictionary which links groupnames and paths
###         of the files you which to combine
###
##################################################################

def createDictionaryGroupPaths(filesToCombine = None, prefix=None, \
                               extension=None):
    
    # The current groups we offer are 
    optionsDict = {
        'CommonEnvelopes'      : str(prefix) + 'Common_Envelopes.' + str(extension),\
        'DoubleCompactObjects' : str(prefix) + 'Double_Compact_Objects.' + str(extension),\
        'Supernovae'           : str(prefix) + 'Supernovae.' + str(extension),\
        'SystemParameters'     : str(prefix) + 'System_Parameters.' + str(extension),\
        'RLOF'                 : str(prefix) + 'RLOF.' + str(extension),\
        'errors'               : str(prefix) + 'errorfile.' + str(extension),\
        'output'               : str(prefix) + 'output.' + str(extension)\
    }

    # Create empty dictionary
    h5GroupDict = {}

    # Fill it in with only the files you want
    if filesToCombine == None:         # For default setting None, use all of the options
        h5GroupDict = optionsDict
    else:                             # If a subset of the data files is specified, use that
        for f in filesToCombine:
            if f in optionsDict.keys():
                h5GroupDict[f] = optionsDict[f]
            else:
                raise ValueError("%s is not a group that exists. \n\
                                 Currently we include %s "%(f, optionsDict.keys()))


    return h5GroupDict





##################################################################
###
### Step 1: Check that the rootDataDir exists 
###         and correct formatting if necessary
###
##################################################################

def verifyPathsFiles(dataRootDir=None, h5GroupDict=None):

    # Ensure the root directory path string ends with a '/'
    if dataRootDir[-1] != "/":
        dataRootDir = dataRootDir + "/"    

    # Throw an error if root directory does not exist
    if not os.path.isdir(dataRootDir):
        raise ValueError("directory not found with path: %s"%(dataRootDir))

    # Return updated dataRootDir
    return dataRootDir
        


##################################################################
###
### Step 2: Create the combined CSV file for each output type
###
##################################################################

def combineOutputsOfFile(dataRootDir=None, h5GroupDict = None, delimiter=None):
    """
    COMPAS has a strict header format which is assumed here
    
    1st line = type data (INT, FLOAT BOOL etc)
    2nd line = unit data (Msol, ergs,Tsol etc)
    3th line = column name
    4th ++++ = data

    This function identifies which categories of output file
    were produced and need to be included in the H5 file. Then
    it collects all the data for a given category and puts it
    in a CSV file with the correct header.
    """

    # Keep track of which data files are picked up by the walker
    setOfUsedDatafiles = set()     # Set of filenames found in output

    for compasDataFilename in h5GroupDict.values():
                         
        ################################################################
        # Initialize variables for safety checks

        isHeaderWritten = False        # Boolean to ensure header is 
                                       # only written once
        nColumnCheck    = None         # Check that column numbers are 
                                       # consistent for each category


        ############################################################################
        # Iterate through each subdirectory to find all the relevant output files 
        #     if none exists, no output file will be produced
        for root,dirs,files in os.walk(dataRootDir):
            for f in files:
                if f == compasDataFilename:

                    # Add to set of discovered datafiles
                    setOfUsedDatafiles.add(f)    

                    # Open the file
                    path = os.path.join(root, f)
                    compasData = open(path)

                    #######################################################################
                    # For the first processed output file, include the header as well

                    if not isHeaderWritten:

                        # Only write the header once
                        isHeaderWritten = True

                        # Create the empty combine file now, so that data can be appended later    
                        with open(dataRootDir+'Combine_'+compasDataFilename, 'w') as combineWrite:
                            pass


                        ######################################################
                        # Read in the appropriate number of header lines
                        for i in range(nLinesInHeader):
                            line = compasData.readline()


                            # Verify that the number of columns is consistent across rows
                            nCols = len(line.split(delimiter))

                            if i == 0:                       # first row - set the required column number
                                nColumnCheck = nCols         

                            else:                            # later rows - verify the column number
                                if nCols != nColumnCheck:    
                                    raise ValueError('wrong number of columns in header=%s'%(i))


                            # Clean up the lines, and write to the output
                            line  = line.replace(" ", "")                # remove whitespace
                            line  = line.replace(delimiter, "\t")        # swap input delimiter with a tab (to simplify writing to h5 later)

                            # Write the header to file
                            with open(dataRootDir + 'Combine_' + compasDataFilename, 'a') as combineWrite:
                                combineWrite.write(line)

                    else:
                        # For later output files, skip the header by reading and not doing anything
                        [compasData.readline() for i in range(nLinesInHeader)]

                    ###################################################
                    # Process the non-header lines of each file
                    for line in compasData:


                        ## Verify that the number of columns is consistent across rows
                        nCols = len(line.split(delimiter))

                        if nCols != nColumnCheck:
                            raise ValueError('wrong number of columns in data')


                        # Clean up the lines, and write to the output
                        line  = line.replace(" ", "")                # remove whitespace         # Is this necessary? Coen didn't have this...
                        line  = line.replace(delimiter, "\t")        # swap input delimiter with a tab (to simplify writing to h5 later)

                        with open(dataRootDir + 'Combine_' + compasDataFilename, 'a') as combineWrite:
                            combineWrite.write(line)

    # Return the set of all used datafiles to be printed at the end
    return setOfUsedDatafiles
    



##################################################################
###
### Step 3: Create a single H5 file for all the data
###
##################################################################

def createH5file(dataRootDir=None, h5GroupDict=None, h5Name='COMPAS_output.h5'): 
    """
    Function to create the h5 file, extract the details of the 
    Combine_ files, and call the functions to fill the h5 with 
    headers and data.
    """
    
    hf = h5.File(dataRootDir + h5Name, 'w')
    
    # Use the h5GroupDict dictionary to find the relevant Combine file for a given group, then fill the h5 group with the header and data
    for group in h5GroupDict.keys():
        combineFilePath = dataRootDir + 'Combine_' + h5GroupDict[group]

    
        # If combine file does not exist, skip it 
        #     This happens if a category of output (e.g RLOF) does not occur in 
        #     any systems of a given run
        if not os.path.isfile(combineFilePath):
            continue
    
        # Create h5 group for the given category, and enter in the header and data
        hf.create_group(group)
        addHdf5HeadersAndAttributes(hf, group, combineFilePath)
        addHdf5Data(hf, group, combineFilePath)

    hf.close()


##################################################################

def addHdf5HeadersAndAttributes(hf,  group, filePath):
    """
    COMPAS has a strict header format which is assumed here
    
    1st line = data type (INT, FLOAT, BOOL, etc)
    2nd line = parameter unit (Msol, ergs, Tsol, etc)
    3th line = parameter name
    """
    
    # Extract header information
    with open(filePath, 'r') as fileRead:
    
        types  = fileRead.readline()[:-1].split('\t')
        units  = fileRead.readline()[:-1].split('\t')
        params = fileRead.readline()[:-1].split('\t')

    # Extract number of systems in file (total file length - header)
    fileLength = int(sp.check_output('wc -l ' + filePath, shell=True).split()[0]) - nLinesInHeader        

    # Need to replace dataType, which is a string, by actual type for h5
    dtypes    = []
    for iType, dataType in enumerate(types):
        if dataType == 'INT':
            dtypes.append(np.int64)
        elif dataType == 'FLOAT':
            dtypes.append(np.float64)
        elif dataType == 'BOOL':
            dtypes.append(bool)
        elif dataType == 'STRING':
            dtypes.append(h5.string_dtype(encoding='utf-8'))
        else:
            raise ValueError("Unrecognised datatype dataType=%s - for column %s in file%s "\
                             %(dataType, params[iType], group))

    # Create the groups in the h5file and add units 
    for param, dtype, unit in zip(params, dtypes, units):
        dset = hf[group].create_dataset(param, dtype=dtype, shape=(fileLength,))
        dset.attrs['units'] = unit
        
    return


##################################################################

def addHdf5Data(hf,  group, filePath):
    """
    Function to append data from Combine_ files
    into the h5 file
    """
    
    # Too slow to go line by line, so load in a modest (in terms of memory) amount at a time
    chunkSize = 500000    # Binary systems
    
    #get the length of the file (minus headers)
    fileLength = int(sp.check_output('wc -l ' + filePath, shell=True).split()[0]) - nLinesInHeader

    # Open the file so that lines can be read in succession
    with open(filePath, 'r') as fileRead:

        # Read the header lines first so they don't interfere with the data
        types  = fileRead.readline()[:-1].split('\t')
        units  = fileRead.readline()[:-1].split('\t')
        params = fileRead.readline()[:-1].split('\t')
        
        # Initialize parameters
        h5group       = hf[group]
        chunkBegin = 0
        chunkEnd = 0    

        # Loop over the file in chunkSize'd chunks
        while chunkEnd < fileLength:
            data = []
    
            chunkEnd = chunkBegin + chunkSize
    
            # Don't try to load in more data than you've got
            if chunkEnd > fileLength:
                chunkEnd = fileLength
    
            # Read in a modest number of lines
            for i in range(chunkEnd-chunkBegin):
                data.append(fileRead.readline()[:-1].split())
                
            data = np.array(data)     # data is now a 2d array where each column is a specific variable
    
            # Add the data array into the h5 file    
            for iParam, param in enumerate(params):
                
                dtype             = type(h5group[param][0])
                h5group[param][chunkBegin:chunkEnd] = np.array(data[:,iParam],dtype=dtype)
        
            # Leapfrog to the next chunk location
            chunkBegin = chunkEnd
    


    
##################################################################
###
### Step 4: Remove the temporary files 
###
##################################################################

def cleanUpInAisleNumber2Please(dataRootDir='./',h5GroupDict=None): 
    """
    Function to remove the temporary Combine_* files 
    """
    
    for group in h5GroupDict.keys():

        combineFilePath = dataRootDir + 'Combine_' + h5GroupDict[group]

        # If combine file does not exist, skip it
        if not os.path.isfile(combineFilePath):
            continue

        shellCommand = 'rm ' + combineFilePath
        sp.Popen(shellCommand, shell=True, executable='/bin/bash')
        



##################################################################
###
### Step 5: Print columns in the h5 file
###
##################################################################

def printH5Columns(dataRootDir='./', h5Name="COMPAS_output.h5"):
    """
    Function to print all files with their column names/length/unit
    Returns nothing, just prints to cell/terminal
    Most of the function is for nice spacing
    """

    # Check that path is correct
    h5data = dataRootDir + h5Name
    if not os.path.isfile(h5data):
        raise ValueError("h5 file not found. Wrong path given?")
    else:
        Data = h5.File(h5data, 'r')

    Files = Data.keys()

    for File in Files:
        print('Filename = %s' %(File))
        print('----------------------')

        # Note:   X*' '   means X spaces in line
        print('\t   column name%sunit%slength'%(29*' ',16*' '))
        print('\t   '+'-----------------'*4)
        
        # In this file give me all the column names
        columns = Data[File].keys()
        
        # For every column in the columns
        for nrc, column in enumerate(columns):

            # Always want the column name printed in 40 char
            spaces = ' '*(40 - len(column))
            length = Data[File][column].shape[0]

            # Always want the unit name printed over 20 chars
            unit   = Data[File][column].attrs['units']
            spaces2 = ' '*(20 - len(unit))

            #--
            length = Data[File][column].shape[0]

            print('\t   %s%s%s%s%s'%(column,spaces, unit,spaces2, length))

            # Every 4 lines print a dashed line to read output easier
            if (nrc%5==4):
                print('\t   '+'-----------------'*4)

    Data.close()


##################################################################
###
### Step 6: Print out which data files were used
###
##################################################################

def printUsedDataFiles(setOfUsedDatafiles={}):
    """
    Last step: print out the set of all data files
    which were included in the H5. Explicitly, this is the
    intersection of the set of desired data files specified
    by the user and the set of outputted data files from 
    COMPAS (since small runs may not produce all the desired
    output files)
    """

    print("\n###########################################################################\n")
    print("\tThe COMPAS datafiles combined into this HDF5 file are:\n")
    [print("\t" + str(datafile)) for datafile in setOfUsedDatafiles]
    print("\n###########################################################################\n")




##################################################################
### 
### With all functions defined, now run the whole script
###
##################################################################

if __name__ == "__main__":
    # If you run this script from a terminal
    # Use the global parameters defined at the top
    # Otherwise call it from your own script with the settings there

    # Only pull the default settings above if this script is run from the command line
    filesToCombine, dataRootDir, prefix, delimiter, extension, h5Name = setDefaults()

    # Run the script above with the settings defined at the top    
    main(filesToCombine=filesToCombine, dataRootDir=dataRootDir, \
         prefix=prefix, delimiter=delimiter, extension=extension, \
         h5Name=h5Name)




