#!/usr/bin/python3
import h5py  as h5      # for reading and writing h5 format
import numpy as np      # for handling arrays
import os                   # for directory walking
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
dataRootDir    = '.'    # Location of root directory of the data - 
                        # user should change this if not running from root dir
prefix         = ''     # Prefix of the data files shown here is default
delimiter      = ','    # Delimeter used in the output csv files -
                        # but can be set differently if desired
extension      = 'csv'  # Extension of the data files

 #Which files to combine current options are
filesToCombine = ['SystemParameters'    ,'CommonEnvelopes',\
                  'DoubleCompactObjects','Supernovae'       ]

h5Name         =  'COMPAS_output.h5' #Name of the output file

###############################################################
#
#     Changing code below this line is at own risk
# 
################################################################


# Number of lines in the data file headers - Probably will never change, 
# but this avoids having "magic numbers" below
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
    combineOutputsOfFile(dataRootDir=dataRootDir, h5GroupDict=h5GroupDict,\
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
    print('Done, :smiling_imp:')

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

    'CommonEnvelopes'     : str(prefix) + 'Common_Envelopes.' + str(extension),\
    'DoubleCompactObjects': str(prefix) + 'Double_Compact_Objects.' + str(extension),\
    'Supernovae'          : str(prefix) + 'Supernovae.' + str(extension),\
    'SystemParameters'    : str(prefix) + 'System_Parameters.' + str(extension)
    }
    ##Future implementation
    #,
    #    'RLOF'                    : str(prefix) + '_RLOF.' + str(extension),\
    #    'errors'                : str(prefix) + 'errorfile.' + str(extension),\
    #    'output'                : str(prefix) + 'output.' + str(extension)\
    #    }

    #Create empty dictionary
    h5GroupDict = {}

    #Fill it in with only the files you want
    for f in filesToCombine:
        if f in optionsDict.keys():
            h5GroupDict[f] = optionsDict[f]
        else:
            raise ValueError("%s is not a group that exists. \n\
                             Currently we include %s "%(f, optionsDict.keys()))


    return h5GroupDict





##################################################################
###
### Step 1: Check that the rootDataDir exists and correct formatting if necessary
###
##################################################################

def verifyPathsFiles(dataRootDir=None, h5GroupDict=None):

    #Test if root directory exists and ensure the path ends with a '/'
    if dataRootDir[-1] != "/":
        dataRootDir = dataRootDir + "/"    

    if not os.path.isdir(dataRootDir):
            raise ValueError("directory not found with path: %s"%(dataRootDir))


    #Test if the files exist in any of the subdirectories

    # for every file that you which to combine and create h5
    for groupName in h5GroupDict.keys():
        fileName = h5GroupDict[groupName]
        # go through all the subdirectories and see if it exists
        Exist  = False
        for root,dirs,files in os.walk(dataRootDir):
            for f in files:
                if f == fileName:
                    Exist = True
        # If it does not exist inform user and stop, user should delete group
        if not Exist:
            raise ValueError("%s not found in rootDir: %s"\
                              %(fileName, dataRootDir))  
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

    for compasDataFilename in h5GroupDict.values():
                         
        ################################################################
        # Initialize variables for safety checks

        isHeaderWritten = False        # Boolean to ensure header is 
                                    # only written once
        nColumnCheck    = None         # Check that column numbers are 
                                    # consistent for each category


        ############################################################################
        # Iterate through each subdirectory to find all the relevant output files 
        #- if none exists, no output file will be produced
        for root,dirs,files in os.walk(dataRootDir):
            for f in files:
                if f == compasDataFilename:

                    path = os.path.join(root, f)

                    #individual output file of run in subfolder
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

                            if i == 0:                         # first row - set the required column number
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

    
        # If combine file does not exist, skip it - 
        #this happens if a category of output (e.g RLOF) does not occur in a given run
        #this should not happen since we have checked for this before
        if not os.path.isfile(combineFilePath):
            raise ValueError("trying to read %s, but does not exist." %(combineFilePath))
    
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
        Data = h5.File(h5data)

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
### With all functions defined, now run the whole script
###
##################################################################

if __name__ == "__main__":
    #If you run this script from a terminal
    #Use the global parameters defined at start script
    #Otherwise call it from your own script with the settings there
    main(filesToCombine=filesToCombine, dataRootDir=dataRootDir, \
         prefix=prefix, delimiter=delimiter, extension=extension, \
         h5Name=h5Name)




