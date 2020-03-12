#!/usr/bin/python3
import h5py  as h5  #for reading and writing h5 format
import numpy as np  #for handling arrays
import os           #for directory walking
import subprocess as sp #for executing terminal command from python


"""
This script turns the different outputs from a COMPAS simulation into 
a single h5file. If you have a simulation in different subfolders due
to a large run, set the path to the data to the main folder it is in.
Note that in order to create the data we combine the data from different
folders in a single csv file first. Hence you need enough space to
store a duplicate of the data. The combined csv files will automatically
be removed afterwards.

The steps:


* - set the path to the data
* - create a list with the number of each file you want
    i.e. giveMeFiles=[1,2,3,4,5]

    The numbers relate to

    1:Prefix_Common_Envelopes.xxx
    2:Prefix_Double_Compact_Objects.xxx
    3:Prefix_Supernovae.xxx
    4:Prefix_System_Parameters.xxx
    5:Prefix_RLOF.txt

    where both the prefix and the extension xxx
    are set in the python submit

* - run the script from terminal in which case you can
    set the prefix/extension/path/files below

  - Or pass the parameters through combineCreateH5file() function
    from notebook or terminal see postProcessing/H5/1_creatingH5.html

The output will be in the same folder as the path or the folder 
the script is placed in.

"""

prefix      = None
extension   = None

pathToData  = './'
giveMeFiles = [1,2,3,4,5]





###############################################################
#
#
#
#     Changing code below this line is at own risk
#
#
#
################################################################







fileNames = {1:'_Common_Envelopes.',\
             2:'_Double_Compact_Objects.',\
             3:'_Supernovae.',\
             4:'_System_Parameters.',\
             5:'_RLOF.'}
#Future implementations
#             5:'Compas_Log_BSE_RLOF.csv',\  
#             6:'errorfile',\
#             7:'output',\
#              }



#How the groups will be called in H5

groupNames = {1:'CommonEnvelopes',\
              2:'DoubleCompactObjects',\
              3:'Supernovae',\
              4:'SystemParameters',\
              5:'RLOF'}
#Future implementations
#             5:'RLOF',\  
#             6:'errors',\
#             7:'output'}
def combineOutputsOfFile(baseDirectoryData='.', groups=[1,2,3,4], \
    prefix=None, extension=None):
    """
    For a simulation in folder baseDirectory
    1 - write a new file name=Combined_'filename'
    2 - go through all the subfolders and find the filename
    3 - Of the first file with filename, copy the header 
        and write the data to Combined_'filename'
    4 - Of other files with filename, copy the data
    5 - close the newly written file

    The numbers listed above are also comments in code here
    #x----

    
    
    COMPAS has a strict header format which is assumed here
    
    1st line = type data (INT, FLOAT BOOL etc)
    2nd line = unit data (Msol, ergs,Tsol etc)
    3th line = column name
    4th ++++ = data
    """
    for group in groups:
        #add prefix and extension
        filename       = prefix+fileNames[group]+extension
        #1-------
        combinedOutput = open(baseDirectoryData+'/Combine'+fileNames[group]+'txt', 'w')
        #boolean to see if we have written the header
        headersWritten = False
        nHeaders       = 3
        nColumnCheck   = None #check each line if nr of entries is the same
                              #if not there is somthing wrong
        #2---- 
        for root,dirs,files in os.walk(baseDirectoryData):

            for f in files:

                if f == filename:

                    path = os.path.join(root, f)

                    #individual output file of run in subfolder
                    outputFile   = open(path)
                    #3--------------
                    if not headersWritten:
                        for i in range(nHeaders):
                            line = outputFile.readline()

                            #clear whitespace
                            line  = line.replace(" ", "")
                            #break line on special character an turn into tab
                            if extension == 'txt':
                                #tab is already special character
                                pass
                            elif extension == 'csv':
                                line  = line.replace(",", "\t")
                            else:
                                raise ValueError("extension is not recognized (csv, txt)")
                            
                            # check number of columns 
                            nCols = len(line.split('\t'))
                            if i ==0: #set the column number check for first line
                                nColumnCheck = nCols
                            if(nCols != nColumnCheck):
                                raise ValueError('wrong number of columns in header=%s'%(i))
                            combinedOutput.write(line)
                        headersWritten = True
                    else:
                        #skip the header by reading and not doing anything
                        [outputFile.readline() for i in range(nHeaders)]

                    #4 -----------
                    for line in outputFile:
                        #break line on special character an turn into tab
                        if extension == 'txt':
                            #tab is already special character
                            pass
                        elif extension == 'csv':
                            line  = line.replace(",", "\t")
                        else:
                            raise ValueError("extension is not recognized (csv, txt)")
                        # check number of columns 
                        nCols = len(line.split('\t'))    
                        if(nCols != nColumnCheck):
                            raise ValueError('wrong number of columns in data')
                        combinedOutput.write(line)
        #5----------
        combinedOutput.close()
    




def addHdf5HeadersAndAttributes(hf,  groupName, filePath):
    """
    COMPAS has a strict header format which is assumed here
    
    1st line = type data (INT, FLOAT BOOL etc)
    2nd line = unit data (Msol, ergs,Tsol etc)
    3th line = column name
    """
    
    file      = open(filePath, 'r')
    #get header, units names
    types     = file.readline()[:-1].split('\t')
    units     = file.readline()[:-1].split('\t')
    headers   = file.readline()[:-1].split('\t')
    #how many entries will a column in the group have?
    #get the length of the file (minus headers)
    fileLength = int(sp.check_output('wc -l ' + filePath, \
                                     shell=True).split()[0]) - 3
    file.close() # only needed the headers here
    #types is strings need to replace by actual type for h5
    dtypes    = []
    for nrt, typ in enumerate(types):
        if typ == 'INT':
            dtypes.append(np.int64)
        elif typ == 'FLOAT':
            dtypes.append(np.float64)
        elif typ == 'BOOL':
            dtypes.append(bool)
        else:
            raise ValueError("Unrecognised datatype typ=%s - for column %s in file%s "\
                             %(typ, headers[nrt], groupName))
    #create the groups in the h5file and add units and explanation string
    for header,dtype,unit in zip(headers,dtypes,units):
        dset = hf[groupName].create_dataset(header,dtype=dtype,shape=(fileLength,))
        dset.attrs['units'] = unit
        #dset.attrs['comment']= columnDescriptions
        
    return

def addHdf5Data(hf,  groupName, filePath):
    
    #too slow to go line by line, so load in a modest 
    #(in term sof memory) amount at a time
    chunkSize = 500000
    
    #get the length of the file (minus headers)
    fileLength = int(sp.check_output('wc -l ' + filePath, \
                                      shell=True).split()[0]) - 3

    file      = open(filePath)
    types     = file.readline()[:-1].split('\t')
    units     = file.readline()[:-1].split('\t')
    headers   = file.readline()[:-1].split('\t')
    
    nrColumns   = len(headers)
    group       = hf[groupName]
    chunkBegin = 0
    chunkEnd = 0    
    while chunkEnd < fileLength:
        data = []

        chunkEnd = chunkBegin + chunkSize

        #dont try to load in more data than you've got
        if chunkEnd > fileLength:
            chunkEnd = fileLength

        #read in a modest number of lines
        for i in range(chunkEnd-chunkBegin):
            data.append(file.readline()[:-1].split())
            
        data = np.array(data)
            
        #data is now a 2d array where each column is a specific variable

        
        for nrcolumn in range(nrColumns):
            #fill in the values in the preshaped array
            #(see dset addHdf5HeadersAndAttributes() )
            
            columnName        = headers[nrcolumn]
            dtype             = type(group[columnName][0])
            group[columnName][chunkBegin:chunkEnd] = np.array(data[:,nrcolumn],dtype=dtype)
            
        chunkBegin = chunkEnd


    
def createH5file(baseDirectoryData='.', groups=[1,2,3,4], h5Name='COMPAS_output.h5'):
    
    hf = h5.File(baseDirectoryData+'/'+h5Name, 'w')
    
    #use the groupNames dictionary to create 
    #the H5file name group. Each of these we will fill with the header and data
    for groupNumber in groups:
        groupName  = groupNames[groupNumber]
        fileName   = '/Combine'+fileNames[groupNumber]+'txt'
        filePath   = baseDirectoryData+fileName
        hf.create_group(groupName)
        addHdf5HeadersAndAttributes(hf, groupName, filePath)
        addHdf5Data(hf, groupName, filePath)
    hf.close()

    

def cleanUpInIsleNumber2Please(baseDirectoryData='.', groups=[1,2,3,4]):
    
    for groupNumber in groups:
        fileName   = '/Combine'+fileNames[groupNumber]+'txt'
        filePath   = baseDirectoryData+fileName
        command    = 'rm '+filePath
        sp.Popen(command, shell=True, executable='/bin/bash')




def printAllColumnsInH5(pathToData):
    """
    Function to print all files with their column names/length/unit
    Returns nothing, just prints to cell/terminal
    Most of the function is for nice spacing

    """

    #Check if a correct path is given
    if not  os.path.isfile(pathToData):
        raise ValueError("h5 file not found. Wrong path given?")
    elif os.path.isfile(pathToData):
        Data  = h5.File(pathToData)


    Files  = Data.keys()

    for File in Files:
        print('Filename = %s' %(File))
        print('----------------------')

        #Everytime you see Xr*' '
        #It means I add X spaces to line it
        print('\t   column name%sunit%slength'%(29*' ',16*' '))
        print('\t   '+'-----------------'*4)
        
        #In this file give me all the column names
        columns = Data[File].keys()
        
        #for every column in the columns
        for nrc,column in enumerate(columns):
            #always want the column name printed in 40 char
            spaces = ' '*(40 - len(column))
            length = Data[File][column].shape[0]
            #always want the unit name printed over 20 chars
            unit   = Data[File][column].attrs['units']
            spaces2 = ' '*(20 - len(unit))
            #--
            length = Data[File][column].shape[0]

            print('\t   %s%s%s%s%s'%(column,spaces, unit,spaces2, length))
            #Every 4 lines print a dashed line to read output easier
            if (nrc%5==4):
                print('\t   '+'-----------------'*4)
    Data.close()




def combineCreateH5file(pathToData=None, files=None, prefix=None, extension=None):
    print('Combining csv files from subdirectories')
    combineOutputsOfFile(baseDirectoryData=pathToData, groups=giveMeFiles, \
                         prefix=prefix, extension=extension)
    print('Creating H5 file')
    createH5file(baseDirectoryData=pathToData, groups=giveMeFiles)
    print('Clean up the combined CSV files')
    cleanUpInIsleNumber2Please(baseDirectoryData=pathToData, groups=giveMeFiles)
    print('Done, :smiling_imp:')
    print()
    i = 0
    while i < 5:
        i+=1
        if i == 3:
            print('#------contents H5 file-----#')
        else:
            print('#---------------------------#')
    print()
    printAllColumnsInH5(pathToData+'/COMPAS_output.h5')



if __name__ == "__main__":
    
    combineCreateH5file(pathToData=pathToData, files=giveMeFiles, prefix=prefix, extension=extension)
