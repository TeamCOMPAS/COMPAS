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
    i.e. giveMeFiles=[1,2,3,4]

    The numbers relate to
    1:Compas_Log_BSE_Common_Envelopes.csv'
    2:Compas_Log_BSE_Double_Compact_Objects.csv'
    3:Compas_Log_BSE_Supernovae.csv'
    4:Compas_Log_BSE_System_Parameters.csv

3 - run the script from terminal (commands are executed at the bottom of this file)

The output will be in the same folder as the path

"""

pathToData  = '/home/cneijssel/Desktop/IlyaData/VignaGomezRepeat'
giveMeFiles = [1,2,3,4]
pathData    = '.' #absolute path :)





###############################################################
#
#
#
#     Changing code below this line is at own risk
#
#
#
################################################################







fileNames = {1:'Compas_Log_BSE_Common_Envelopes.csv',\
             2:'Compas_Log_BSE_Double_Compact_Objects.csv',\
             3:'Compas_Log_BSE_Supernovae.csv',\
             4:'Compas_Log_BSE_System_Parameters.csv'}
#Future implementations
#             5:'Compas_Log_BSE_RLOF.csv',\  
#             6:'errorfile',\
#             7:'output',\
#              }



#How the groups will be called in H5

groupNames = {1:'CommonEnvelopes',\
              2:'DoubleCompactObjects',\
              3:'Supernovae',\
              4:'SystemParameters'}
#Future implementations
#             5:'RLOF',\  
#             6:'errors',\
#             7:'output'}
def combineOutputsOfFile(baseDirectoryData='.', groups=[1,2,3,4]):
    """
    For a simulation in folder baseDirectory
    1 - write a new file name=Combined_'filename'
    2 - go through all the subfolders and find the filename
    3 - Of the first file with filename, copy the header 
        and write the data to Combined_'filename'
    4 - Of other files with filename, copy the data
    5 - close the newly written file
    
    input = filename , string(with extension)
    
    COMPAS has a strict header format which is assumed here
    
    1st line = type data (INT, FLOAT BOOL etc)
    2nd line = unit data (Msol, ergs,Tsol etc)
    3th line = column name
    """
    for group in groups:
        filename       = fileNames[group]
                         
        #1-------
        combinedOutput = open(baseDirectoryData+'/Combine_'+filename, 'w')
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


                            line  = line.replace(" ", "")
                            line  = line.replace(",", "\t")
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
                        nCols = len(line.split(','))
                        line  = line.replace(",", "\t")
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
        fileName   = '/Combine_'+fileNames[groupNumber]
        filePath   = baseDirectoryData+fileName
        hf.create_group(groupName)
        addHdf5HeadersAndAttributes(hf, groupName, filePath)
        addHdf5Data(hf, groupName, filePath)
    hf.close()

    

def cleanUpInIsleNumber2Please(baseDirectoryData='.', groups=[1,2,3,4]):
    
    for groupNumber in groups:
        fileName   = '/Combine_'+fileNames[groupNumber]
        filePath   = baseDirectoryData+fileName
        command    = 'rm '+filePath
        sp.Popen(command, shell=True, executable='/bin/bash')


combineOutputsOfFile(baseDirectoryData=pathToData, groups=giveMeFiles)
print('Combining csv files from subdirectories')
createH5file(baseDirectoryData=pathToData, groups=giveMeFiles)
print('Creating H5 file')
cleanUpInIsleNumber2Please(baseDirectoryData=pathToData, groups=giveMeFiles)
print('Clean up the combined CSV files')
print('Done, :smiling_imp:')
