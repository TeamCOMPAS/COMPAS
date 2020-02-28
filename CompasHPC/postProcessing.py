import numpy as np
import os
import gc
import time
import subprocess as sp
import h5py as h5

###Options####
collectAndCombineOutput = True
runPlottingRoutines     = False #See plottingRoutines.py for options
cleanUp                 = False #True # Whether to tar individual output folders
h5OutputFileName        = "COMPASOutput.h5"
giveMeFiles             = [1,2,3,4] # See below for file names

#-- Base directory where results are
baseDirectory = os.getcwd()

#-- Where to output post-processed results
rootOutputDir = os.getcwd()

##############

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--masterFolderDir", help="Directory to masterFolder")
args = parser.parse_args()

if(runPlottingRoutines):
    from plottingRoutines import plottingRoutines

compasRootDir = os.environ.get('COMPAS_ROOT_DIR')

sourceDirectory = os.path.join(compasRootDir, 'sr/COMPAS')

#make the HDF5 file
h5FilePath = os.path.join(rootOutputDir, h5OutputFileName)
h5File = h5.File(h5FilePath, 'w')

def copyPythonSubmit(baseDirectory="./", h5FilePath=h5FilePath):
    """
    """
    #put the contents of the python submit in there
    ps = open(os.path.join(baseDirectory, "pythonSubmit.py")) #os.path.join(baseDirectory, 'masterFolder/pythonSubmit.py'))
    #makes one long string and stores it, with \n linebreaks included
    pythonSubmit = ''
    for l in ps:
        pythonSubmit += l
    ps.close()

    h5File = h5.File(h5FilePath, 'a')

    #psDset = h5File.create_dataset('pythonSubmit', dtype="S" + str(len(pythonSubmit)), shape=(1,))
    #psDset[0] = pythonSubmit 
    h5File.attrs.create('pythonSubmit', pythonSubmit, dtype=h5.special_dtype(vlen=str))

    h5File.close()

    return

def collectSource():
    """
    """

    h5File = h5.File(h5FilePath, 'a')

    #put the c++ source in there
    sourceGroup = h5File.create_group('cppSource')

    for root,dirs,files in os.walk(sourceDirectory):
        for f in files:
            if f[-2:] == '.h' or f[-4:] == '.cpp' or f[-4:] == '.hpp':

                #makes one long string and stores it, with \n linebreaks included
                src = ''
                path = os.path.join(root,f)
                fil = open(path)
                for l in fil:
                    src += l
                fil.close()
                #srcDset = sourceGroup.create_dataset(f,dtype="S"+str(len(src)),shape=(1,))
                #srcDset[0] = src
                sourceGroup.attrs.create(f, src, dtype=h5.special_dtype(vlen=str))
    
    h5File.close()

    return
 
def cleanUpOutputFolders():
    """
    Tar up individual output folders to save space (in particular, inodes)
    """
    os.chdir(rootOutputDir)
    
    tarCommand = "tar czf tarredFolders.tar.gz output*"
    
    print(tarCommand)
    
    os.system(tarCommand)

    rmCommand = "rm -r ./output*"
    
    print(rmCommand)

    os.system(rmCommand)

    return

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

    return
    
def createH5file(baseDirectoryData='.', groups=[1,2,3,4], h5FilePath='COMPAS_output.h5'):

    hf = h5.File(h5FilePath, 'a')
    
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

    return

def cleanUpInIsleNumber2Please(baseDirectoryData='.', groups=[1,2,3,4]):
    
    for groupNumber in groups:
        fileName   = 'Combine_' + fileNames[groupNumber]
        filePath   = os.path.join(baseDirectoryData, fileName)
        command    = 'rm ' + filePath

        print("fileName = ", fileName)
        print("filePath = ", filePath)
        print("Removing ", filePath)
        print(command)
        proc = sp.Popen(command, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, executable='/bin/bash')
        
        ## Execute command
        if (sys.version_info > (3, 0)):
            proc.stdin.write(command.encode('utf-8'))
        else:
            proc.stdin.write(command)

        out, err = proc.communicate()

    return

#########################################################################################

if collectAndCombineOutput:
    combineOutputsOfFile(baseDirectoryData=rootOutputDir, groups=giveMeFiles)
    print('Combining csv files from subdirectories')
    createH5file(baseDirectoryData=rootOutputDir, groups=giveMeFiles, h5FilePath=h5FilePath)
    print('Creating H5 file')
    cleanUpInIsleNumber2Please(baseDirectoryData=rootOutputDir, groups=giveMeFiles)
    print('Clean up the combined CSV files')
    print('Done, :smiling_imp:')

copyPythonSubmit(baseDirectory=args.masterFolderDir)
collectSource()

if(runPlottingRoutines):
    plottingRoutines()

if(cleanUp):
    cleanUpOutputFolders()

print("post-processing completed successfully")
