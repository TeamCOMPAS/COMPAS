#!/usr/bin/python3
import h5py  as h5
import numpy as np
import os



def sanityChecks(h_old, filename, columns, seeds):

    data = h_old[filename]

    if columns == ['All']:
        columns = list(data.keys())

    for column in columns:
        try:
            test = data[column]
        except:
            raise ValueError("column %s does not exist in %s"\
                             %(column, filename))

    seedData = data['SEED'][()]

    mask      = np.in1d(seeds, seedData)
    not_exist = np.logical_not(mask)
    if np.sum(not_exist) !=0:
        raise ValueError("seed(s) %s do not exist in %s"\
                         %(seeds[not_exist], filename))


def createDataInNewH5(h_old, h_new, filename, columns, seeds):

    h_new.create_group(filename)

    dataOld = h_old[filename]
    mask    = np.in1d(dataOld['SEED'][()],seeds )

    if columns == ['All']:
        columns = list(dataOld.keys())

    for column in columns:

        unit       = dataOld[column].attrs['units']
        columnOld  = dataOld[column][()]
        data       = columnOld[mask]
        dataNew    = h_new[filename].create_dataset(column, data=data)
        dataNew.attrs['units'] = unit
        #add attribute comment
    
def reduceH5(pathToOld = None, pathToNew = None, dictFiles=None, dictColumns=None, dictSeeds=None):


        if pathToOld is None:
            raise ValueError("pathToOld not given")
        if pathToNew is None:
            raise ValueError("pathToNew not given")
        if dictFiles is None:
            raise ValueError("dictionary of files not given is None")
        if dictColumns is None:
            raise ValueError("dictionary of columns not given is None")
        if dictColumns is None:
            raise ValueError("dictionary of seeds not given is None")

        if((len(dictFiles) != len(dictColumns)) or\
           (len(dictSeeds) != len(dictColumns)) or\
           (len(dictFiles) != len(dictSeeds))):
            raise ValueError("the 3 dictionaries are not of same lengt!")

        h_old = h5.File(pathToOld)
        h_new = h5.File(pathToNew, 'w')


        keys = dictFiles.keys()

        for key in keys:
            fileName = dictFiles[key]          
            sanityChecks(h_old, dictFiles[key], dictColumns[key], dictSeeds[key])
            createDataInNewH5(h_old, h_new, dictFiles[key], dictColumns[key], dictSeeds[key])
            
        h_old.close()
        h_new.close()


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
