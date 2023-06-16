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
    
def reduceH5(pathToOld = None, pathToNew = None, dictColumns=None, dictSeeds=None):


        if pathToOld is None:
            raise ValueError("pathToOld not given")
        if pathToNew is None:
            raise ValueError("pathToNew not given")
        #if dictFiles is None:
            #raise ValueError("dictionary of files is not given")
        if dictColumns is None:
            raise ValueError("dictionary of columns is not given")
        if dictSeeds is None:
            raise ValueError("dictionary of seeds is not given")
        if set(dictColumns.keys()) != set(dictSeeds.keys()):
            raise ValueError("The column and seed dictionaries do not agree on the set of files (keys)")


        # This check is not needed if the last check above passes
        #if((len(dictFiles) != len(dictColumns)) or\
        #   (len(dictSeeds) != len(dictColumns)) or\
        #   (len(dictFiles) != len(dictSeeds))):									# RTW 3/25/20 - do we need the third statement?
        #    raise ValueError("The 3 dictionaries are not of same length!")

        h_old = h5.File(pathToOld, 'r')
        h_new = h5.File(pathToNew, 'w')


        #keys = dictFiles.keys()

        for filename in dictColumns.keys():
            #fileName = key
            sanityChecks(h_old, filename, dictColumns[filename], dictSeeds[filename])
            createDataInNewH5(h_old, h_new, filename, dictColumns[filename], dictSeeds[filename])
            # RTW
            
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
        Data  = h5.File(pathToData, 'r')


    Files  = Data.keys()

    for File in Files:
        print()
        print('Filename = %s' %(File))
        print('----------------------')

        #Every time you see Xr*' '
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
