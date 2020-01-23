#!/usr/bin/env python3

"""
Python 3 Based 

http://compas.science

Part of COMPAS-post-Processing

"""
import h5py as h5
import os


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
