'''

h5view.py 

This program displays summary, header and content information for specified COMPAS HDF5 file(s).
It's fairly rudimetary - the HDF5 package provide h5dump and h5ls which have far more options
than this program - but this program is somewhat COMPAS aware.


Usage
=====

h5view.py [-h] [-f FILENAME_FILTER] [-r [RECURSION_DEPTH]] [-S] [-H]
                 [-C [CONTENTS]] [-s] [-x EXCLUDE_GROUP [EXCLUDE_GROUP ...]]
                 [-V SEED_LIST [SEED_LIST ...]]
                 input [input ...]

HDF5 file content viewer.

positional arguments:
  input
    input directory and/or file name(s)

optional arguments:
  -h, --help
    show this help message and exit
  -f FILENAME_FILTER, --filter FILENAME_FILTER
    input filename filter (default = *)
  -r [RECURSION_DEPTH], --recursive [RECURSION_DEPTH]
    recursion depth (default is no recursion)
  -S, --summary
    display summary output for HDF5 file (default is not to displat summary)
  -H, --headers
    display file headers for HDF5 file (default is not to display headers)
  -C [CONTENTS], --contents [CONTENTS]
    display file contents for HDF5 file: argument is number of entries (+ve from top, -ve
    from bottom) (default is not to display contents)
  -s, --stop-on-error
    stop all copying if an error occurs (default is skip to next file and continue)
  -x EXCLUDE_GROUP [EXCLUDE_GROUP ...], --exclude EXCLUDE_GROUP [EXCLUDE_GROUP ...]
    list of input groups to be excluded (default is all groups will be copied)
  -V SEED_LIST [SEED_LIST ...], --seeds SEED_LIST [SEED_LIST ...]
    list of seeds to be printed (for content printing) (default is print all seeds)



Functionality overview
======================

Displays summary, header and content information for specified COMPAS HDF5 file(s).  If none
of --summary [-S], --headers [-H], or --contents [-C] are specified, --summary is assumed.
If any of --summary [-S], --headers [-H], or --contents [-C] are specified, then only the 
option(s) specified are actioned.


Summary information displays, for each COMPAS file in the HDF5 file:
   - the name of the COMPAS file
   - the number of columns in the COMPAS file
   - the number of entries in the COMPAS file (actually, the maximum number of entries in any column in the COMPAS file)


Header information displays, for each COMPAS file in the HDF5 file:
   - the name of each column in the COMPAS file
   - the number of entries in each column of the COMPAS file
   - the data type of each column of the COMPAS file
   - the units associated with each column of the COMPAS file
     (with the exception of the Run_Details file - there are no units associated with columns in the Run_Details file)


Contents information displays, for each COMPAS file in the HDF5 file:
   - a header showing the column names in the COMPAS file
   - a row for each entry in the COMPAS file, showing the column values for that row (comma delimited)

   The contents display can be limited in two ways:

      (a) The --contents option takes and optional argument: an integer number of rows to display.
          The argement to --contents can be positive or negative: a positive value indicates that
          the number of rows specified by the argument should be displayed from the start of the file;
          a negative value indicates that the number of fows specified by the (absolute value of the)
          argument should be displayed from the end of the file.  The +ve and -ve arguments to the
          --contents option are akin the the Unix 'head' and 'tail' commands.

      (b) The --seeds option allows the user to specify a list of SEED values that should be printed.
          If the --seeds option is specified, only rows containing the seeds specified by the user will
          be printed - and only if the are in the entries printed if limited by the --contents argument
          described in (a).

          Note that printing only seeds specified in a list of seeds could be slow - we effectively have
          to look through the entire dataset looking for the seeds required...


JR, April 2021
'''

#!/usr/bin/python3
import sys
import os
import math
import datetime
import numpy as np
import h5py as h5
import argparse
from fnmatch import fnmatch

# h5view assumes the COMPAS run details filename is 'Run_Details'
# change the filename on the next line if necessary
Run_Details_Filename = 'Run_Details'


# getDataType()
#
# Determine COMPAS datatype from hdf5 datatype

def getDataType(dType = ''):

    dataType = 'UNKNOWN'
    
    if   dType == 'uint8'  : dataType = 'BOOL'
    elif dType == 'uint16' : dataType = 'UNSIGNED_SHORT_INT'
    elif dType == 'uint32' : dataType = 'UNSIGNED_INT'
    elif dType == 'uint64' : dataType = 'UNSIGNED_LONG_INT'
    elif dType == 'int16'  : dataType = 'SHORT_INT'
    elif dType == 'int32'  : dataType = 'INT'
    elif dType == 'int64'  : dataType = 'LONG_INT'
    elif dType == 'float32': dataType = 'FLOAT' 
    elif dType == 'float64': dataType = 'DOUBLE'
    elif dType[0:2] == '|S': 
        strLen = int(dType[2:len(dType)]) - 1       # strip null terminator
        dataType = 'STRING(' + str(strLen) + ')'

    return dataType


# printSummary()
#
# Prints HDF5 file summary
#
# For each COMPAS file (group) in the HDF5 file, prints:
#
#    - COMPAS filename
#    - Number of columns in the file (datsets in the group)
#    - Maximum number of entries of any column (entries per column should be equal)

def printSummary(h5name = None, h5file = None, excludeList = ''):

    ok = True                                                                                       # result
  
    try:

        mtime = os.path.getmtime(h5name)                                                            # time file last modified
        lastModified = datetime.datetime.fromtimestamp(mtime)                                       # ... formatted

        fileSize = os.path.getsize(h5name)                                                          # file size (in bytes)
        strFileSize = ('{:<11.4f}').format(fileSize / 1024.0 / 1024.0 / 1024.0)                     # ... formatted in GB

        print('\n\nSummary of HDF5 file', h5name)
        print('='*(21 + len(h5name)))

        print('\nFile size    :', strFileSize.strip(), 'GB')
        print('Last modified:', lastModified)

        # get widths for columns to be displayed - it's a bit more work, 
        # and a bit redundant, but it is neater... (and we don't have thousands of files...)
        maxFilenameLen  = -1
        maxColumns      = -1
        groupMaxEntries = []

        keyList = list(h5file.keys())
        detailedOutput = isinstance(h5file[keyList[0]], h5.Dataset)                                 # processing detailed output file?

        if detailedOutput:                                                                          # detailed output file

            maxColumns = len(keyList)
            groupMaxEntries.append(-1)
            for dIdx, dataset in enumerate(h5file.keys()):
                if h5file[dataset].shape[0] > groupMaxEntries[0]: groupMaxEntries[0] = h5file[dataset].shape[0]

        else:                                                                                       # assume COMPAS_Output file

            for gIdx, group in enumerate(h5file.keys()):

                groupMaxEntries.append(-1)

                if len(group) > maxFilenameLen: maxFilenameLen = len(group)
                if len(h5file[group].keys()) > maxColumns: maxColumns = len(h5file[group].keys())

                columns = h5file[group].keys()
                for idx, column in enumerate(columns):
                    if h5file[group][column].shape[0] > groupMaxEntries[gIdx]: groupMaxEntries[gIdx] = h5file[group][column].shape[0]



        for widthColumns in range(10):                                                              # better not be more than 10**10 columns!
            if maxColumns < 10**widthColumns: break                                                 # 'columns' width

        maxEntries = max(groupMaxEntries)
        for widthEntries in range(10):                                                              # better not be more than 10**10 entries per column!
            if maxEntries < 10**widthEntries: break                                                 # 'entries' width

        if detailedOutput:                                                                          # detailed output file

            print(('\n{:<' + str(max(7, widthColumns)) + '}   {:<' + str(max(7, widthEntries)) + '}')
                  .format('Columns', 'Entries'))
            print('-'*(max(7, widthColumns)), ' ', '-'*(max(7, widthEntries)))

            print(('{:>' + str(max(7, widthColumns)) + '}   {:>' + str(max(7, widthEntries)) + '}')
                  .format(len(keyList), groupMaxEntries[0]))

        else:                                                                                       # assume COMPAS_Output file

            print(('\n{:<' + str(maxFilenameLen) + '}   {:<' + str(max(7, widthColumns)) + '}   {:<' + str(max(7, widthEntries)) + '}   {:<' + str(max(12, widthEntries)) + '}')
                  .format('COMPAS Filename', 'Columns', 'Entries', 'Unique SEEDs'))
            print('-'*(maxFilenameLen), ' ', '-'*(max(7, widthColumns)), ' ', '-'*(max(7, widthEntries)), ' ', '-'*(max(12, widthEntries)))

            # do Run_Details file first
            if not Run_Details_Filename in excludeList:                                             # ... if not excluded
                print(('{:<' + str(maxFilenameLen) + '}   {:>' + str(max(7, widthColumns)) + '}   {:>' + str(max(7, widthEntries)) + '}')
                      .format(Run_Details_Filename, len(h5file[Run_Details_Filename].keys()), len(h5file[Run_Details_Filename][list(h5file[Run_Details_Filename].keys())[0]])))

            # do remaining files (groups)
            for gIdx, group in enumerate(h5file.keys()):
                if group in excludeList: continue                                                   # skip if excluded
                if group == Run_Details_Filename: continue                                          # Run_details already done (or excluded)

                try:
                    uniqueSeedsStr = str(len(np.unique(h5file[group]['SEED'])))
                except Exception as e:
                    uniqueSeedsStr = " "

                print(('{:<' + str(maxFilenameLen) + '}   {:>' + str(max(7, widthColumns)) + '}   {:>' + str(max(7, widthEntries)) + '}   {:>' + str(max(12, widthEntries)) + '}')
                      .format(group, len(h5file[group].keys()), groupMaxEntries[gIdx], uniqueSeedsStr))

        print('\n')

    except Exception as e:                                                                          # error occurred accessing the input file
        print('printSummary: Error accessing HDF5 file', h5name, ':', str(e))
        ok = False

    return ok


# printHeaders()
#
# Prints headers for each COMPAS file (group) in the HDF5 file
#
# For each COMPAS file (group) in the HDF5 file:
#
#    - for each column (dataset) in the COMPAS file (group), prints:
#       - column (dataset) name
#       - actual number of entries in the column (dataset)
#       - the COMPAS data type of the column (dataset)
#       - the units for the column (dataset) (except for the Run_Deatils file - 
#         columns (datasets) in the Run-Details file have no units)

def printHeaders(h5name = None, h5file = None, excludeList = ''):

    ok = True                                                                                       # result
  
    try:

        print('\n\nHeaders for HDF5 file', h5name)
        print('='*(22 + len(h5name)), '\n')

        keyList = list(h5file.keys())
        detailedOutput = isinstance(h5file[keyList[0]], h5.Dataset)                                 # processing detailed output file?

        if detailedOutput:

            # get widths for columns to be displayed - it's a bit more
            # work, and a bit redundant, but it is neater... (and we don't have thousands of columns...)
            maxDatasetNameLen = -1
            maxEntries        = -1
            maxDataTypeLen    = -1
            maxUnitsLen       = -1
            for column in h5file.keys():
                if len(column) > maxDatasetNameLen: maxDatasetNameLen = len(column)
                if h5file[column].shape[0] > maxEntries: maxEntries = h5file[column].shape[0]

                dataType = getDataType(str(h5file[column].dtype))
                if len(dataType) > maxDataTypeLen: maxDataTypeLen = len(dataType)

                # get units - not for run details file
                try:
                    units = h5file[column].attrs['units'].decode('utf-8')
                except AttributeError:
                    units = h5file[column].attrs['units']

                if len(units) > maxUnitsLen: maxUnitsLen = len(units)

            for widthEntries in range(10):                                                          # better not be more than 10**10 entries per column!
                if maxEntries < 10**widthEntries: break                                             # 'entries' width

            # print data

            print(('\n{:<' + str(maxDatasetNameLen) + '}   {:<' + str(max(7, widthEntries)) + '}   {:<' + str(max(8, maxDataTypeLen)) + '}   {:<' + str(max(5, maxUnitsLen)) + '}')
                  .format('Datatset Name', 'Entries', 'Data Type', 'Units'))
            print('.'*(maxDatasetNameLen), ' ', '.'*(max(7, widthEntries)), ' ', '.'*(max(9, maxDataTypeLen)), ' ', '.'*(max(5, maxUnitsLen)))

            for column in h5file.keys():
                dataType = getDataType(str(h5file[column].dtype))

                 # get units
                try:
                    units = h5file[column].attrs['units'].decode('utf-8')
                except AttributeError:
                    units = h5file[column].attrs['units']

                print(('{:<' + str(maxDatasetNameLen) + '}   {:>' + str(max(7, widthEntries)) + '}   {:<' + str(max(9, maxDataTypeLen)) + '}   {:<' + str(max(5, maxUnitsLen)) + '}')
                      .format(column, h5file[column].shape[0], dataType, units))

            print()

        else:

            # do Run_Details file first
            if not Run_Details_Filename in excludeList:                                             # ... if not excluded

                # get widths for columns to be displayed - it's a bit more
                # work, and a bit redundant, but it is neater... (and we don't have thousands of columns...)
                maxDatasetNameLen = -1
                maxEntries        = -1
                maxDataTypeLen    = -1
                for column in h5file[Run_Details_Filename].keys():
                    if len(column) > maxDatasetNameLen: maxDatasetNameLen = len(column)
                    if h5file[Run_Details_Filename][column].shape[0] > maxEntries: maxEntries = h5file[Run_Details_Filename][column].shape[0]

                    dataType = getDataType(str(h5file[Run_Details_Filename][column].dtype))
                    if len(dataType) > maxDataTypeLen: maxDataTypeLen = len(dataType)

                for widthEntries in range(10):                                                      # better not be more than 10**10 entries per column!
                    if maxEntries < 10**widthEntries: break                                         # 'entries' width

                # print data
                print('COMPAS file:', Run_Details_Filename)
                print('-'*(13 + len(Run_Details_Filename)))

                print(('\n{:<' + str(maxDatasetNameLen) + '}   {:<' + str(max(7, widthEntries)) + '}   {:<' + str(max(8, maxDataTypeLen)) + '}')
                      .format('Datatset Name', 'Entries', 'Data Type'))
                print('.'*(maxDatasetNameLen), ' ', '.'*(max(7, widthEntries)), ' ', '.'*(max(9, maxDataTypeLen)))

                for column in h5file[Run_Details_Filename].keys():
                    dataType = getDataType(str(h5file[Run_Details_Filename][column].dtype))

                    print(('{:<' + str(maxDatasetNameLen) + '}   {:>' + str(max(7, widthEntries)) + '}   {:<' + str(max(8, maxDataTypeLen)) + '}')
                          .format(column, h5file[Run_Details_Filename][column].shape[0], dataType))

                print()

            # do remaining files (groups)
            for group in h5file.keys():
                if group in excludeList: continue                                                   # skip if excludedd
                if group == Run_Details_Filename: continue                                          # Run_details already done (or excluded)

                if not Run_Details_Filename in excludeList: print()

                # get widths for columns to be displayed - it's a bit more
                # work, and a bit redundant, but it is neater... (and we don't have thousands of columns...)
                maxDatasetNameLen = -1
                maxEntries        = -1
                maxDataTypeLen    = -1
                maxUnitsLen       = -1
                for column in h5file[group].keys():
                    if len(column) > maxDatasetNameLen: maxDatasetNameLen = len(column)
                    if h5file[group][column].shape[0] > maxEntries: maxEntries = h5file[group][column].shape[0]

                    dataType = getDataType(str(h5file[group][column].dtype))
                    if len(dataType) > maxDataTypeLen: maxDataTypeLen = len(dataType)

                    # get units - not for run details file
                    try:
                        units = h5file[group][column].attrs['units'].decode('utf-8')
                    except AttributeError:
                        units = h5file[group][column].attrs['units']

                    if len(units) > maxUnitsLen: maxUnitsLen = len(units)

                for widthEntries in range(10):                                                      # better not be more than 10**10 entries per column!
                    if maxEntries < 10**widthEntries: break                                         # 'entries' width

                # print data
                print('COMPAS file:', group)
                print('-'*(13 + len(group)))

                print(('\n{:<' + str(maxDatasetNameLen) + '}   {:<' + str(max(7, widthEntries)) + '}   {:<' + str(max(8, maxDataTypeLen)) + '}   {:<' + str(max(5, maxUnitsLen)) + '}')
                      .format('Datatset Name', 'Entries', 'Data Type', 'Units'))
                print('.'*(maxDatasetNameLen), ' ', '.'*(max(7, widthEntries)), ' ', '.'*(max(9, maxDataTypeLen)), ' ', '.'*(max(5, maxUnitsLen)))

                for column in h5file[group].keys():
                    dataType = getDataType(str(h5file[group][column].dtype))

                    # get units
                    try:
                        units = h5file[group][column].attrs['units'].decode('utf-8')
                    except AttributeError:
                        units = h5file[group][column].attrs['units']

                    print(('{:<' + str(maxDatasetNameLen) + '}   {:>' + str(max(7, widthEntries)) + '}   {:<' + str(max(9, maxDataTypeLen)) + '}   {:<' + str(max(5, maxUnitsLen)) + '}')
                          .format(column, h5file[group][column].shape[0], dataType, units))

                print()

        print('\n')

    except Exception as e:                                                                          # error occurred accessing the input file
        print('printHeaders: Error accessing HDF5 file', h5name, ':', str(e))
        ok = False

    return ok


# printContents()
#
# Prints the contents of each COMPAS file (group) in the HDF5 file
#
# Prints the values of the entries in each column for each COMPAS file (group) in the HDF5 file
#
# What is printed can be limited by the 'count' and 'seed' parameters.  If count is specified
# it is the number of entries to print.  A +ve count indicates the entries printed are from the
# head (start) of the file; a -ve count indicates the entries are from the tail (end) of the file
# (akin to the Unix 'head' and 'tail' commands)
#
# If a list of seeds is specified via the seed parameter, only entries matching those seed values
# will be printed - and only if the are in the entries printed if limited by count.
#
# Printing only seeds specified in a list of seeds could be slow - we effectively have to look
# through the entire dataset looking for the seeds required...

def printContents(h5name = None, h5file = None, excludeList = '', count = sys.maxsize, seeds = []):

    ok = True                                                                                   # result
  
    if count == 0: return ok                                                                    # nothing to do

    try:

        print('\n\nContents of HDF5 file', h5name)
        print('='*(22 + len(h5name)), '\n')

        keyList = list(h5file.keys())
        detailedOutput = isinstance(h5file[keyList[0]], h5.Dataset)                                 # processing detailed output file?

        if detailedOutput:


            columns = list(h5file.keys())

            hdr = ''
            for column in columns: 
                if hdr == '': hdr += column                                                     # add value to hdr (first value)
                else        : hdr += ', ' + column                                              # add value to hdr (subsequent values)

            print(hdr)  
            
            rows = len(h5file[columns[0]])
            rowsPrinted = 0
            for idx in range(rows):           

                printIt = True

                index = rows - idx - 1 if count < 0 else idx                                    # correct index - head or tail of file       

                row = ''
                for column in h5file.keys():
                    try:
                        value = h5file[column][index].decode('utf-8')
                    except AttributeError:
                        value = h5file[column][index]

                    if column == 'Run-End' or column == 'Run-Start': value = value[:-1]         # strip trailing \n

                    if seeds and column == 'SEED' and not(value in seeds):
                        printIt = False
                        break

                    dataType = getDataType(str(h5file[column].dtype))                    # data type
                    if dataType[0:6] == 'STRING':                                               # string?
                        if value[0:1] == "'": value = value[1:]                                 # yes - strip leading quote if present
                        if value[len(value)-1:len(value)] == "'": value = value[:-1]            # strip trailing quote if present
                        value = value.strip()                                                   # strip leading and trailing blanks 

                    if row == '': row += str(value)                                             # add value to output row (first value)
                    else        : row += ', ' + str(value)                                      # add value to outputrow (subsequent values)

                if printIt:
                    print(row)

                    if count < sys.maxsize:                                                     # count limited?
                        rowsPrinted += 1                                                        # yes - increment number printed
                        if rowsPrinted >= abs(count): break                                     # check for done

            print()




        else:

            # do Run_Details file first
            if not Run_Details_Filename in excludeList:                                             # ... if not excluded

                # print data
                print('COMPAS file:', Run_Details_Filename)
                print('-'*(13 + len(Run_Details_Filename)))

                columns = list(h5file[Run_Details_Filename].keys())

                hdr = ''
                for column in columns: 
                    if hdr == '': hdr += column                                                     # add value to hdr (first value)
                    else        : hdr += ', ' + column                                              # add value to hdr (subsequent values)

                print(hdr)  
            
                rows = len(h5file[Run_Details_Filename][columns[0]])
                rowsPrinted = 0
                for idx in range(rows):

                    index = rows - idx - 1 if count < 0 else idx                                    # correct index - head or tail of file       

                    row = ''
                    for column in h5file[Run_Details_Filename].keys():
                        try:
                            value = h5file[Run_Details_Filename][column][index].decode('utf-8')
                        except AttributeError:
                            value = h5file[Run_Details_Filename][column][index]

                        if column == 'Run-End' or column == 'Run-Start': value = value[:-1]         # strip trailing \n

                        dataType = getDataType(str(h5file[Run_Details_Filename][column].dtype))     # data type
                        if dataType[0:6] == 'STRING':                                               # string?
                            if value[0:1] == "'": value = value[1:]                                 # yes - strip leading quote if present
                            if value[len(value)-1:len(value)] == "'": value = value[:-1]            # strip trailing quote if present
                            value = value.strip()                                                   # strip leading and trailing blanks 

                        if row == '': row += str(value)                                             # add value to output row (first value)
                        else        : row += ', ' + str(value)                                      # add value to outputrow (subsequent values)

                    print(row)

                    if count < sys.maxsize:                                                         # count limited?
                        rowsPrinted += 1                                                            # yes - increment number printed
                        if rowsPrinted >= abs(count): break                                         # check for done

                print()

        
            # do remaining files (groups)
            for group in h5file.keys():
                if group in excludeList: continue                                                   # skip if excludedd
                if group == Run_Details_Filename: continue                                          # Run_details already done (or excluded)

                if not Run_Details_Filename in excludeList: print()

                # print data
                print('COMPAS file:', group)
                print('-'*(13 + len(group)))

                columns = list(h5file[group].keys())

                hdr = ''
                for column in columns: 
                    if hdr == '': hdr += column                                                     # add value to hdr (first value)
                    else        : hdr += ', ' + column                                              # add value to hdr (subsequent values)

                print(hdr)  
            
                rows = len(h5file[group][columns[0]])
                rowsPrinted = 0
                for idx in range(rows):           

                    printIt = True

                    index = rows - idx - 1 if count < 0 else idx                                    # correct index - head or tail of file       

                    row = ''
                    for column in h5file[group].keys():
                        try:
                            value = h5file[group][column][index].decode('utf-8')
                        except AttributeError:
                            value = h5file[group][column][index]

                        if column == 'Run-End' or column == 'Run-Start': value = value[:-1]         # strip trailing \n

                        if seeds and column == 'SEED' and not(value in seeds):
                            printIt = False
                            break

                        dataType = getDataType(str(h5file[group][column].dtype))                    # data type
                        if dataType[0:6] == 'STRING':                                               # string?
                            if value[0:1] == "'": value = value[1:]                                 # yes - strip leading quote if present
                            if value[len(value)-1:len(value)] == "'": value = value[:-1]            # strip trailing quote if present
                            value = value.strip()                                                   # strip leading and trailing blanks 

                        if row == '': row += str(value)                                             # add value to output row (first value)
                        else        : row += ', ' + str(value)                                      # add value to outputrow (subsequent values)

                    if printIt:
                        print(row)

                        if count < sys.maxsize:                                                     # count limited?
                            rowsPrinted += 1                                                        # yes - increment number printed
                            if rowsPrinted >= abs(count): break                                     # check for done

                print()

        print('\n')
        
    except Exception as e:                                                                      # error occurred accessing the input file
        print('printContents: Error accessing HDF5 file', h5name, ':', str(e))
        ok = False

    return ok


# viewHDF5File()
#
# Displays the contents of the file passed in 'path' parameter
# The amount of data displayed is governed by the parameters:
#
# - summary
# - headers
# - count
#
# See relevant functions above for detail

def viewHDF5File(path, excludeList = '', summary = True, headers = False, count = sys.maxsize, seeds = []):

    ok = True                                                                                                                       # result
  
    try:
        with h5.File(path, 'r') as srcFile:                                                                                         # open the input HDF5 file

            srcFname = os.path.abspath(srcFile.filename)                                                                            # fully-qualified source filename

            if summary:                                                                                                             # summary?
                ok = printSummary(h5name = srcFname, h5file = srcFile, excludeList = excludeList)                                   # yes - print summary

            if ok and headers:                                                                                                      # headers?
                ok = printHeaders(h5name = srcFname, h5file = srcFile, excludeList = excludeList)                                   # yes - print headers

            if ok and count != 0:                                                                                                   # contents?
                ok = printContents(h5name = srcFname, h5file = srcFile, excludeList = excludeList, count = count, seeds = seeds)    # yes - print contents


    except Exception as e:                                                                                                          # error occurrd accessing the input file
        print('Error accessing HDF5 file', path, ':', str(e))
        ok = False

    return ok


# processDirectory()
#
# Processes file and directories in directory specified by 'path' parameter
# Recursion is controlled by 'recursive' and 'depth' parameters

def processDirectory(path, 
                     recursive   = 0, 
                     fileFilter  = '*.h5', 
                     stopOnError = False, 
                     depth       = 0, 
                     excludeList = '',
                     summary     = True, 
                     headers     = False, 
                     count       = sys.maxsize,
                     seeds       = []):

    ok = True                                                           # result

    thisPath = os.path.abspath(path)                                    # absolute path

    try:
        for dirpath, dirnames, filenames in os.walk(thisPath):          # walk directory
            absDirpath = os.path.abspath(dirpath)                       # absolute path
            print('Processing directory', absDirpath)                   # announce directory being processed

            for filename in filenames:                                  # for each filename
                if fnmatch(filename, fileFilter):                       # filename matches filter?
                    ok = viewHDF5File(absDirpath + '/' + filename, 
                                      excludeList = excludeList,
                                      summary     = summary,
                                      headers     = headers,
                                      count       = count,
                                      seeds       = seeds)              # yes - view the file

                if stopOnError and not ok: break                        # check error - stop if required

            if stopOnError and not ok: break                            # check error - stop if required

            if depth < recursive:                                       # recurse?
                depth += 1                                              # increment recursion depth

                for dirname in dirnames:                                # for each directory
                    processDirectory(absDirpath + '/' + dirname, 
                                     recursive   = recursive, 
                                     fileFilter  = fileFilter, 
                                     stopOnError = stopOnError, 
                                     depth       = depth, 
                                     excludeList = excludeList,
                                     summary     = summary,
                                     headers     = headers,
                                     count       = count,
                                     seeds       = seeds)               # process the directory

            break                                                       # control recursion

    except Exception as e:                                              # error occurred accessing directory
        print('Error accessing directory', thisPath, ':', str(e))       # announce error
        ok = False                                                      # we're done

    return ok


def main():

    # setup argument parser
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position = 4, width = 90)
    parser = argparse.ArgumentParser(description = 'HDF5 file content viewer.', formatter_class = formatter)

    # define arguments
    parser.add_argument('inputPaths', metavar = 'input', type = str, nargs = '+', help = 'input directory and/or file name(s)')
    parser.add_argument('-f', '--filter', dest = 'filename_filter', type = str, action = 'store',  default = '*', help = 'input filename filter (default = *)')
    parser.add_argument('-r', '--recursive', dest = 'recursion_depth', type = int, nargs = '?', action = 'store', default = 0, const = sys.maxsize,  help = 'recursion depth (default is no recursion)')
    parser.add_argument('-S', '--summary', dest = 'summary', action = 'store_true',  default = False, help = 'display summary output for HDF5 file (default is not to display summary)')
    parser.add_argument('-H', '--headers', dest = 'headers', action = 'store_true',  default = False, help = 'display file headers for HDF5 file (default is not to display headers)')
    parser.add_argument('-C', '--contents', dest = 'contents', type = int, nargs = '?', action = 'store', default = 0, const = sys.maxsize, help = 'display file contents for HDF5 file: argument is number of entries (+ve from top, -ve from bottom) (default is not to display contents)')
    parser.add_argument('-s', '--stop-on-error', dest = 'stop_on_error', action = 'store_true',  default = False, help = 'stop all copying if an error occurs (default is skip to next file and continue)')
    parser.add_argument('-x', '--exclude', dest = 'exclude_group', type = str, nargs = '+', action = 'store', default = '', help = 'list of input groups to be excluded (default is all groups will be copied)')
    parser.add_argument('-V', '--seeds', dest = 'seed_list', type = int, nargs = '+', action = 'store', default = [], help = 'list of seeds to be printed (for content printing) (default is print all seeds)')

    # parse arguments
    args = parser.parse_args()

    fileFilter = args.filename_filter + '.h5'                                                               # add file extension to filter

    excludeList = ' '.join(args.exclude_group)                                                              # construct exclude list for groups

    if not args.summary and not args.headers and not args.contents: args.summary = True                     # summary is default if nothing else specified

    # process input files and directories
    for thisPath in args.inputPaths:                                                                        # for each input path
        thisFullPath = os.path.abspath(thisPath)                                                            # fully-qualified filename
        if os.path.exists(thisFullPath):                                                                    # path exists?
            if os.path.isfile(thisFullPath):                                                                # yes - is it a file?
                if fnmatch(thisPath, fileFilter):                                                           # yes - filename matches filter?
                    ok = viewHDF5File(thisFullPath, 
                                      excludeList = excludeList,
                                      summary     = args.summary, 
                                      headers     = args.headers, 
                                      count       = args.contents,
                                      seeds       = args.seed_list)                                         # yes - view it
                else:                                                                                       # no - does not match filter
                    print('Warning:', thisPath, 'does not match file filter (', fileFilter, '): ignored')   # show warning
            elif os.path.isdir(thisFullPath):                                                               # not a file - directory?
                                                                                                            # yes - process directory
                ok = processDirectory(thisFullPath, 
                                      recursive   = args.recursion_depth, 
                                      fileFilter  = fileFilter, 
                                      stopOnError = args.stop_on_error, 
                                      excludeList = excludeList,
                                      summary     = args.summary, 
                                      headers     = args.headers, 
                                      count       = args.contents,
                                      seeds       = args.seed_list)
            else:                                                                                           # not a file or directory
                print('Warning:', thisFullPath, 'is not a file or a directory: ignored')                    # show warning
        else:                                                                                               # path does not exist
            print('Warning:', thisFullPath, 'does not exist: ignored')                                      # show warning

        if args.stop_on_error and not ok:                                                                   # error occurred, and stop-on-error specified?
            print('Error encountered: view stopped')                                                        # yes - announce error
            break                                                                                           # and stop


if __name__ == "__main__":
    main()