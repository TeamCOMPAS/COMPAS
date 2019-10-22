import sys
import numpy as np
import argparse


# process file - read, parse and return data
#
# reads the file passed as parameter 'filename'
# parses record using delimiter passed as parameter 'delimiter'
#
# returns lists:
#
#    names  - 1-d list of column names (headings)
#    types  - 1-d list of column data types
#    units  - 1-d list of column units
#    values - 2-d list of data values: values[column][row]
#
#    all returned lists will be empty on failure (e.g. unable to open file)

def ProcessFile(filename, delimiter):

    names  = []
    types  = []
    units  = []
    values = []

    # read and parse the data from the file

    try:
        fi = open(filename)

        types = fi.readline().split(delimiter)                      # read types header record and split on delimiter
        types = [t.strip() for t in types]                          # strip leading and trailing whitespace (including \n) from each element in list

        units = fi.readline().split(delimiter)                      # read units header record and split on delimiter
        units = [u.strip() for u in units]                          # strip leading and trailing whitespace (including \n) from each element in list

        # rewind file to beginning, read and parse the data from the file
        # data values are extracted to 2-d list ([column][row])

        fi.seek(0)
        data = np.genfromtxt(fi, skip_header = 2, delimiter = delimiter, deletechars = " !#$%&'*+,.:;=?[\]{|}~", names = True)

        try:
            fi.close()

            names = data.dtype.names                                                                # record column names

            values = [[0.0 for y in range(len(data))] for x in range(len(data.dtype.names))]        # size the 2-d values list

            # read the data from the file into 2-d values list
            for row in range(len(data)):

                # extract the data values to a 2-d list
                for column in range(min(len(data.dtype.names), len(data[row]))):
                    values[column][row] = data[row][column]

        except (OSError, IOError):
            print('Close of file ', filename ,' failed')

    except (OSError, IOError):
        print('Unable to open file ', filename)

    return (names, types, units, values)


#
# Main program

# check we have one file to process

parser = argparse.ArgumentParser(description = 'COMPAS data file parse example.')

parser.add_argument('filename', metavar = 'filename', type = str, nargs = 1, help = 'name of COMPAS data file')

args = parser.parse_args()


# process file

print('\nProcessing file', args.filename[0], '\n')

names, types, units, values = ProcessFile(args.filename[0], ',')            # use this line for COMMA delimited files

#names, types, units, values = ProcessFile(args.filename[0], '\t')           # use this line for TAB delimited files


# check for failure

if len(values) < 2:
    print('Parse failed')
    sys.exit(1)


# print counts

print('column names  count =', len(names))
print('column types  count =', len(types))
print('values column count =', len(values))
print('values row    count =', len(values[0]))


# print column names

print('\nColumn names:\n')
column = 0
for n in names:
    print('Column', column, ', Name =', n)
    column += 1


# print column types

print('\nColumn types:\n')
column = 0
for t in types:
    print('Column', column, ', Type =', t)
    column += 1


# print column units

print('\nColumn units:\n')
column = 0
for u in units:
    print('Column', column, ', Units =', u)
    column += 1


# the data can be accessed by value[column][row]
# print a couple of values

if len(values) > 4 and len(values[0]) > 14:

    print()

    column = 3
    row    = 5
    print('Column', column, ', name =', names[column], ', type =', types[column], ', units =', units[column], ', value at row', row, '=', values[column][row])

    column = 4
    row    = 14
    print('Column', column, ', name =', names[column], ', type =', types[column], ', units =', units[column], ', value at row', row, '=', values[column][row])


# look for a specific column name

columnName = 'Radius_1'
try:
    whichColumn = names.index(columnName)
except ValueError:
    print('\nColumn name', columnName, 'not found!')
    whichColumn = -1

if whichColumn >= 0:
    print('\nColumn name', columnName, 'is column index', whichColumn, '(0-based)')

# print a couple of values from the column just found

if len(values[whichColumn]) > 17:

    print()

    column = whichColumn

    row = 3
    print('Column', column, ', name =', names[column], ', type =', types[column], ', units =', units[column], ', value at row', row, '=', values[column][row])

    row = 16
    print('Column', column, ', name =', names[column], ', type =', types[column], ', units =', units[column], ', value at row', row, '=', values[column][row])


# look for a specific column name that won't (well, shouldn't) be found

columnName = 'This_Should_Not_Be_Found'
try:
    whichColumn = names.index(columnName)
except ValueError:
    print('\nColumn name', columnName, 'not found!')
