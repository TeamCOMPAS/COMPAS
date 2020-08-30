import numpy as np               #for handling arrays
import h5py as h5                #for reading the COMPAS data
import matplotlib.pyplot as plt  #for plotting
import math

# Get some data to plot
pathToData = 'COMPAS_Output.h5'
Data  = h5.File(pathToData)

seedOfInterest = 1

### Print all the parameter values next to their names

output_array = []

for strCategory in list(Data.keys()):       #['CommonEnvelopes', 'DoubleCompactObjects', 'Supernovae', 'SystemParameters']

    #output_col = []
    # Find occurences of the seed in this category
    category = Data[strCategory]
    allSeeds = category['SEED'][()]         # all the seeds which exist in this category
    maskSeeds = allSeeds == seedOfInterest     # boolean array which is true only where the seedOfInterest lives
    iSeeds = np.where(maskSeeds)[0]

    print("\n", strCategory, " has ", sum(maskSeeds), " line(s) which match")
    #output_col.append( "\n" + strCategory + " has " + str(sum(maskSeeds)) + " line(s) which match" )

    # If no occurences, skip this category
    if sum(maskSeeds) != 0: 
    
        for iSeed in iSeeds:
             for parameter in list(category.keys()):     # The subparams per category
                 print("\t" + parameter + " = " + str(category[parameter][()][iSeed])) 
                 #output_col.append("\t" + parameter + " = " + str(category[parameter][()][iSeed])) 
             print() 

    #output_array.append(output_col)

#print(output_array)

#toTranspose = list(zip_longest(*output_array))
#print(transpose(toTranspose))
