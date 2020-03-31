###################################################################
#                                                                 #                                                               
#  Example of plotting detailed output COMPAS with python         #
#                                                                 #
###################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read file and create dataframe.
data_path = './COMPAS_Output/Detailed_Output/BSE_Detailed_Output_0.csv'
df = pd.read_table(data_path, sep=',', header=2, skipinitialspace=True)
df.columns = df.columns.str.strip() # Strip white space around headers

# Creating new column for total mass
df['Mtot']=df['Mass_1']+df['Mass_2']

# Create subplots arranged in a 2x3 table, sharing x-axes whenever plots have the same x-axis range.
fig, axes = plt.subplots(2,3,  sharex=True)

# Plot primary, secondary, and total mass
axes[0][0].plot(df['Time'], df['Mass_1'], linestyle='-', c='r', label='M1-Primary')
axes[0][0].plot(df['Time'], df['Mass_2'], linestyle='-', c='b', label='M2-Secondary')
axes[0][0].plot(df['Time'], df['Mtot'], linestyle='-', c='k', label='Mtot_C')
axes[0][0].set_ylabel('[Msun]')
axes[0][0].set_xlabel('time [Myr]')
axes[0][0].grid(linestyle=':', c='gray')
axes[0][0].legend(prop={'size':8})

# Plot fraction of Roche Lobe filled
axes[0][1].plot(df['Time'], df['Radius_1/RL'], linestyle='-', c='r', label='R1/RL1')
axes[0][1].plot(df['Time'], df['Radius_2/RL'], linestyle='-', c='b', label='R2/RL2')
axes[0][1].set_ylabel('[Rsun]')
axes[0][1].set_yscale('log')
axes[0][1].set_xlabel('time [Myr]')
axes[0][1].grid(linestyle=':', c='gray')
axes[0][1].legend(prop={'size':8})

# Plot radius
axes[0][2].plot(df['Time'], df['Radius_1'], linestyle='-', c='r', label='Radius1')
axes[0][2].plot(df['Time'], df['Radius_2'], linestyle='-', c='b', label='Radius2')
axes[0][2].set_ylabel('[Rsun]')
axes[0][2].set_xlabel('time [Myr]')
axes[0][2].grid(linestyle=':', c='gray')
axes[0][2].legend(prop={'size':8})

# Plot separation
axes[1][0].plot(df['Time'], df['Separation'], linestyle='-', c='k', label='Separation')
axes[1][0].set_ylabel('[Rsun]')
axes[1][0].set_xlabel('time [Myr]')
axes[1][0].legend(prop={'size':8})


# Plot eccentricity
axes[1][1].plot(df['Time'], df['Eccentricity'], linestyle='-', c='k', label= 'Eccentricity')
axes[1][1].set_ylabel('')
axes[1][1].set_xlabel('time [Myr]')
axes[1][1].grid(linestyle=':', c='gray')
axes[1][1].legend(prop={'size':8})

# Plot stellar types
axes[1][2].plot(df['Time'], df['Stellar_Type_1'], linestyle='-', c='r', label='StellarType1')
axes[1][2].plot(df['Time'], df['Stellar_Type_2'], linestyle='-', c='b', label='StellarType2')
axes[1][2].set_ylabel('Type')
axes[1][2].set_xlabel('time [Myr]')
axes[1][2].grid(linestyle=':', c='gray')
axes[1][2].legend(prop={'size':8})

# This is to replace the ticks in the last plot with stellar types with words 
types = ['MS<0','MS', 'HG', 'FGB', 'CHeB', 'EAGB', 'TPAGB', 'HeMS', 'HeHG',\
             'HeGB', 'HeWD', 'COWD', 'ONeWD', 'NS', 'BH', 'MR']
axes[1][2].set_yticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]) 
axes[1][2].set_yticklabels(types)
axes[1][2].tick_params(width=2, labelsize=10)#, labelleft='off')
axes[1][2].yaxis.grid(True)

fig.subplots_adjust(left=0.05, right=0.99) #adjusting boundaries of the plotter

# Below I automattically set the title to print the parameters of the first row in the fil
# -> df.iloc[0] and then the parameter of interest
fig.suptitle('M1=%s    M2=%s    a=%s Rsun   e=%s Z=0.02 \n  red= primary, blue = secondary, black=binary system' %(df.iloc[0]['Mass_1'], df.iloc[0]['Mass_2'], df.iloc[0]['Separation'], df.iloc[0]['Eccentricity']))

fig.subplots_adjust(wspace=.3)
plt.show()
