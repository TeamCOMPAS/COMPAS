###################################################################
#                                                                 #                                                               
#  Example of plotting detailed output COMPAS with python         #
#                                                                 #
# ##################################################################

import pandas as pd
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import matplotlib.ticker as ticker  
import os, sys

def main():

    ### Read file and create dataframe.
    data_path = './COMPAS_Output/Detailed_Output/BSE_Detailed_Output_0.h5'
    Data = h5.File(data_path, 'r')

    ### Creating new column for total mass
    Mtot = Data['Mass(1)'][()] + Data['Mass(2)'][()]
    
    ### Create subplots arranged in a 2x3 table, sharing x-axes whenever plots have the same x-axis range.
    fig, axes = plt.subplots(2,2, sharex=True, figsize=(10, 6)) # Dimensions chosen for aesthetics in the paper
    
    ### Plot mass attributes 
    lines = []  # This is used to split the legend into 2 pieces
    lines += axes[0][0].plot(Data['Time'][()], Mtot, linestyle='-', c='k', label='System Mass')
    lines += axes[0][0].plot(Data['Time'][()], Data['Mass(1)'][()], linestyle='-', c='r', label='Total Mass 1')
    lines += axes[0][0].plot(Data['Time'][()], Data['Mass_He_Core(1)'][()], linestyle='--', c='r', label='He Core 1')
    lines += axes[0][0].plot(Data['Time'][()], Data['Mass_CO_Core(1)'][()], linestyle=':', c='r', label='CO Core 1')
    lines += axes[0][0].plot(Data['Time'][()], Data['Mass(2)'][()], linestyle='-', c='b', label='Total Mass 2')
    lines += axes[0][0].plot(Data['Time'][()], Data['Mass_He_Core(2)'][()], linestyle='--', c='b', label='He Core 2')
    lines += axes[0][0].plot(Data['Time'][()], Data['Mass_CO_Core(2)'][()], linestyle=':', c='b', label='CO Core 2')
    axes[0][0].set_ylabel(r'Mass $/ \; M_{\odot}$')
    axes[0][0].grid(linestyle=':', c='gray')
    axes[0][0].legend(lines[:4], axes[0][0].get_legend_handles_labels()[1][:4], prop={'size':8}, loc=(0.05, 0.55))
    axes[0][0].add_artist(Legend(axes[0][0], lines[4:], axes[0][0].get_legend_handles_labels()[1][4:], prop={'size':8}, loc=(0.05, 0.15)))
    axes[0][0].tick_params(labelbottom=True)
    axes[0][0].xaxis.set_major_locator(ticker.MultipleLocator(2.0))
    axes[0][0].xaxis.set_minor_locator(ticker.MultipleLocator(1.0))
    axes[0][0].set_xlim(-0.5, 9)
    
    ### Plot radius attributes 
    axes[0][1].plot(Data['Time'][()], Data['SemiMajorAxis'][()], linestyle='-', c='k', label='Semi-Major Axis')
    axes[0][1].plot(Data['Time'][()], Data['Radius(1)'][()]/Data['Radius(1)|RL'][()], linestyle='--', c='r', label='Roche Radius 1')
    axes[0][1].plot(Data['Time'][()], Data['Radius(2)'][()]/Data['Radius(2)|RL'][()], linestyle='--', c='b', label='Roche Radius 2')
    axes[0][1].plot(Data['Time'][()], Data['Radius(1)'][()], linestyle='-', c='r', label='Stellar Radius 1')
    axes[0][1].plot(Data['Time'][()], Data['Radius(2)'][()], linestyle='-', c='b', label='Stellar Radius 2')
    axes[0][1].set_ylabel(r'Radius $/ \; R_{\odot}$')
    axes[0][1].set_yscale('log')
    axes[0][1].grid(linestyle=':', c='gray')
    axes[0][1].legend(prop={'size':8})
    axes[0][1].tick_params(labelbottom=True)
    
    ### Plot eccentricity
    axes[1][0].plot(Data['Time'][()], Data['Eccentricity'][()], linestyle='-', c='k') #, label= 'Eccentricity')
    axes[1][0].set_ylabel('Eccentricity')
    axes[1][0].set_xlabel('Time / Myr')
    axes[1][0].set_ylim(-0.05, 1.05)
    axes[1][0].grid(linestyle=':', c='gray')
    
    ### Plot stellar types
    stellarTypes, useTypes, typeNameMap = getStellarTypes(Data)
    
    axes[1][1].plot(Data['Time'][()], typeNameMap(Data['Stellar_Type(1)'][()]), linestyle='-', c='r', label='Stellar Type 1')
    axes[1][1].plot(Data['Time'][()], typeNameMap(Data['Stellar_Type(2)'][()]), linestyle='-', c='b', label='Stellar Type 2')
    axes[1][1].set_ylabel('Stellar Type')
    axes[1][1].set_xlabel('Time / Myr')
    axes[1][1].grid(linestyle=':', c='gray')
    axes[1][1].legend(prop={'size':8}) #, loc='lower left')
    axes[1][1].set_yticks(range(useTypes.shape[0]))
    axes[1][1].set_yticklabels([stellarTypes[typeNum] for typeNum in useTypes])
    axes[1][1].tick_params(width=2, labelsize=10)
    axes[1][1].yaxis.grid(True)
    
    ### Add in plot letters # The specific values were chosen for aesthetics in the COMPAS methods paper
    axes[0][0].text( axes[0][0].get_xlim()[1]*-.18, axes[0][0].get_ylim()[1]*0.87,  'a)', fontweight='bold')
    axes[0][1].text( axes[0][1].get_xlim()[1]*-.25, axes[0][1].get_ylim()[1]**0.75, 'b)', fontweight='bold') # exponent because it's on a log axis
    axes[1][0].text( axes[1][0].get_xlim()[1]*-.18, axes[1][0].get_ylim()[1]*0.87,  'c)', fontweight='bold')
    axes[1][1].text( axes[1][1].get_xlim()[1]*-.25, axes[1][1].get_ylim()[1]*0.87,  'd)', fontweight='bold')
    
    ### Print out a synopsys of the evolutionary history to the command line
    printEvolutionaryHistory(Data)
    
    ### Finalize the boundaries, save, and show
    fig.subplots_adjust(left=0.05)  #adjusting boundaries of the plotter
    fig.subplots_adjust(wspace=.3)
    plt.savefig('gw151226evol.eps', format='eps')
    plt.savefig('gw151226evol.png', format='png')
    plt.show()




### Helper functions

def getStellarTypes(Data):
    """
    This function extracts only the stellar types which actually arise in the binary's evolution,
    and produces a map between the used type numbers and names.
    """

    # List of Hurley stellar types
    stellarTypes = [r'MS$<0.7M_{\odot}$', r'MS$\geq0.7M_{\odot}$', 'HG', 'FGB', 'CHeB', 'EAGB', 'TPAGB', 'HeMS', 'HeHG', 'HeGB', 'HeWD', 'COWD', 'ONeWD', 'NS', 'BH', 'MR']

    useTypes = np.unique(np.append(Data['Stellar_Type(1)'][()], Data['Stellar_Type(2)'][()]))
    if (0 in useTypes) != (1 in useTypes): # XOR
        stellarTypes[0] = stellarTypes[1] = 'MS'

    def typeNameMap(x):
        return np.digitize(x, useTypes, right=True) 

    return stellarTypes, useTypes, typeNameMap


def printFormattedEvolutionLine(time, event, m1, t1, m2, t2, a, e):
    # All values are floats except event which is a string and t1, t2 which are ints
    print("{:10.6f}   {:21}  {:7.3f}    {:2}    {:7.3f}    {:2}   {:8.3f}  {:5.3f}" .format(time, event, m1, t1, m2, t2, a, e))

def printEvolutionaryHistory(Data):
    """
    This function prints a synopsys of the evolutionary history to the command line; it can eventually include cartoons as well.
    """

    print('Time (Myr), Event,                  M1 (M_o), type1, M2 (M_o), type2, a (R_o),   e')

    printFormattedEvolutionLine( Data['Time'][0], 'start:Z='+str(Data['Metallicity@ZAMS(1)'][0]),  
                                 Data['Mass(1)'][0], Data['Stellar_Type(1)'][0],  
                                 Data['Mass(2)'][0], Data['Stellar_Type(2)'][0],  
                                 Data['SemiMajorAxis'][0],  Data['Eccentricity'][0])

    for i in range(1,Data['Time'].size):
        if Data['MT_History'][i]>0:         #mass transfer happened
            if Data['MT_History'][i]==Data['MT_History'][i-1]:
                continue    #repeated entry
            if Data['MT_History'][i]==1:
                MTstring='Stable MT: 1 to 2'
            elif Data['MT_History'][i]==2:
                MTstring='Stable MT: 2 to 1'
            elif Data['MT_History'][i]==3:
                MTstring='CE: 1 to 2'
            elif Data['MT_History'][i]==4:
                MTstring='CE: 2 to 1'
            elif Data['MT_History'][i]==5:
                MTstring='CE: Double Core'
            elif Data['MT_History'][i]==6:
                MTstring='CE: both MS'
            elif Data['MT_History'][i]==7:
                MTstring='CE: MS with CO'
            elif Data['MT_History'][i]==8:
                MTstring='CE merger, rad env'
            elif Data['MT_History'][i]==9:
                MTstring='CE merger, >RLOF'
            else:
                MTstring='Unknown MT'

            printFormattedEvolutionLine( Data['Time'][i], MTstring, 
                                         Data['Mass(1)'][i], Data['Stellar_Type(1)'][i], 
                                         Data['Mass(2)'][i], Data['Stellar_Type(2)'][i], 
                                         Data['SemiMajorAxis'][i], Data['Eccentricity'][i])

        if Data['Stellar_Type(1)'][i]!=Data['Stellar_Type(1)'][i-1]:    #type of star 1 changed
            printFormattedEvolutionLine( Data['Time'][i], 'Star 1: '+str(Data['Stellar_Type(1)'][i-1])+'->'+str(Data['Stellar_Type(1)'][i]),
                                         Data['Mass(1)'][i], Data['Stellar_Type(1)'][i], 
                                         Data['Mass(2)'][i], Data['Stellar_Type(2)'][i],
                                         Data['SemiMajorAxis'][i], Data['Eccentricity'][i])
                        
        if Data['Stellar_Type(2)'][i]!=Data['Stellar_Type(2)'][i-1]:    #type of star 2 changed
            printFormattedEvolutionLine( Data['Time'][i], 'Star 2: '+str(Data['Stellar_Type(2)'][i-1])+'->'+str(Data['Stellar_Type(2)'][i]),
                                         Data['Mass(1)'][i], Data['Stellar_Type(1)'][i], 
                                         Data['Mass(2)'][i], Data['Stellar_Type(2)'][i],
                                         Data['SemiMajorAxis'][i], Data['Eccentricity'][i])

        if Data['Eccentricity'][i]>1 or Data['SemiMajorAxis'][i]<0:     #unbound
            print(Data['Time'][i], '   Unbound binary')

    isDCO = (Data['Stellar_Type(1)'][i]==13 or Data['Stellar_Type(1)'][i]==14) and (Data['Stellar_Type(2)'][i]==13 or Data['Stellar_Type(2)'][i]==14)
    isDoubleWD = (Data['Stellar_Type(1)'][i]==10 or Data['Stellar_Type(1)'][i]==11 or Data['Stellar_Type(1)'][i]==12) and (Data['Stellar_Type(2)'][i]==10 or Data['Stellar_Type(2)'][i]==11 or Data['Stellar_Type(2)'][i]==12)
    if Data['Time'][i]<14000 and not isDCO and not isDoubleWD and not (Data['Eccentricity'][i]>1 or Data['SemiMajorAxis'][i]<0):     #proxy for unrecorded meregr
        print(Data['Time'][i], '   Stellar merger')

    #Merger time calculation follows https://iopscience.iop.org/article/10.3847/2515-5172/ac2d35
    if ((Data['Stellar_Type(1)'][i]==13 or Data['Stellar_Type(1)'][i]==14) and (Data['Stellar_Type(2)'][i]==13 or Data['Stellar_Type(2)'][i]==14)):
        Msunkg=1.98892e30
        c=299792458
        G=6.67428e-11
        Rsun = 695500000
        beta=64/5*G**3*Data['Mass(1)'][i]*Data['Mass(2)'][i]*(Data['Mass(1)'][i]+Data['Mass(2)'][i])*Msunkg**3/c**5
        T0=(Data['SemiMajorAxis'][i]*Rsun)**4/4/beta
        e=Data['Eccentricity'][i]
        Tdelay=T0*(1-e**2)**(7/2)*(1+0.27*e**10 + 0.33*e**20 +  0.2*e**1000)/3.15e7/1e6
        printFormattedEvolutionLine( Data['Time'][i]+Tdelay, 'GW merger in {:.1f} Myr'.format(Tdelay),
                                     Data['Mass(1)'][i], Data['Stellar_Type(1)'][i],  
                                     Data['Mass(2)'][i], Data['Stellar_Type(2)'][i],  
                                     Data['SemiMajorAxis'][i],  Data['Eccentricity'][i])


if __name__ == "__main__":
    main()

