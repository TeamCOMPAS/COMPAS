###################################################################
#                                                                 #                                                               
#  Example of plotting detailed output COMPAS with python         #
#                                                                 #
# ##################################################################

import os, sys
import math
import numpy as np
import h5py as h5
import pandas as pd
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker  
from matplotlib.legend import Legend
from matplotlib import rcParams, transforms, patches
import matplotlib.gridspec as gridspec

def main():
    ### Read file and create dataframe.
    # data_path = '/home/rwillcox/astro/compas/COMPAS/output/detailed_evol_vanDenHeuval_plots/COMPAS_Output_2/Detailed_Output/BSE_Detailed_Output_0.h5'
    data_path = '/Users/13lauy1/git/compas2/COMPAS/postProcessing/BSE_Detailed_Output_0.h5'

    Data = h5.File(data_path, 'r')

    ### Collect the important events in the detailed evolution
    events = getAllEvents(Data) # Calculate the events here, for use in plot sizing parameters
    printEvolutionaryHistory(events=events)
    events = [event for event in events if event.eventClass != 'Stype'] # want to ignore simple stellar type changes

    ### Produce the two plots
    # makeDetailedPlots(Data, events)
    plotVanDenHeuval(events=events)
    plt.savefig('vanDenHeuvalPlot.eps', format='eps')
    plt.show()


fontparams = {
    "font.serif": "Times New Roman",
    "text.usetex": "True",
    "axes.grid": "True",
    "grid.color": "gray",
    "grid.linestyle": ":",
    "axes.titlesize": "18",
    "axes.labelsize": "12",
    "xtick.labelsize": "10",
    "ytick.labelsize": "10",
    "xtick.labelbottom": "True", 
    "legend.framealpha": "1",
    #"legend.fontsize": "None",
}



####### Functions to organize and call the plotting functions 

def makeDetailedPlots(Data=None, events=None):

    #listOfPlots = [ plotMassAttributes, plotLengthAttributes, plotStellarTypeAttributesAndEccentricity ]
    listOfPlots = [ plotMassAttributes, plotLengthAttributes, plotStellarTypeAttributes, plotEccentricity]

    events = [event for event in events if event.eventClass != 'Stype'] # want to ignore simple stellar type changes
    num_events = len(events)
    event_times = [event.time for event in events]


    rcParams.update(fontparams) # Set configurations for uniform plot output
    fig, axes = plt.subplots(nrows=len(listOfPlots), figsize=(10, 20)) # W, H
    #gs = gridspec.GridSpec(ncols=nEventsColumns, nrows=nRows, figure=fig, height_ratios=height_ratios)

    for ii, specificPlot in enumerate(listOfPlots): # exclude the last one

        #ax = fig.add_subplot(gs[ii,:])     
        ax = axes[ii]
        if ii == 0:
            ax0 = ax
        else:
            ax.sharex(ax0)

        # TODO: Set the reverse log scale for time

        # Plot the data
        handles, labels = specificPlot(fig=fig, ax=ax, Data=Data)

        # Add some breathing space at the top of the plot
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax*1.1)

        # Add vertical lines for specific event times
        [ax.axvline(time, ymin=0.975, zorder=0) for time in event_times]

        # On all plots, add the event letters
        # TODO: find a way to not let them overlap
        #if ii == 0:
        for jj in range(num_events):
            ax.text(x=event_times[jj], y=ax.get_ylim()[1]*1.05, s=chr(ord('@')+1+jj)) # The unicode representation of the capital letters - works as long as there are less than 26 images to show

        if handles is not None:
            ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1.03,0.5), fancybox=True)

        if (ii == len(listOfPlots)-1): # last of the regular plots
            ax.set_xlabel('Time / Myr')



    #### Finalize the boundaries, save, and show
    fig.suptitle('Detailed evolution for seed={}'.format(Data['SEED'][()][0]), fontsize=24) #, y=1)
    fig.tight_layout(h_pad=8, rect= (0, .08, 1, .95)) #, h_pad=8) # (left, bottom, right, top) 
    plt.savefig('detailedEvolutionPlot.eps', format='eps')
    #plt.show()    




######## Plotting functions


def plotMassAttributes(fig=None, ax=None, Data=None):

    ### Plot mass attributes 
    # Create new column for total mass
    Mtot = Data['Mass(1)'][()] + Data['Mass(2)'][()]
    ax.plot(Data['Time'][()], Mtot, linestyle='-', c='k', label='System Mass')
    ax.plot(Data['Time'][()], Data['Mass(1)'][()], linestyle='-', c='r', label='Total Mass 1')
    ax.plot(Data['Time'][()], Data['Mass_He_Core(1)'][()], linestyle='--', c='r', label='He Core 1')
    ax.plot(Data['Time'][()], Data['Mass_CO_Core(1)'][()], linestyle=':', c='r', label='CO Core 1')
    ax.plot(Data['Time'][()], Data['Mass(2)'][()], linestyle='-', c='b', label='Total Mass 2')
    ax.plot(Data['Time'][()], Data['Mass_He_Core(2)'][()], linestyle='--', c='b', label='He Core 2')
    ax.plot(Data['Time'][()], Data['Mass_CO_Core(2)'][()], linestyle=':', c='b', label='CO Core 2')

    ax.set_ylabel(r'Mass $/ \; M_{\odot}$')

    return ax.get_legend_handles_labels()
    
    
def plotLengthAttributes(fig=None, ax=None, Data=None):
          
    ### Plot radius attributes 
    ax.plot(Data['Time'][()], Data['SemiMajorAxis'][()], linestyle='-', c='k', label='Semi-Major Axis')
    ax.plot(Data['Time'][()], Data['Radius(1)'][()], linestyle='-', c='r', label='Stellar Radius 1')
    ax.plot(Data['Time'][()], Data['Radius(2)'][()], linestyle='-', c='b', label='Stellar Radius 2')
    # Need to mask out when the denominator is 0
    mask1 = Data['Radius(1)|RL'][()] != 0
    ax.plot(Data['Time'][()][mask1], Data['Radius(1)'][()][mask1]/Data['Radius(1)|RL'][()][mask1], linestyle='--', c='r', label='Roche Radius 1')
    mask2 = Data['Radius(2)|RL'][()] != 0
    ax.plot(Data['Time'][()][mask2], Data['Radius(2)'][()][mask2]/Data['Radius(2)|RL'][()][mask2], linestyle='--', c='b', label='Roche Radius 2')

    ax.set_ylabel(r'Radius $/ \; R_{\odot}$')
    ax.set_yscale('log')

    return ax.get_legend_handles_labels()
    

def plotEccentricity(fig=None, ax=None, Data=None):


    ### Plot eccentricity
    ax.plot(Data['Time'][()], Data['Eccentricity'][()], linestyle='-', c='k') #, label= 'Eccentricity')
    ax.set_ylabel('Eccentricity')

    ax.set_ylim(-0.05, 1.05)
    #ax.legend(framealpha=1, prop={'size':8} ) 
    ax.grid(linestyle=':', c='gray')
    
    #return ax.get_legend_handles_labels()
    return None, None
    
def plotStellarTypeAttributes(fig=None, ax=None, Data=None):


    ### Plot stellar types
    stellarTypes, useTypes, typeNameMap = getStellarTypes(Data)
    
    ax.plot(Data['Time'][()], typeNameMap(Data['Stellar_Type(1)'][()]), linestyle='-', c='r', label='Stellar Type 1')
    ax.plot(Data['Time'][()], typeNameMap(Data['Stellar_Type(2)'][()]), linestyle='-', c='b', label='Stellar Type 2')
    ax.set_ylabel('Stellar Type')
    ax.legend(prop={'size':8}) #, loc='lower left')
    ax.set_yticks(range(useTypes.shape[0]))
    ax.set_yticklabels([stellarTypes[typeNum] for typeNum in useTypes])
    #ax.tick_params(width=2, labelsize=10)
    ax.yaxis.grid(True)

    ax.legend(framealpha=1, prop={'size':8} ) 
    ax.grid(linestyle=':', c='gray')

    return ax.get_legend_handles_labels()
    

def plotStellarTypeAttributesAndEccentricity(fig=None, ax=None, Data=None):

    ax1 = ax
    ax2 = ax.twinx()
    
    ### Plot stellar types
    stellarTypes, useTypes, typeNameMap = getStellarTypes(Data)
    
    handle1 = ax1.plot(Data['Time'][()], typeNameMap(Data['Stellar_Type(1)'][()]), linestyle='-', c='r', label='Stellar Type 1')
    handle2 = ax1.plot(Data['Time'][()], typeNameMap(Data['Stellar_Type(2)'][()]), linestyle='-', c='b', label='Stellar Type 2') 
    ax1.set_ylabel('Stellar Type')
    ax1.set_yticks(range(useTypes.shape[0]))
    ax1.set_yticklabels([stellarTypes[typeNum] for typeNum in useTypes])

    ### Plot eccentricity
    handle3 = ax2.plot(Data['Time'][()], Data['Eccentricity'][()]-.01, linestyle='-', c='k', label= 'Eccentricity') # the minor subtraction makes the curve easier to find
    ax2.set_ylabel('Eccentricity', labelpad=10)
    ax2.set_yticks([0, .25, .5, .75, 1.0])
    ax2.set_ylim(-0.05, 1.05)
    ax2.tick_params(axis='y', left=False, right=True, direction='out', labelleft=False, labelright=True, pad=-25)

    # Legend
    handles, labels = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles.extend(handles2)
    labels.extend(labels2)
    ax.legend( handles=handles, labels=labels) 

    # Grid
    ax2.yaxis.grid(False)

    return handles, labels


def plotVanDenHeuval(events=None):
    num_events = len(events)
    fig, axs = plt.subplots(num_events, 1)
    plt.rcParams["text.usetex"] = True  # Use latex
    
    for ii in range(num_events):
        img = events[ii].eventImage
        axs[ii].imshow(img)
        axs[ii].set_xticks([])
        axs[ii].set_yticks([])
        axs[ii].yaxis.set_label_position("right")
        plt.subplots_adjust(hspace=0)

        pltString = "$t$ = {:.1f} Myr, $a = {:.1f}$ $R_\odot$ \n $M_1$ = {:.1f} $M_\odot$, $M_2$ = {:.1f} $M_\odot$ \n"+events[ii].eventString
        pltString = pltString.format(events[ii].time,events[ii].a,events[ii].m1,events[ii].m2)
        
        pad = 5
        axs[ii].annotate(pltString, xy=(0,0.5), xytext=(-axs[ii].yaxis.labelpad + pad,0),xycoords=axs[ii].yaxis.label,fontsize=8,textcoords='offset points', ha='left', va='center')
        axs[ii].annotate(chr(ord('@')+1+ii), xy=(-0.15,0.8),xycoords='axes fraction',fontsize=8,fontweight='bold')






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




###########################################################
### 
### Evolutionary Events
### 
###########################################################


def testing(Data):
    event = Event(Data, 0, 'Beg')
    img = event.getEventImage(2, True)

    fig, ax = plt.subplots()

    imgAspectWH = 1.7786666666666666 # inverse of above
    ax.set_xlim([0, imgAspectWH])
    ax.set_ylim([0, 1])
    x0, xF = ax.get_xlim()
    y0, yF = ax.get_ylim()
    ax.imshow(img, extent=(x0, xF, y0, yF)) # l, r, b, t 
    plt.show()


class Event(object):

    def __init__(self, Data, index, eventClass, **kwargs):

        self.Data   = Data
        self.index  = index
        self.eventClass = eventClass # Can be any of 'Beg', 'End', 'MT', 'SN', 'Stype'

        ii = index
        self.time   = Data['Time'][ii] 
        self.m1     = Data['Mass(1)'][ii]
        self.m2     = Data['Mass(2)'][ii] 
        self.stype1 = Data['Stellar_Type(1)'][ii] 
        self.stype2 = Data['Stellar_Type(2)'][ii] 
        self.a      = Data['SemiMajorAxis'][ii]
        self.e      = Data['Eccentricity'][ii]


        self.eventString = self.getEventDetails(**kwargs)
        #self.eventImage = getEventImage(image_num)
        #self.eventImage = self.getEventImage()


    def getEventDetails(self, **kwargs):
        """
        Use the event class and timestep, and possibly additional kwargs,
        to define the event string in a systematic way
        """
        eventClass = self.eventClass
        Data = self.Data
        ii = self.index
        rotate_image = False # Set to True if event goes from 2->1
        image_num = None

        if eventClass == 'Beg': # TODO
            eventString = r'Zero-age main-sequence'
            image_num = 2 

        elif eventClass == 'MT':
            mtValue = Data['MT_History'][ii]
            self.eventSubClass = mtValue
            
            if mtValue == 1:
                eventString = r'Stable mass transfer: 1 to 2'
                if self.stype2 < 13:
                    image_num = 26
                else:
                    image_num = 44
                    rotate_image = True
            elif mtValue == 2:
                eventString = r'Stable mass transfer: 2 to 1'
                if self.stype1 < 13:
                    image_num = 26
                    rotate_image = True
                else:
                    image_num = 44
            elif mtValue == 3:
                eventString = r'Common envelope initiated by 1'
                if (self.stype1 < 13) & (self.stype2 < 13):
                    image_num = 28
                else:
                    image_num = 49
            elif mtValue == 4:
                eventString = r'Common envelope initiated by 2'
                if (self.stype1 < 13) & (self.stype2 < 13):
                    image_num = 28
                else:
                    image_num = 49
            elif mtValue == 5:
                eventString = r'Double-core common envelope'
                image_num = 28
            elif mtValue == 6:
                eventString = r'CE: both MS'
                image_num = 28
            elif mtValue == 7:
                eventString = r'CE: MS with CO'
                image_num = 49
            else:
                raise ValueError("Unknown MT: {}".format(mtValue))

        elif eventClass == 'SN':
            whichStar = kwargs['whichStar']
            remnantType = 'NS' if (Data['Stellar_Type({})'.format(whichStar)][ii] == 13) else 'BH'
            compType = Data['Stellar_Type({})'.format(2 if whichStar==1 else 1)][ii]
            status = 'unbound' if (Data['Eccentricity'][ii]>1 or Data['SemiMajorAxis'][ii]<0) else 'intact'
            eventString = r'Star {} undergoes supernova and forms a {}, {}'.format(whichStar, remnantType, status)
            if compType < 13:
                image_num = 13 # 13 for normal companion
            else:
                image_num = 15 # 15 for CO companion

        elif eventClass == 'Stype':
            whichStar = kwargs['whichStar']
            stypePre = str(Data['Stellar_Type({})'.format(whichStar)][ii-1])
            stypePost= str(Data['Stellar_Type({})'.format(whichStar)][ii])
            eventString = r'Star {}: {}-$>${}'.format(whichStar, stypePre, stypePost)

        elif eventClass == 'End':
            state = kwargs['state']
            stype1 = Data['Stellar_Type(1)'][-1]
            stype2 = Data['Stellar_Type(2)'][-1]
            m1 = Data['Mass(1)'][-1]
            m2 = Data['Mass(2)'][-1]

            if state == 'DCO':
                Msunkg=1.98892e30
                c=299792458
                G=6.67428e-11
                Rsun = 695500000
                a = Data['SemiMajorAxis'][-1]*Rsun
                e = Data['Eccentricity'][-1]
                beta=64/5*G**3*m1*m2*(m1+m2)*Msunkg**3/c**5
                T0=a**4/4/beta 
                Tdelay=T0*(1-e**2)**(7/2)*(1+0.31*e**10 + 0.27*e**20 +  0.2*e**1000)/3.15e7/1e6
                eventString = r'Double compact object ({}+{}) merging in {:.1f} Myr'.format(stype1, stype1, Tdelay)

                if (stype1 == 13) & (stype2 == 13):
                    image_num = 55
                elif (stype1 == 14) & (stype2 == 14):
                    image_num = 51
                else:
                    image_num = 53

            elif state == "Unbound":
                eventString = r'Unbound: {}+{}'.format(stype1, stype2)
                image_num = None

            elif state == "Merger":
                mTot = m1+m2 
                mHe = Data['Mass_He_Core(1)'][-1] + Data['Mass_He_Core(2)'][-1]
                mCO = Data['Mass_CO_Core(1)'][-1] + Data['Mass_CO_Core(2)'][-1]
                eventString = r'Stellar Merger: {}+{}'.format(stype1, stype2)
                image_num = 37

            else:
                #raise ValueError("Unknown event state: {}".format(state))
                eventString = r'Unspecified endstate: {}+{}'.format(stype1, stype2)

        else:
            raise ValueError("Unknown event class: {}".format(self.eventClass))

        if image_num != None:
            self.eventImage = self.getEventImage(image_num, rotate_image)

        return eventString  

    def getEventImage(self, image_num, rotate_image):
        """
        Map the eventClass and possibly eventSubClass, with information
        on the stellar types, to get the van Den Heuval diagrams.
        """

        # self.imgFile = '/home/rwillcox/astro/compas/COMPAS/docs/media/vanDenHeuval_figures/{}.png'.format(image_num)
        self.imgFile = '/Users/13lauy1/git/compas2/COMPAS/docs/media/vanDenHeuval_figures/{}.png'.format(image_num)
        img = plt.imread(self.imgFile) # import image
        if rotate_image:
            img = img[:,::-1,:] # flip across y-axis
        return img







### Collecting events

def addEvent(allEvents, Data, ii, eventClass, **kwargs):

    newEvent = Event(Data, ii, eventClass, **kwargs)
    allEvents.append(newEvent)


def getAllEvents(Data):

    allEvents = []

    ### Add first timestep
    addEvent(allEvents, Data, 0, eventClass='Beg')

    ### Get all intermediary events
    for ii in range(Data['Time'].size):

        # Ignore first timestep, it's accounted for above
        if ii == 0:
            continue 

        # Note: These should all be if clauses, not elif/else, because they are not mutually exclusive

        ### Mass transfer happened
        if (Data['MT_History'][ii]>0) and not (Data['MT_History'][ii]==Data['MT_History'][ii-1]): # Not a repeated entry
            addEvent(allEvents, Data, ii, eventClass='MT')

        ### Type of star 1 changed
        if Data['Stellar_Type(1)'][ii]!=Data['Stellar_Type(1)'][ii-1]:    
            if (Data['Stellar_Type(1)'][ii] in [13, 14]): # SN star 1
                addEvent(allEvents, Data, ii, eventClass='SN', whichStar=1)
            else:
                addEvent(allEvents, Data, ii, eventClass='Stype', whichStar=1)

        ### Type of star 2 changed
        if Data['Stellar_Type(2)'][ii]!=Data['Stellar_Type(2)'][ii-1]:    
            if (Data['Stellar_Type(2)'][ii] in [13, 14]): # SN star 2
                addEvent(allEvents, Data, ii, eventClass='SN', whichStar=2)
            else:
                addEvent(allEvents, Data, ii, eventClass='Stype', whichStar=2)

    ### Add an event for final state of the binary
    isDCO = (Data['Stellar_Type(1)'][-1] in np.arange(10, 15)) and (Data['Stellar_Type(2)'][-1] in np.arange(10, 15)) # Both stars are WDs, NSs, or BHs
    isUnbound = (Data['Eccentricity'][-1]>1 or Data['SemiMajorAxis'][-1]<0)
    isMerger = (Data['Time'][-1]<14000) and not isDCO and not isUnbound     # System must have merged with at least one standard component

    if isDCO:
        state = "DCO" 
    elif isUnbound:
        state = "Unbound" 
    elif isMerger:
        state = "Merger"
    else:
        state = "Undef"
    addEvent(allEvents, Data, -1, eventClass='End', state=state)
    
    return allEvents


    
### Printing events


def printEvolutionaryHistory(Data=None, events=None):
    """
    This function prints a synopsys of the evolutionary history to the command line; it can eventually include cartoons as well.
    """
    if events != None:
        Data = events[0].Data
    elif Data != None:
        events = getAllEvents(Data)
    else:
        raise ValueError("No usable input given")


    print('Time (Myr), Event,                            M1 (M_o), type1, M2 (M_o), type2, a (R_o),   e')

    for event in events:
        ii = event.index
        printFormattedEvolutionLine( Data['Time'][ii], event.eventString.replace('$', ''), 
                                     Data['Mass(1)'][ii], Data['Stellar_Type(1)'][ii], 
                                     Data['Mass(2)'][ii], Data['Stellar_Type(2)'][ii], 
                                     Data['SemiMajorAxis'][ii], Data['Eccentricity'][ii])

def printFormattedEvolutionLine(time, event, m1, t1, m2, t2, a, e):
    # All values are floats except event which is a string and t1, t2 which are ints (stellar types)
    print("{:10.6f}   {:31}  {:7.3f}    {:2}    {:7.3f}    {:2}   {:8.3f}  {:5.3f}" .format(time, event, m1, t1, m2, t2, a, e))


if __name__ == "__main__":
    main()

