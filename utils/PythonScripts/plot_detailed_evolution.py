###################################################################
#                                                                 #                                                               
#  Example of plotting detailed output COMPAS with python         #
#                                                                 #
###################################################################

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

compasRootDir = os.path.expandvars(os.environ['COMPAS_ROOT_DIR'])

def main():
    ### Read file and create dataframe.
    try:
        optional_input = sys.argv[1] 
        if optional_input is not None:
            data_path = optional_input
    except IndexError: # default
        data_path = 'COMPAS_Output/Detailed_Output/BSE_Detailed_Output_0.h5'

    Data = h5.File(data_path, 'r')

    ### Collect the important events in the detailed evolution
    events = allEvents(Data).allEvents # Calculate the events here, for use in plot sizing parameters
    printEvolutionaryHistory(events=events)
    events = [event for event in events if event.eventClass != 'Stype'] # want to ignore simple stellar type changes

    ### Produce the two plots
    makeDetailedPlots(Data, events)
    plotVanDenHeuval(events=events)
    plt.savefig('vanDenHeuvalPlot.eps', bbox_inches='tight',pad_inches = 0, format='eps')
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
}



####### Functions to organize and call the plotting functions 

def makeDetailedPlots(Data=None, events=None):

    listOfPlots = [ plotMassAttributes, plotLengthAttributes, plotStellarTypeAttributes, plotEccentricity]

    events = [event for event in events if event.eventClass != 'Stype'] # want to ignore simple stellar type changes
    num_events = len(events)
    event_times = [event.time for event in events]


    rcParams.update(fontparams) # Set configurations for uniform plot output
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 8)) # W, H

    for ii, specificPlot in enumerate(listOfPlots): # exclude the last one

        ax = axes.flatten()[ii]

        # TODO: Set the reverse log scale for time

        # Plot the data
        handles, labels = specificPlot(fig=fig, ax=ax, Data=Data)

        # Add some breathing space at the top of the plot
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax*1.1)

        # Add vertical lines for specific event times
        [ax.axvline(time, ymin=0.975, zorder=0) for time in event_times]

        ### Top plots should have event letters spaced out, bottom plots should have Time label and tick labels
        # Add the event letters to the first plot
        if ii in [0, 1]: # top plots
            spaced_out_event_times = space_out(event_times, min_separation=ax.get_xlim()[1]/75) # min_separation of xmax/50 was found to fit the letter sizes well
            for jj in range(num_events):
                yOffsetFactor = 1.5 if (ax.get_yscale() == 'log') else 1.02
                ax.text(x=spaced_out_event_times[jj], y=ax.get_ylim()[1]*yOffsetFactor, s=chr(ord('@')+1+jj)) # The unicode representation of the capital letters - works as long as there are less than 26 images to show
            ax.axes.xaxis.set_ticklabels([])
        else: # bottom plots
            ax.set_xlabel('Time / Myr')
        
        if handles is not None:
            ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1.03,0.5), fancybox=True)


    #### Finalize the boundaries, save, and show
    fig.suptitle('Detailed evolution for seed = {}'.format(Data['SEED'][()][0]), fontsize=18) 
    fig.tight_layout(h_pad=1, w_pad=1, rect= (0.08, 0.08, .98, .98), pad=0.) # (left, bottom, right, top) 
    plt.savefig('detailedEvolutionPlot.eps', bbox_inches='tight',pad_inches = 0, format='eps')





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
    ax.grid(linestyle=':', c='gray')
    
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
    # Only want events with an associated image
    events = [event for event in events if (event.eventImage is not None)]
    num_events = len(events)
    fig, axs = plt.subplots(num_events, 1)
    fig.set_figwidth(9)
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


def space_out(original_vals, min_separation=None):
    """
    This function takes a sorted array of floats (in this case, event times)
    and spaces them out from each other to have a minimum separation min_separation.
    
    The purpose of this is so that the event letters don't overlap on the plot.

    Idea of this function: for each pair which is too close, subtract off 
    some amount (nudge) from the lower, add the same amount to the upper, 
    and do this over the whole range. For big clumps, the middle ones won't move 
    (+/- will cancel out), but as the outer ones move away, the inner ones will 
    have room to stretch out.
    """

    vals = np.array(original_vals).copy()
    if min_separation == None:
        min_separation = np.max(vals)/50 # is this a good value?
    nudge = min_separation/10 # keep the nudge small so that you don't overdo the jump

    while(np.min(np.diff(vals))) < min_separation:
        maskTooClose = np.diff(vals) < min_separation
        vals[:-1][maskTooClose] -= nudge
        vals[1:][maskTooClose]  += nudge

    return vals
        


###########################################################
### 
### Evolutionary Events
### 
###########################################################


class Event(object):

    def __init__(self, Data, index, eventClass, stellarTypeMap, **kwargs):

        self.Data   = Data
        self.index  = index
        self.eventClass = eventClass # Can be any of 'Beg', 'End', 'MT', 'SN', 'Stype'
        self.stellarTypeMap = stellarTypeMap

        ii = index
        self.time   = Data['Time'][ii] 
        self.m1     = Data['Mass(1)'][ii]
        self.m2     = Data['Mass(2)'][ii] 
        self.stype1 = Data['Stellar_Type(1)'][ii] 
        self.stype2 = Data['Stellar_Type(2)'][ii] 
        self.stypeName1 = stellarTypeMap[self.stype1]
        self.stypeName2 = stellarTypeMap[self.stype2]
        self.a      = Data['SemiMajorAxis'][ii]
        self.e      = Data['Eccentricity'][ii]
        self.Z1     = Data['Metallicity@ZAMS(1)'][ii]

        self.eventImage = None
        self.eventString = self.getEventDetails(**kwargs)


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

        if eventClass == 'Beg': 
            eventString = r'Zero-age main-sequence, metallicity Z={:5.4f}'.format(self.Z1)
            image_num = 2 

        elif eventClass == 'MT':
            mtValue = Data['MT_History'][ii]
            self.eventSubClass = mtValue
            
            if mtValue == 1:
                eventString = r'Stable mass transfer from 1 to 2'
                if self.stype2 < 13:
                    image_num = 26
                else:
                    image_num = 44
                    rotate_image = True
            elif mtValue == 2:
                eventString = r'Stable mass transfer from 2 to 1'
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
            remnantTypeName = self.stellarTypeMap[Data['Stellar_Type({})'.format(whichStar)][ii]] 
            compType = Data['Stellar_Type({})'.format(2 if whichStar==1 else 1)][ii]
            status = '. Orbit becomes unbound' if (Data['Eccentricity'][ii]>1 or Data['SemiMajorAxis'][ii]<0) else ''
            eventString = r'Star {} undergoes supernova and forms a {}{}'.format(whichStar, remnantTypeName, status)
            if compType < 13:
                image_num = 13 # 13 for normal companion
            else:
                image_num = 15 # 15 for CO companion

        elif eventClass == 'Stype':
            whichStar = kwargs['whichStar']
            stypePre = self.stellarTypeMap[Data['Stellar_Type({})'.format(whichStar)][ii-1]]
            stypePost= self.stellarTypeMap[Data['Stellar_Type({})'.format(whichStar)][ii]]
            eventString = r'Star {}: {}-$>${}'.format(whichStar, stypePre, stypePost)

        elif eventClass == 'End':
            state = kwargs['state']
            stype1 = self.stype1 
            stype2 = self.stype2 
            m1     = self.m1 
            m2     = self.m2 

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
                eventString = r'Double compact object ({}+{}) merging in {:.1f} Myr'.format(self.stypeName1, self.stypeName2, Tdelay)

                if (stype1 == 13) & (stype2 == 13):
                    image_num = 55
                elif (stype1 == 14) & (stype2 == 14):
                    image_num = 51
                else:
                    image_num = 53

            elif state == "Unbound":
                eventString = r'Unbound: {}+{}'.format(self.stypeName1, self.stypeName2)
                image_num = None

            elif state == "Merger":
                mTot = m1+m2 
                mHe = Data['Mass_He_Core(1)'][-1] + Data['Mass_He_Core(2)'][-1]
                mCO = Data['Mass_CO_Core(1)'][-1] + Data['Mass_CO_Core(2)'][-1]
                eventString = r'Stellar Merger: {}+{}'.format(self.stypeName1, self.stypeName2)
                image_num = 37

            else:
                #raise ValueError("Unknown event state: {}".format(state))
                eventString = r'Unspecified endstate: {}+{}'.format(self.stypeName1, self.stypeName2) 

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

        self.imgFile = compasRootDir + 'utils/media/vanDenHeuval_figures/{}.png'.format(image_num)
        img = plt.imread(self.imgFile) # import image
        if rotate_image:
            img = img[:,::-1,:] # flip across y-axis
        return img



class allEvents(object):
    def __init__(self, Data):

        self.Data   = Data
        self.stellarTypeMap, _, _ = getStellarTypes(Data)

        # Collect all events into allEvents list
        self.allEvents = []
        self.getAllEvents()


    def getAllEvents(self):
    
        Data = self.Data
    
        ### Add first timestep
        self.addEvent(0, eventClass='Beg')
    
        ### Get all intermediary events
        for ii in range(Data['Time'].size):
    
            # Ignore first timestep, it's accounted for above
            if ii == 0:
                continue 
    
            # Note: These should all be if clauses, not elif/else, because they are not mutually exclusive
    
            ### Mass transfer happened
            if (Data['MT_History'][ii]>0) and not (Data['MT_History'][ii]==Data['MT_History'][ii-1]): # Not a repeated entry
                self.addEvent(ii, eventClass='MT')
    
            ### Type of star 1 changed
            if Data['Stellar_Type(1)'][ii]!=Data['Stellar_Type(1)'][ii-1]:    
                if (Data['Stellar_Type(1)'][ii] in [13, 14]): # SN star 1
                    self.addEvent(ii, eventClass='SN', whichStar=1)
                else:
                    self.addEvent(ii, eventClass='Stype', whichStar=1)
    
            ### Type of star 2 changed
            if Data['Stellar_Type(2)'][ii]!=Data['Stellar_Type(2)'][ii-1]:    
                if (Data['Stellar_Type(2)'][ii] in [13, 14]): # SN star 2
                    self.addEvent(ii, eventClass='SN', whichStar=2)
                else:
                    self.addEvent(ii, eventClass='Stype', whichStar=2)
    
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
        self.addEvent(-1, eventClass='End', state=state)
        
        return allEvents


    def addEvent(self, ii, eventClass, **kwargs):
    
        newEvent = Event(self.Data, ii, eventClass, self.stellarTypeMap, **kwargs)
        self.allEvents.append(newEvent)



###########################################################
### 
### Printing events
### 
###########################################################


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

