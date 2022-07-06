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
    plotVanDenHeuvel(events=events)
    plt.savefig('vanDenHeuvelPlot.eps', bbox_inches='tight',pad_inches = 0, format='eps')
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

    listOfPlots = [ plotMassAttributes, plotLengthAttributes, plotStellarTypeAttributes, plotHertzsprungRussell]

    events = [event for event in events if event.eventClass != 'Stype'] # want to ignore simple stellar type changes
    if events[-1].eventClass == "End":
        events.pop()                                                    # don't include the 'End' eventClasses, they can compress the rest of the evolution too much
    num_events = len(events)
    event_times = [event.time for event in events]
    stopTimeAt = event_times[-1] * 1.05 # End time at the last event, plus 5% for convenience.
    if num_events == 1:
        stopTimeAt = Data['Time'][-1] * 1.05 # plot all the way to the end of the run if no events beyond ZAMS
    mask = Data['Time'][()] < stopTimeAt # Mask the data to not include the 'End' events

    rcParams.update(fontparams) # Set configurations for uniform plot output

    # Configure 2x2 subplots, for masses, lengths, stellar types, and HR diagram (in order top to bottom left to right)
    fig = plt.figure(figsize=(15,8)) # W, H
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2, sharex = ax1)
    ax3 = fig.add_subplot(2, 2, 3, sharex = ax1)
    ax4 = fig.add_subplot(2, 2, 4)
    axes = [ax1, ax2, ax3, ax4]

    for ii, specificPlot in enumerate(listOfPlots): # exclude the last one

        ax = axes[ii]

        # TODO: Set the reverse log scale for time

        # Plot the data
        handles, labels = specificPlot(ax=ax, Data=Data, events=events, mask=mask)

        # Add some breathing space at the top of the plot
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax*1.1)

        # Add vertical lines for specific event times
        [ax.axvline(time, ymin=0.975, zorder=0) for time in event_times]

        ### Time plots should have event letters spaced out, all plots should have both axes labelled and ticked
        # Add the event letters to the first plot
        if ii in [0, 1, 2]: # top plots
            spaced_out_event_times = space_out(event_times, min_separation=ax.get_xlim()[1]/75) # min_separation of xmax/50 was found to fit the letter sizes well
            for jj in range(num_events):
                yOffsetFactor = 1.5 if (ax.get_yscale() == 'log') else 1.02
                ax.text(x=spaced_out_event_times[jj], y=ax.get_ylim()[1]*yOffsetFactor, s=chr(ord('@')+1+jj)) # The unicode representation of the capital letters - works as long as there are less than 26 images to show
            ax.set_xlabel('Time / Myr')

        
        if handles is not None:
            ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1.03,0.5), fancybox=True)


    #### Finalize the boundaries, save, and show
    fig.suptitle('Detailed evolution for seed = {}'.format(Data['SEED'][()][0]), fontsize=18) 
    fig.tight_layout(h_pad=1, w_pad=1, rect= (0.08, 0.08, .98, .98), pad=0.) # (left, bottom, right, top) 
    plt.savefig('detailedEvolutionPlot.eps', bbox_inches='tight',pad_inches = 0, format='eps')





######## Plotting functions


def plotMassAttributes(ax=None, Data=None, mask=None, **kwargs):

    ### Plot mass attributes 
    # Create new column for total mass
    Mtot = Data['Mass(1)'][()][mask] + Data['Mass(2)'][()][mask]
    ax.plot(Data['Time'][()][mask], Mtot, linestyle='-', c='k', label='System Mass')
    ax.plot(Data['Time'][()][mask], Data['Mass(1)'][()][mask], linestyle='-', c='r', label='Total Mass 1')
    ax.plot(Data['Time'][()][mask], Data['Mass_He_Core(1)'][()][mask], linestyle='--', c='r', label='He Core 1')
    ax.plot(Data['Time'][()][mask], Data['Mass_CO_Core(1)'][()][mask], linestyle=':', c='r', label='CO Core 1')
    ax.plot(Data['Time'][()][mask], Data['Mass(2)'][()][mask], linestyle='-', c='b', label='Total Mass 2')
    ax.plot(Data['Time'][()][mask], Data['Mass_He_Core(2)'][()][mask], linestyle='--', c='b', label='He Core 2')
    ax.plot(Data['Time'][()][mask], Data['Mass_CO_Core(2)'][()][mask], linestyle=':', c='b', label='CO Core 2')

    ax.set_ylabel(r'Mass $/ \; M_{\odot}$')

    return ax.get_legend_handles_labels()
    
    
def plotLengthAttributes(ax=None, Data=None, mask=None, **kwargs):
          
    ### Plot radius attributes 
    ax.plot(Data['Time'][()][mask], Data['SemiMajorAxis'][()][mask], linestyle='-', c='k', label='Semi-Major Axis')
    ax.plot(Data['Time'][()][mask], Data['SemiMajorAxis'][()][mask]*(1-Data['Eccentricity'][()][mask]), linestyle=':', c='k', label='Periapsis')
    ax.plot(Data['Time'][()][mask], Data['Radius(1)'][()][mask], linestyle='-', c='r', label='Stellar Radius 1')
    ax.plot(Data['Time'][()][mask], Data['Radius(2)'][()][mask], linestyle='-', c='b', label='Stellar Radius 2')
    # Need to mask out when the denominator is 0
    mask1 = mask & (Data['Radius(1)|RL'][()] != 0)
    ax.plot(Data['Time'][()][mask1], Data['Radius(1)'][()][mask1]/Data['Radius(1)|RL'][()][mask1], linestyle='--', c='r', label='Roche Radius 1')
    mask2 = mask & (Data['Radius(2)|RL'][()] != 0)
    ax.plot(Data['Time'][()][mask2], Data['Radius(2)'][()][mask2]/Data['Radius(2)|RL'][()][mask2], linestyle='--', c='b', label='Roche Radius 2')

    ax.set_ylabel(r'Radius $/ \; R_{\odot}$')
    ax.set_yscale('log')

    return ax.get_legend_handles_labels()
    

def plotEccentricity(ax=None, Data=None, mask=None, **kwargs):

    ### Plot eccentricity
    ax.plot(Data['Time'][()], Data['Eccentricity'][()], linestyle='-', c='k') #, label= 'Eccentricity')
    ax.set_ylabel('Eccentricity')

    ax.set_ylim(-0.05, 1.05)
    ax.grid(linestyle=':', c='gray')
    
    return None, None
    
def plotStellarTypeAttributes(ax=None, Data=None, mask=None, **kwargs):

    ### Plot stellar types
    stellarTypes, useTypes, typeNameMap = getStellarTypes(Data)
    
    ax.plot(Data['Time'][()][mask], typeNameMap(Data['Stellar_Type(1)'][()][mask]), linestyle='-', c='r', label='Stellar Type 1')
    ax.plot(Data['Time'][()][mask], typeNameMap(Data['Stellar_Type(2)'][()][mask]), linestyle='-', c='b', label='Stellar Type 2')
    ax.set_ylabel('Stellar Type')
    ax.legend(prop={'size':8}) #, loc='lower left')
    ax.set_yticks(range(useTypes.shape[0]))
    ax.set_yticklabels([stellarTypes[typeNum] for typeNum in useTypes])
    ax.yaxis.grid(True)

    ax.legend(framealpha=1, prop={'size':8} ) 
    ax.grid(linestyle=':', c='gray')

    return ax.get_legend_handles_labels()
    

def plotStellarTypeAttributesAndEccentricity(ax=None, Data=None, mask=None, **kwargs):

    ax1 = ax
    ax2 = ax.twinx()
    
    ### Plot stellar types
    stellarTypes, useTypes, typeNameMap = getStellarTypes(Data)
    
    handle1 = ax1.plot(Data['Time'][()][mask], typeNameMap(Data['Stellar_Type(1)'][()][mask]), linestyle='-', c='r', label='Stellar Type 1')
    handle2 = ax1.plot(Data['Time'][()][mask], typeNameMap(Data['Stellar_Type(2)'][()][mask]), linestyle='-', c='b', label='Stellar Type 2') 
    ax1.set_ylabel('Stellar Type')
    ax1.set_yticks(range(useTypes.shape[0]))
    ax1.set_yticklabels([stellarTypes[typeNum] for typeNum in useTypes])

    ### Plot eccentricity
    handle3 = ax2.plot(Data['Time'][()][mask], Data['Eccentricity'][()][mask]-.01, linestyle='-', c='k', label= 'Eccentricity') # the minor subtraction makes the curve easier to find
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




def plotHertzsprungRussell(ax=None, Data=None, events=None, mask=None, **kwargs):

    ### Plot HR diagram: L vs Teff

    #Data['Teff(1)'][()] #K
    #Data['Luminosity(1)'][()] # Lsol
    
    ax.plot(Data['Teff(1)'][()][mask], Data['Luminosity(1)'][()][mask], linestyle='-', c='r', label='Star 1')
    ax.plot(Data['Teff(2)'][()][mask], Data['Luminosity(2)'][()][mask], linestyle='-', c='b', label='Star 2')
    ax.set_xlabel(r'Temperature [log(T/K)]')
    ax.set_ylabel(r'Luminosity [log($L/L_\odot$)]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim([min(1e3,xlim[0]), max(3e6,xlim[1])])
    ax.set_ylim([min(1e-4,ylim[0]), max(1e6,ylim[1])])
    xlim = ax.get_xlim() 
    ylim = ax.get_ylim()
    ax.invert_xaxis()

    # Add lines of const radii
    for R in np.logspace(-9, 5, 15):
        exp = "{:.1e}".format(R)
        exp = exp[-3] + exp[-1]
        if ((int(exp)%2)==1): # skip odd ones to remove clutter
            continue
        T_K = np.logspace(3, 7, 41) # in K
        T = T_K/6e3  # Tsol=6e3K
        def get_L(t): # assumes K
            return R*R *t*t*t*t
        L = get_L(T)
        ax.plot(T_K,L, '--k', alpha=0.2)
        # Plot the Rsol text at the bottom and right
        Lbot = ylim[0]*8 #Lsun  -2
        Trgt = xlim[0]*2 #3e3
        Tbot = np.sqrt(np.sqrt(Lbot/(R*R)))*6e3 # K
        Lrgt = get_L(Trgt/6e3)
        alpha=0.4
        if (Tbot > Trgt) and (Tbot < xlim[1]):
            ax.text(x=Tbot, y=Lbot, s=r"$R_\odot^{{{exp}}}$".format(exp=exp), alpha=alpha)
        elif (Lrgt > Lbot) and (Lrgt < ylim[1]):
            ax.text(x=Trgt, y=Lrgt, s=r"$R_\odot^{{{exp}}}$".format(exp=exp), alpha=alpha)

    # Add in the letters corresponding to various events
    event_times = [event.time for event in events]
    mask2 = mask & (np.in1d(Data['Time'][()], event_times))
    Tmsk = Data['Teff(1)'][()][mask2]
    Lmsk = Data['Luminosity(1)'][()][mask2]
    for jj in range(np.sum(mask2)):
        ax.text(x=Tmsk[jj], y=Lmsk[jj], s=chr(ord('@')+1+jj)) # The unicode representation of the capital letters - works as long as there are less than 26 images to show

    ax.legend(framealpha=1, prop={'size':8} ) 
    ax.grid(linestyle=':', c='gray')

    return ax.get_legend_handles_labels()



def plotVanDenHeuvel(events=None):
    # Only want events with an associated image
    events = [event for event in events if (event.eventImage is not None)]
    num_events = len(events)
    fig, axs = plt.subplots(num_events, 1)
    if num_events == 1:
        axs = [axs]
    fig.set_figwidth(9)
    plt.rcParams["text.usetex"] = True  # Use latex
    
    for ii in range(num_events):
        img = events[ii].eventImage
        axs[ii].imshow(img)
        axs[ii].set_xticks([])
        axs[ii].set_yticks([])
        axs[ii].yaxis.set_label_position("right")
        plt.subplots_adjust(hspace=0)

        if ii==0:
            pltString = "$t$ = {:.1f} Myr, $a$ = {:.1f} $R_\odot$ \n $M_1$ = {:.1f} $M_\odot$, $M_2$ = {:.1f} $M_\odot$ \n"+events[ii].eventString
            pltString = pltString.format(events[ii].time, events[ii].a,events[ii].m1,events[ii].m2)
        else:
            pltString = "$t$ = {:.1f} Myr, $a$ = {:.1f} to {:.1f} $R_\odot$ \n $M_1$ = {:.1f} to {:.1f} $M_\odot$, $M_2$ = {:.1f} to {:.1f} $M_\odot$ \n"+events[ii].eventString
            pltString = pltString.format(events[ii].time,events[ii].aprev, events[ii].a,events[ii].m1prev,events[ii].m1,events[ii].m2prev,events[ii].m2)
        
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
    stellarTypes = [r'MS$<0.7M_\odot$', r'MS$\geq0.7M_\odot$', 'HG', 'FGB', 'CHeB', 'EAGB', 'TPAGB', 'HeMS', 'HeHG', 'HeGB', 'HeWD', 'COWD', 'ONeWD', 'NS', 'BH', 'MR']

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

    if len(vals) > 1:
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
        self.a      = Data['SemiMajorAxis'][ii]
        if ii==0:
            self.m1prev=0
            self.m2prev=0
            self.aprev=0
        else:
            self.m1prev = Data['Mass(1)'][ii-1]
            self.m2prev = Data['Mass(2)'][ii-1]
            self.aprev  = Data['SemiMajorAxis'][ii-1]
        # Cheap kludge for SN separations - later, should clean up when detailed printing is called
        if eventClass == 'SN':
            try: # Bad form to do a try except, but this works for now
                self.a = Data['SemiMajorAxis'][ii+1]
                self.aprev = Data['SemiMajorAxis'][ii]
            except:
                self.a = Data['SemiMajorAxis'][ii]
                self.aprev = Data['SemiMajorAxis'][ii-1]
        self.stype1 = Data['Stellar_Type(1)'][ii]
        self.stype2 = Data['Stellar_Type(2)'][ii] 
        self.stypeName1 = stellarTypeMap[self.stype1]
        self.stypeName2 = stellarTypeMap[self.stype2]
        self.e      = Data['Eccentricity'][ii]
        self.Z1     = Data['Metallicity@ZAMS(1)'][ii]

        self.eventImage = None
        self.endState = None # sets the endstate - only relevant if eventClass=='End'
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
                self.eventClass = 'End'
                self.endState = 'Merger'
                eventString = r'Stellar Merger: {}+{}'.format(self.stypeName1, self.stypeName2)
                image_num = 37
            else:
                raise ValueError("Unknown MT: {}".format(mtValue))

        elif eventClass == 'SN':
            whichStar = kwargs['whichStar']
            remType = Data['Stellar_Type({})'.format(whichStar)][ii]
            remnantTypeName = self.stellarTypeMap[remType] 
            compType = Data['Stellar_Type({})'.format(2 if whichStar==1 else 1)][ii]
            disrupted = (Data['Eccentricity'][ii]>1 or Data['SemiMajorAxis'][ii]<0)
            status = '. Orbit becomes unbound' if disrupted else ''
            eventString = r'Star {} undergoes supernova and forms a {}{}'.format(whichStar, remnantTypeName, status)
            if disrupted: 
                if compType < 13:        # normal companion
                    if remType == 13:    # with NS
                        image_num = 19 
                    else:                # with BH
                        image_num = 20   
                elif compType == 13:     # NS companion
                    if remType == 13:    # with NS
                        image_num = 22 
                    else:                # with BH
                        image_num = 21   
                else:                    # BH companion 
                    if remType == 13:    # with NS
                        image_num = 24   
                    else:                # with BH
                        image_num = 23   
            else:
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
            self.endState = state
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
                eventString = r'Double compact object ({}+{}) merging in {:.2e} Myr'.format(self.stypeName1, self.stypeName2, Tdelay)

                if (stype1 == 13) & (stype2 == 13):
                    image_num = 55
                elif (stype1 == 14) & (stype2 == 14):
                    image_num = 51
                else:
                    image_num = 53

            elif state == "Unbound":
                eventString = r'Unbound: {}+{}'.format(self.stypeName1, self.stypeName2)
                image_num = None

            else:
                eventString = r'Evolution ended by run duration: {}+{}'.format(self.stypeName1, self.stypeName2) 
                image_num = 2

        else:
            raise ValueError("Unknown event class: {}".format(self.eventClass))

        if image_num != None:
            self.eventImage = self.getEventImage(image_num, rotate_image)

        return eventString  

    def getEventImage(self, image_num, rotate_image):
        """
        Map the eventClass and possibly eventSubClass, with information
        on the stellar types, to get the van den Heuvel diagrams.
        """

        self.imgFile = compasRootDir + 'utils/media/vanDenHeuvel_figures/{}.png'.format(image_num)
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
        isMerger = False
        for ii in range(Data['Time'].size):
    
            # Ignore first timestep, it's accounted for above
            if ii == 0:
                continue 
            # if ending has come, break out
            if isMerger:
                break
    
            # Note: These should all be if clauses, not elif/else, because they are not mutually exclusive
    
            ### Mass transfer happened
            if (Data['MT_History'][ii]>0) and not (Data['MT_History'][ii]==Data['MT_History'][ii-1]): # Not a repeated entry
                isMerger = self.addEvent(ii, eventClass='MT') # if a stellar merger, the eventClass changes
    
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
        if not isMerger: # set if a merger was flagged earlier
            isUnbound = (Data['Eccentricity'][-1]>1 or Data['SemiMajorAxis'][-1]<0)
            isDCO = (Data['Stellar_Type(1)'][-1] in np.arange(10, 15)) and (Data['Stellar_Type(2)'][-1] in np.arange(10, 15)) # Both stars are WDs, NSs, or BHs
        
            if isUnbound:
                state = "Unbound" 
            elif isDCO:
                state = "DCO" 
            else:
                state = "Undef"
            self.addEvent(-1, eventClass='End', state=state)
        
        return allEvents


    def addEvent(self, ii, eventClass, **kwargs):
    
        newEvent = Event(self.Data, ii, eventClass, self.stellarTypeMap, **kwargs)
        self.allEvents.append(newEvent)
        return newEvent.endState == 'Merger'



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

