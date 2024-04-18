import h5py as h5
import numpy as np
import pandas as pd


########################################################################
### 
### Function to print the data from a given COMPAS HDF5 group 
### in a readable pandas template
### 
########################################################################

def printCompasDetails(data, *seeds, mask=()):
    """
    Function to print the full Compas output for given seeds, optionally with an additional mask
    """
    list_of_keys = list(data.keys())

    # Check if seed parameter exists - if not, just print without (e.g RunDetails)
    if ('SEED' in list_of_keys) | ('SEED>MT' in list_of_keys): # Most output files 
        #SEED>MT is a relic from older versions, but we leave this in for backwards compatibility

        # Set the seed name parameter, mask on seeds as needed, and set the index
        seedVariableName='SEED' if ('SEED' in list_of_keys) else 'SEED>MT'
        list_of_keys.remove(seedVariableName) # this is the index above, don't want to include it
    
        allSeeds = data[seedVariableName][()]
        seedsMask = np.in1d(allSeeds, seeds)
        if len(seeds) == 0: # if any seed is included, do not reset the mask
            seedsMask = np.ones_like(allSeeds).astype(bool)
        if mask == ():
            mask = np.ones_like(allSeeds).astype(bool)
        mask &= seedsMask

        df = pd.DataFrame.from_dict({param: data[param][()][mask] for param in list(data.keys())}).set_index(seedVariableName).T

    else: # No seed parameter, so do custom print for Run Details

        # Get just the keys without the -Derivation suffix - those will be a second column
        keys_not_derivations = []
        for key in list_of_keys:
            if '-Derivation' not in key:
                keys_not_derivations.append(key)
        
        # Some parameter values are string types, formatted as np.bytes_, need to convert back
        def convert_strings(param_array):
            if isinstance(param_array[0], np.bytes_):
                return param_array.astype(str)
            else:
                return param_array

        df_keys = pd.DataFrame.from_dict({param: convert_strings(data[param][()]) for param in keys_not_derivations }).T
        nCols = df_keys.shape[1] # Required only because if we combine RDs, we get many columns (should fix later)
        df_keys.columns = ['Parameter']*nCols
        df_drvs = pd.DataFrame.from_dict({param: convert_strings(data[param+'-Derivation'][()]) for param in keys_not_derivations }).T
        df_drvs.columns = ['Derivation']*nCols
        df = pd.concat([df_keys, df_drvs], axis=1)

    # Add units as first col
    units_dict = {key:data[key].attrs['units'].astype(str) for key in list_of_keys}
    df.insert(loc=0, column='(units)', value=pd.Series(units_dict))
    return df



########################################################################
### 
### Get event histories of MT data, SN data, and combined MT, SN data
### 
########################################################################

def getMtEvents(MT):                                     
    """
    This function takes in the `BSE_RLOF` output category from COMPAS, and returns the information
    on the Mass Transfer (MT) events that happen for each seed. The events do not have to be in order, 
    either chronologically or by seed, this function will reorder them as required.
    
    OUT:
        tuple of (returnedSeeds, returnedEvents, returnedTimes)
        returnedSeeds (list): ordered list of the unique seeds in the MT file
        returnedEvents (list): list of sublists, where each sublist contains all the MT events for a given seed.
            MT event tuples take the form :
            (stellarTypePrimary, stellarTypeSecondary, isRlof1, isRlof2, isCEE)
        returnTimes (list): is a list of sublists of times of each of the MT events
    """

    mtSeeds = MT['SEED'][()]
    mtTimes = MT['Time<MT'][()]
    mtPrimaryStype = MT['Stellar_Type(1)<MT'][()]
    mtSecondaryStype = MT['Stellar_Type(2)<MT'][()]
    mtIsRlof1 = MT['RLOF(1)>MT'][()] == 1
    mtIsRlof2 = MT['RLOF(2)>MT'][()] == 1
    mtIsCEE = MT['CEE>MT'][()] == 1

    # We want the return arrays sorted by seed, so sort here.
    mtSeedsInds = np.lexsort((mtTimes, mtSeeds)) # sort by seeds then times - lexsort sorts by the last column first...
    mtSeeds = mtSeeds[mtSeedsInds]  
    mtTimes = mtTimes[mtSeedsInds]  
    mtPrimaryStype = mtPrimaryStype[mtSeedsInds]
    mtSecondaryStype = mtSecondaryStype[mtSeedsInds]
    mtIsRlof1 = mtIsRlof1[mtSeedsInds]
    mtIsRlof2 = mtIsRlof2[mtSeedsInds]
    mtIsCEE = mtIsCEE[mtSeedsInds]

    # Process the MT events

    returnedSeeds = []                                      # array of seeds - will only contain seeds that have MT events
    returnedEvents = []                                     # array of MT events for each seed in returnedSeeds
    returnedTimes = []                                      # array of times for each event in returnedEvents (for each seed in returnedSeeds)

    lastSeed = -1                                           # initialize most recently processed seed

    for seedIndex, thisSeed in enumerate(mtSeeds):          # iterate over all RLOF file entries
        thisTime = mtTimes[seedIndex]                       # time for this RLOF file entry
        thisEvent = (mtPrimaryStype[seedIndex], mtSecondaryStype[seedIndex], 
                     mtIsRlof1[seedIndex], mtIsRlof2[seedIndex], mtIsCEE[seedIndex])  # construct event tuple

        # If this is an entirely new seed:
        if thisSeed != lastSeed:                            # same seed as last seed processed?
            returnedSeeds.append(thisSeed)                  # no - new seed, record it
            returnedTimes.append([thisTime])                # initialize the list of event times for this seed
            returnedEvents.append([thisEvent])              # initialize the list of events for this seed
            lastSeed = thisSeed                             # update the latest seed

        # Add event, if it is not a duplicate
        try:
            eventIndex = returnedEvents[-1].index(thisEvent)  # find eventIndex of this particular event tuple in the array of events for this seed 
            if thisTime > returnedTimes[-1][eventIndex]:      # ^ if event is not a duplicate, this will throw a ValueError 
                returnedTimes[-1][eventIndex] = thisTime      # if event is duplicate, update time to the later of the duplicates
        except ValueError:                                    # event is not a duplicate:
            returnedEvents[-1].append(thisEvent)              # record new event tuple for this seed
            returnedTimes[-1].append(thisTime)                # record new event time for this seed

    return returnedSeeds, returnedEvents, returnedTimes       # see above for description


def getSnEvents(SN):                                     
    """
    This function takes in the `BSE_Supernovae` output category from COMPAS, and returns the information
    on the Supernova (SN) events that happen for each seed. The events do not have to be in order chronologically,
    this function will reorder them as required.
    
    OUT:
        tuple of (returnedSeeds, returnedEvents, returnedTimes)     
        returnedSeeds (list): ordered list of all the unique seeds in the SN file
        returnedEvents (list): list of sublists, where each sublist contains all the SN events for a given seed.
            SN event tuples take the form :
            (stellarTypeProgenitor, stellarTypeRemnant, whichStarIsProgenitor, isBinaryUnbound)
        returnedTimes (list): is a list of sublists of times of each of the SN events
    """
    
    snSeeds = SN['SEED'][()]
    snTimes = SN['Time'][()]
    snProgStype = SN['Stellar_Type_Prev(SN)'][()]
    snRemnStype = SN['Stellar_Type(SN)'][()]
    snWhichProg = SN['Supernova_State'][()]
    snIsUnbound = SN['Unbound'][()] == 1

    # We want the return arrays sorted by seed, so sort here.
    snSeedsInds = np.lexsort((snTimes, snSeeds)) # sort by seeds then times - lexsort sorts by the last column first...
    snSeeds = snSeeds[snSeedsInds]  
    snTimes = snTimes[snSeedsInds]  
    snProgStype = snProgStype[snSeedsInds]
    snRemnStype = snRemnStype[snSeedsInds]
    snWhichProg = snWhichProg[snSeedsInds]
    snIsUnbound = snIsUnbound[snSeedsInds]
    
    # Process the SN events

    returnedSeeds = []                                      # array of seeds - will only contain seeds that have SN events
    returnedEvents = []                                     # array of SN events for each seed in returnedSeeds
    returnedTimes = []                                      # array of times for each event in returnedEvents (for each seed in returnedSeeds)

    lastSeed = -1                                           # initialize most recently processed seed

    for seedIndex, thisSeed in enumerate(snSeeds):          # iterate over all SN file entries
        thisTime = snTimes[seedIndex]                       # time for this SN file entry
        thisEvent = (snProgStype[seedIndex], snRemnStype[seedIndex], 
                     snWhichProg[seedIndex], snIsUnbound[seedIndex]) # construct event tuple
               
        # If this is an entirely new seed:
        if thisSeed != lastSeed:                            # same seed as last seed processed?
            returnedSeeds.append(thisSeed)                  # no - new seed, record it
            returnedTimes.append([thisTime])                #   initialize the list of event times for this seed
            returnedEvents.append([thisEvent])              #   initialize the list of events for this seed
            lastSeed = thisSeed                             #   update the latest seed
        else:                                               # yes - second SN event for this seed
            returnedTimes[-1].append(thisTime)              #   append time at end of array
            returnedEvents[-1].append(thisEvent)            #   append event at end of array
                
    return returnedSeeds, returnedEvents, returnedTimes     # see above for description


def getEventHistory(h5file, exclude_null=False):
    """
    Get the event history for all seeds, including both RLOF and SN events, in chronological order.
    IN:
        h5file (h5.File() type): COMPAS HDF5 output file
        exclude_null (bool): whether or not to exclude seeds which undergo no RLOF or SN events
    OUT:
        tuple of (returnedSeeds, returnedEvents)
        returnedSeeds (list): ordered list of all seeds in the output
        returnedEvents (list): a list of the collected SN and MT events from the 
            getMtEvents and getSnEvents functions above
    """

    SP = h5file['BSE_System_Parameters']
    MT = h5file['BSE_RLOF']
    SN = h5file['BSE_Supernovae']
    allSeeds = SP['SEED'][()]                                              # get all seeds
    mtSeeds, mtEvents, mtTimes = getMtEvents(MT)                           # get MT events
    snSeeds, snEvents, snTimes = getSnEvents(SN)                           # get SN events

    numMtSeeds = len(mtSeeds)                                               # number of MT events
    numSnSeeds = len(snSeeds)                                               # number of SN events

    if numMtSeeds < 1 and numSnSeeds < 1: return []                         # no events - return empty history

    returnedSeeds = []                                                      # array of seeds - will only contain seeds that have events (of any type)
    returnedEvents = []                                                     # array of events - same size as returnedSeeds (includes event times)

    eventOrdering = ['MT', 'SN']                                            # order of preference for simultaneous events

    mtIndex = 0                                                             # index into MT events arrays
    snIndex = 0                                                             # index into SN events arrays

    if exclude_null:
        seedsToIterate = np.sort(np.unique(np.append(mtSeeds, snSeeds)))    # iterate over all the seeds that have either MT or SN events
    else:
        seedsToIterate = allSeeds
        
    idxOrdered = np.argsort(seedsToIterate)
    returnedSeeds = [None] * np.size(seedsToIterate)                        # array of seeds - will only contain seeds that have events (of any type)
    returnedEvents = [None] * np.size(seedsToIterate)                       # array of events - same size as returnedSeeds (includes event times)

    for idx in idxOrdered:
        seed = seedsToIterate[idx]
        seedEvents = []                                                     # initialise the events for the seed being processed

        # Collect any MT events for this seed, add the time of the event and the event type
        while mtIndex < numMtSeeds and mtSeeds[mtIndex] == seed:
            for eventIndex, event in enumerate(mtEvents[mtIndex]):
                seedEvents.append(('MT', mtTimes[mtIndex][eventIndex], *mtEvents[mtIndex][eventIndex]))
            mtIndex += 1

        # Collect any SN events for this seed, add the time of the event and the event type 
        while snIndex < numSnSeeds and snSeeds[snIndex] == seed:
            for eventIndex, event in enumerate(snEvents[snIndex]):
                seedEvents.append(('SN', snTimes[snIndex][eventIndex], *snEvents[snIndex][eventIndex]))
            snIndex += 1

        seedEvents.sort(key=lambda ev:(ev[1], eventOrdering.index(ev[0])))  # sort the events by time and event type (MT before SN if at the same time)

        returnedSeeds[idx] = seed                                           # record the seed in the seeds array being returned
        returnedEvents[idx] = seedEvents                                    # record the events for this seed in the events array being returned

    return returnedSeeds, returnedEvents                                    # see above for details



###########################################
### 
### Produce strings of the event histories
### 
###########################################

def buildEventString(events):
    """
    Function to produce a string representing the event history of a single binary for quick readability.
    IN:
        events (list of tuples): events output from getEventHistory()
    OUT:
        eventString (string): string representing the event history of the binary
    
    MT strings look like: 
        P>S, P<S, or P=S where P is primary type, S is secondary type, 
        and >, < is RLOF (1->2 or 1<-2) or = for CEE

    SN strings look like:
        P*SR for star1 the SN progenitor,or 
        R*SP for star2 the SN progenitor,
        where P is progenitor type, R is remnant type, 
        S is state (I for intact, U for unbound)

    Event strings for the same seed are separated by the undesrcore character ('_')
    """
    
    # Empty event
    if len(events) == 0:
        return 'NA'
    
    eventStr = ''                                                          # event string for this star
    for event in events:
        if event[0] == 'MT':                                               # MT event
            eventStr += str(event[2])                                      # primary stellar type
            eventStr += '=' if event[6] else ('>' if event[4] else '<')    # event type: CEE, RLOF 1->2, RLOF 2->1
            eventStr += str(event[3])                                      # secondary stellar type
        else:                                                              # assume SN event (until other event types are added...)
            eventStr += str(event[2]) if event[4] == 1 else str(event[3])  # Progenitor or Remnant depending upon which star is the SN
            eventStr += '*U' if event[5] else '*I'                         # unbound or intact
            eventStr += str(event[3]) if event[4] == 1 else str(event[2])  # Progenitor or Remnant depending upon which star is the SN 

        eventStr += '_'                                                    # event separator

    return eventStr[:-1]                                                   # return event string for this star (pop the last underscore first)



def getEventStrings(h5file=None, allEvents=None):
    """
    Function to calculate the event history strings for either the entire Compas output, or some list of events
    IN: One of
        h5file (h5.File() type): COMPAS HDF5 output file
        allEvents (list of tuples)
    OUT:
        eventStrings (list): list of strings of the event history of each seed
    """

    # If output is 
    if (h5file == None) & (allEvents == None):
        return 
    elif (allEvents == None):
        _, allEvents = getEventHistory(h5file)
    
    eventStrings = []                                                        # array of event strings to be returned
    for eventsForGivenSeed in allEvents:  
        eventString = buildEventString(eventsForGivenSeed)
        eventStrings.append(eventString)                                     # append event string for this star (pop the last underscore first)

    return eventStrings
