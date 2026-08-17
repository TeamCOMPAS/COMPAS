import h5py as h5
import numpy as np
import pandas as pd
from typing import NewType

try:
    from numpy.dtypes import StringDType
except ImportError:
    StringDType = getattr(np.dtypes, "StringDType", np.str_)

### New Types
MaskNdarray = NewType('MaskNdarray', np.ndarray[bool])
H5File = NewType('H5File', h5._hl.files.File)
H5Group = NewType('H5Group', h5._hl.group.Group)
MtEventTuple = NewType('MtEventTuple',
                       tuple[np.ndarray[int],
                             np.ndarray[int],
                             np.ndarray[bool],
                             np.ndarray[bool],
                             np.ndarray[bool],
                             np.ndarray[bool]])
SnEventTuple = NewType('SnEventTuple',
                       tuple[np.ndarray[int],
                             np.ndarray[int],
                             np.ndarray[bool],
                             np.ndarray[bool]])
EventHistoryString = NewType('EventHistoryString', np.array(np.str_))


### Function to print the data from a given COMPAS HDF5 group in a readable pandas template

def convert_bytes_array_to_strings(param_array):
    """Check and convert np.bytes_ array to strings.

    Parameters
    ----------
    param_array : np.ndarray[np.bytes_ | str]
        a numpy array containing either np.bytes_ or strings

    Returns
    -------
    param_array : np.ndarray[str]
        the same array, but converted to strings if not already
    """
    if np.issubdtype(param_array.dtype, np.bytes_):
        return param_array.astype(str)
    else:
        return param_array


def print_compas_details_dataframe(data: H5Group,
                                   *seeds: int,
                                   mask: MaskNdarray = ()) -> pd.DataFrame:
    """Return a pd.DataFrame for the COMPAS output contained in an H5Group.

    Parameters
    ----------
    data : H5Group
        an H5Group containing COMPAS data
    seeds : int or array_like of ints, optional
        seeds of specific systems to include (default is to include all seeds)
    mask : array_like of booleans, optional
        boolean mask to filter out specific systems (default is no mask)
        Must be same length as data['SEED']

    Returns
    -------
    compas_details : pd.DataFrame
        The values of each parameter in the H5Group, excluding those removed by `seeds` or `mask`.

    Notes
    -----
    Designed for easier visualization in a jupyter notebook, for inspection / debugging purposes.
    If both `seeds` and `mask` are supplied, resulting systems must pass both filters.

    Examples
    --------
    >>> mt_data = h5.File('COMPAS_Output.h5')['BSE_RLOF']
    >>> print_compas_details_dataframe(mt_data)
    [output of all BSE_RLOF events]

    >>> mt_seeds = mt_data['SEED'][()]
    >>> print_compas_details_dataframe(mt_data, mt_seeds[:50])
    [output of all BSE_RLOF events occurring in the first 50 seeds]

    >>> cee_events = mt_data['CEE>MT'][()] == 1 # needed to convert to boolean mask
    >>> print_compas_details_dataframe(mt_data, mt_seeds[:50], mask=cee_events)
    [output of all Common Envelope events occurring in the first 50 seeds]
    """
    # Check if SEED parameter exists in data
    if 'SEED' in data or 'SEED>MT' in data:
        # SEED>MT is a relic from older versions, but we leave this in for
        # backwards compatibility
        seed_variable_name = 'SEED' if 'SEED' in data else 'SEED>MT'

        # If `seeds` or `mask` arguments supplied, create the relevant mask
        all_seeds = np.asarray(data[seed_variable_name][()])
        seeds_mask = np.isin(all_seeds, seeds)
        if len(seeds) == 0:  # If `seeds` argument is not supplied, set the default mask
            seeds_mask = np.ones_like(all_seeds, dtype=bool)

        if isinstance(mask, tuple) and len(mask) == 0:
            mask = np.ones_like(all_seeds, dtype=bool)
        else:
            mask = np.asarray(mask, dtype=bool)
            if mask.size == 0:
                mask = np.ones_like(all_seeds, dtype=bool)
            if mask.shape != all_seeds.shape:
                raise ValueError(
                    "mask must have the same length as the data['SEED'] array"
                )

        mask = mask & seeds_mask
        df = pd.DataFrame.from_dict(
            {param: data[param][()][mask] for param in data}).set_index(seed_variable_name).T

    else:  # No seed parameter, currently only applies to the RunDetails H5Group
        # Get just the keys without the -Derivation suffix
        keys_not_derivations = [
            key for key in data if '-Derivation' not in key]

        df_keys = pd.DataFrame.from_dict({param: convert_bytes_array_to_strings(
            data[param][()]) for param in keys_not_derivations}).T
        # Required only because if we have combined RunDetails output from a
        # previous h5copy, we get many columns (should fix later)
        nCols = df_keys.shape[1]
        df_keys.columns = ['Parameter'] * nCols
        df_drvs = pd.DataFrame.from_dict({param: convert_bytes_array_to_strings(
            data[param + '-Derivation'][()]) for param in keys_not_derivations}).T
        df_drvs.columns = ['Derivation'] * nCols
        df = pd.concat([df_keys, df_drvs], axis=1)

    # Add units as first col
    units_dict = {
        key: str(data[key].attrs['units']) if 'units' in data[key].attrs else ''
        for key in data
    }
    df.insert(loc=0, column='(units)', value=pd.Series(units_dict))
    return df


### Get event histories of MT data, SN data, and combined MT, SN data

def get_mt_data_tuple(mt_data: H5Group) -> tuple[list, list, list]:
    """Calculates the EventTuple for the BSE_RLOF output H5Group.

    Parameters
    ----------
    mt_data : H5Group
        the COMPAS output H5Group corresponding to BSE_RLOF

    Returns
    -------
    returned_seeds : list
        an ordered list of the unique seeds in the mt_data file
    returned_events : list
        a list of sublists, one sublist per seed, where each sublist
        contains all the MtEventTuples for the given seed
    returned_times : list
        a list of sublists of times of each of the mt_data events

    Notes
    -----
    MtEventTuples take the form: (stellar_type_primary, stellar_type_secondary, is_rlof1, is_rlof2, is_cee, is_mrg).
    The events in the input do not have to be ordered chronologically, this function orders them, in the event
    that the input was coallated from a previous h5copy command. The output is not meant to be used generally,
    it is just a feeder for get_event_history below.
    """
    mt_seeds = mt_data['SEED'][()]
    mt_times = mt_data['Time<MT'][()]
    mt_primary_stype = mt_data['Stellar_Type(1)<MT'][()]
    mt_secondary_stype = mt_data['Stellar_Type(2)<MT'][()]
    mt_is_rlof1 = mt_data['RLOF(1)>MT'][()] == 1
    mt_is_rlof2 = mt_data['RLOF(2)>MT'][()] == 1
    mt_is_cee = mt_data['CEE>MT'][()] == 1
    mt_is_mrg = mt_data['Merger'][()] == 1

    # We want the return arrays sorted by seed, so sort here.
    # sort by seeds then times - lexsort sorts by the last column first...
    mt_seeds_idx = np.lexsort((mt_times, mt_seeds))
    mt_seeds = mt_seeds[mt_seeds_idx]
    mt_times = mt_times[mt_seeds_idx]
    mt_primary_stype = mt_primary_stype[mt_seeds_idx]
    mt_secondary_stype = mt_secondary_stype[mt_seeds_idx]
    mt_is_rlof1 = mt_is_rlof1[mt_seeds_idx]
    mt_is_rlof2 = mt_is_rlof2[mt_seeds_idx]
    mt_is_cee = mt_is_cee[mt_seeds_idx]
    mt_is_mrg = mt_is_mrg[mt_seeds_idx]

    # Process the mt_data events
    # array of seeds - will only contain seeds that have mt_data events
    returned_seeds = []
    # array of mt_data events for each seed in returned_seeds
    returned_events = []
    # array of times for each event in returned_events (for each seed in
    # returned_seeds)
    returned_times = []

    # initialize most recently processed seed
    last_seed = -1
    for seed_index, this_seed in enumerate(
            mt_seeds):          # iterate over all RLOF file entries
        # time for this RLOF file entry
        this_time = mt_times[seed_index]
        this_event = (
            mt_primary_stype[seed_index],
            mt_secondary_stype[seed_index],
            mt_is_rlof1[seed_index],
            mt_is_rlof2[seed_index],
            mt_is_cee[seed_index],
            mt_is_mrg[seed_index])  # construct event tuple

        # If this is an entirely new seed:
        if this_seed != last_seed:                            # same seed as last seed processed?
            # no - new seed, record it
            returned_seeds.append(this_seed)
            # initialize the list of event times for this seed
            returned_times.append([this_time])
            # initialize the list of events for this seed
            returned_events.append([this_event])
            last_seed = this_seed                             # update the latest seed

        # Add event, if it is not a duplicate
        try:
            # find event_index of this particular event tuple in the array of
            # events for this seed
            event_index = returned_events[-1].index(this_event)
            # ^ if event is not a duplicate, this will throw a ValueError
            if this_time > returned_times[-1][event_index]:
                # if event is duplicate, update time to the later of the
                # duplicates
                returned_times[-1][event_index] = this_time
        except ValueError:                                    # event is not a duplicate:
            # record new event tuple for this seed
            returned_events[-1].append(this_event)
            # record new event time for this seed
            returned_times[-1].append(this_time)

    # see above for description
    return returned_seeds, returned_events, returned_times


def get_sn_data_tuple(sn_data: H5Group) -> tuple[list, list, list]:
    """Calculates the EventTuple for the BSE_Supernovae output H5Group.

    Parameters
    ----------
    sn_data : H5Group
        the COMPAS output H5Group corresponding to BSE_Supernovae

    Returns
    -------
    returned_seeds : list
        an ordered list of the unique seeds in the sn_data file
    returned_events : list
        a list of sublists, one sublist per seed, where each sublist
        contains all the SnEventTuples for the given seed
    returned_times : list
        a list of sublists of times of each of the sn_data events

    Notes
    -----
    SnEventTuples take the form: (stellar_type_progenitor, stellar_type_remnant, which_starIsProgenitor, is_binary_unbound).
    The events in the input do not have to be ordered chronologically, this function orders them, in the event
    that the input was coallated from a previous h5copy command. The output is not meant to be used generally,
    it is just a feeder for get_event_history below.
    """
    sn_seeds = sn_data['SEED'][()]
    sn_times = sn_data['Time'][()]
    sn_prog_stype = sn_data['Stellar_Type_Prev(SN)'][()]
    sn_remn_stype = sn_data['Stellar_Type(SN)'][()]
    sn_which_prog = sn_data['Supernova_State'][()]
    sn_is_unbound = sn_data['Unbound'][()] == 1

    # We want the return arrays sorted by seed, so sort here.
    # sort by seeds then times - lexsort sorts by the last column first...
    sn_seeds_idx = np.lexsort((sn_times, sn_seeds))
    sn_seeds = sn_seeds[sn_seeds_idx]
    sn_times = sn_times[sn_seeds_idx]
    sn_prog_stype = sn_prog_stype[sn_seeds_idx]
    sn_remn_stype = sn_remn_stype[sn_seeds_idx]
    sn_which_prog = sn_which_prog[sn_seeds_idx]
    sn_is_unbound = sn_is_unbound[sn_seeds_idx]

    # Process the sn_data events
    # array of seeds - will only contain seeds that have sn_data events
    returned_seeds = []
    # array of sn_data events for each seed in returned_seeds
    returned_events = []
    # array of times for each event in returned_events (for each seed in
    # returned_seeds)
    returned_times = []

    # initialize most recently processed seed
    last_seed = -1
    for seed_index, this_seed in enumerate(
            sn_seeds):          # iterate over all sn_data file entries
        # time for this sn_data file entry
        this_time = sn_times[seed_index]
        this_event = (
            sn_prog_stype[seed_index],
            sn_remn_stype[seed_index],
            sn_which_prog[seed_index],
            sn_is_unbound[seed_index])  # construct event tuple

        # If this is an entirely new seed:
        if this_seed != last_seed:                            # same seed as last seed processed?
            # no - new seed, record it
            returned_seeds.append(this_seed)
            # initialize the list of event times for this seed
            returned_times.append([this_time])
            # initialize the list of events for this seed
            returned_events.append([this_event])
            last_seed = this_seed  # update the latest seed
        else:                                               # yes - second sn_data event for this seed
            returned_times[-1].append(this_time)  # append time at end of array
            # append event at end of array
            returned_events[-1].append(this_event)

    # see above for description
    return returned_seeds, returned_events, returned_times


def get_event_history(
        data: H5File, include_null: bool = True) -> tuple[list, list]:
    """Get the event history for all seeds, including both MT and SN events, in chronological order.

    Parameters
    ----------
    data : H5File
        the COMPAS output H5file
    include_null : bool
        whether to include seeds which undergo no MT or SN events (default is True)

    Returns
    -------
    returned_seeds : list
        an ordered list of all the unique seeds in the output file
    returned_events : list
        a list where each element is a chronological set of events corresponding to the associated seed
    """
    sp_data = data['BSE_System_Parameters']
    mt_data = data['BSE_RLOF']
    sn_data = data['BSE_Supernovae']
    # get all seeds
    all_seeds = sp_data['SEED'][()]
    mt_seeds, mt_events, mt_times = get_mt_data_tuple(
        mt_data)                # get mt data tuple
    sn_seeds, sn_events, sn_times = get_sn_data_tuple(
        sn_data)                # get sn data tuple

    # number of mt_data events
    num_mt_seeds = len(mt_seeds)
    # number of sn_data events
    num_sn_seeds = len(sn_seeds)

    if num_mt_seeds < 1 and num_sn_seeds < 1:
        return []                       # no events - return empty history

    # array of seeds - will only contain seeds that have events (of any type)
    returned_seeds = []
    # array of events - same size as returned_seeds (includes event times)
    returned_events = []

    # order of preference for simultaneous events
    event_ordering = ['RL', 'CE', 'SN', 'MG']

    # index into mt_data events arrays
    mt_index = 0
    # index into sn_data events arrays
    sn_index = 0

    if include_null:
        seeds_to_iterate = all_seeds
    else:
        # iterate over all the seeds that have either mt_data or sn_data events
        seeds_to_iterate = np.sort(np.unique(np.append(mt_seeds, sn_seeds)))

    idx_ordered = np.argsort(seeds_to_iterate)
    # array of seeds - will only contain seeds that have events (of any type)
    returned_seeds = [None] * np.size(seeds_to_iterate)
    # array of events - same size as returned_seeds (includes event times)
    returned_events = [None] * np.size(seeds_to_iterate)

    for idx in idx_ordered:
        seed = seeds_to_iterate[idx]
        # initialise the events for the seed being processed
        seed_events = []

        # Collect any mt_data events for this seed, add the time of the event
        # and the event type
        while mt_index < num_mt_seeds and mt_seeds[mt_index] == seed:
            for event_index, event in enumerate(mt_events[mt_index]):
                _, _, is_rl1, is_rl2, is_cee, is_mrg = event
                # event type: Merger, CEE, RLOF 1->2, RLOF 2->1
                event_key = 'MG' if is_mrg else 'CE' if is_cee else 'RL'
                seed_events.append(
                    (event_key,
                     mt_times[mt_index][event_index],
                        *mt_events[mt_index][event_index]))
            mt_index += 1

        # Collect any sn_data events for this seed, add the time of the event
        # and the event type
        while sn_index < num_sn_seeds and sn_seeds[sn_index] == seed:
            for event_index, event in enumerate(sn_events[sn_index]):
                event_key = 'SN'
                seed_events.append(
                    (event_key,
                     sn_times[sn_index][event_index],
                        *sn_events[sn_index][event_index]))
            sn_index += 1

        # sort the events by time and event type (mt_data before sn_data if at
        # the same time)
        seed_events.sort(key=lambda ev: (ev[1], event_ordering.index(ev[0])))
        # record the seed in the seeds array being returned
        returned_seeds[idx] = seed
        # record the events for this seed in the events array being returned
        returned_events[idx] = seed_events
    return returned_seeds, returned_events


###########################################
# ## Produce strings of the event histories

simplified_stellar_type_dict = {
    0: 'MS',
    1: 'MS',
    2: 'HG',
    3: 'GB',
    4: 'GB',
    5: 'GB',
    6: 'GB',
    7: 'HE',
    8: 'HE',
    9: 'HE',
    10: 'WD',
    11: 'WD',
    12: 'WD',
    13: 'NS',
    14: 'BH',
    15: 'MR',
    16: 'MS',
}


def build_event_string(
        events: list,
        use_int_stypes: bool = False) -> EventHistoryString:
    """Produce a string representation of the event history of a single binary.

    Parameters
    ----------
    events : list
        list of mixed MtEventTuple and SnEventTuple's
    use_int_stypes : bool, optional
        whether to report stellar types as integers
        (default is False, report them as 2-char stellar types)

    Returns
    -------
    event_str : EventHistoryString
        condensed string representation of the event history of a binary, in array form

    Notes
    -----
    MT strings look like:
        P>S, P<S, or P=S where P is primary type, S is secondary type,
        and '>', '<' is RLOF (star1 -> star2 or star1 <- star2)
        or otherwise '=' for CEE or '&' for Merger.
    SN strings look like:
        P*SR if star1 is the SN progenitor, or
        R*SP if star2 is the SN progenitor,
        where P is progenitor type, R is remnant type,
        S is state ('i' for intact, 'u' for unbound)
    Event strings for the same seed are separated by an undesrcore '_'

    TODO
    ----
    Currently unvectorized, do the vectorization later for a speed up
    """
    def _remap_stype(int_stype):
        if use_int_stypes:
            return str(int_stype)
        else:
            return simplified_stellar_type_dict[int_stype]

    # Empty event
    if len(events) == 0:
        return 'NA'

    event_str = ''
    for event in events:
        if event[0] == 'SN':                                                # SN event
            _, time, stype_p, stype_c, which_sn, is_unbound = event
            # Progenitor or Remnant depending upon which star is the SN
            if which_sn == 1:
                char_l = _remap_stype(stype_p)
                char_r = _remap_stype(stype_c)
            else:
                char_l = _remap_stype(stype_c)
                char_r = _remap_stype(stype_p)
            char_m = '*u' if is_unbound else '*i'                           # unbound or intact
        else:                                                               # Any of the MT events
            _, time, stype_1, stype_2, is_rl1, is_rl2, is_cee, is_mrg = event
            # primary stellar type
            char_l = _remap_stype(stype_1)
            # secondary stellar type
            char_r = _remap_stype(stype_2)
            # event type: CEE, RLOF 2->1, RLOF 1->2
            char_m = '&' if is_mrg else '=' if is_cee else '<' if is_rl2 else '>'
        # event string for this star, _ is event separator
        event_str += "{}{}{}_".format(char_l, char_m, char_r)
    # return event string for this star (pop the last underscore first)
    event_str = np.array(event_str[:-1], dtype=np.str_)
    return event_str


def get_event_strings(
        data: H5File = None,
        all_events: list = None,
        use_int_stypes: bool = False):
    """Calculate the event history strings for a COMPAS population.

    Parameters
    ----------
    data : H5File, optional
        the COMPAS output H5file
    all_events : list, optional
        a list where each element is a chronological set of events corresponding to the associated seed,
        one of the outputs of get_event_history(data)
    use_int_stypes : bool, optional

    Returns
    -------
    event_strings : array of EventHistoryString's
        strings of the event history of each system

    Notes
    -----
    Exactly one of data or all_events must be included. If neither, nothing is returned.
    """

    # If output is
    if (data is None) & (all_events is None):
        return
    elif (all_events is None):
        _, all_events = get_event_history(data)

    # Keep the full string values instead of truncating to a single character.
    event_strings = np.empty(len(all_events), dtype=object)
    for ii, events_for_given_seed in enumerate(all_events):
        event_string = build_event_string(
            events_for_given_seed,
            use_int_stypes=use_int_stypes)
        # append event string for this star (pop the last underscore first)
        event_strings[ii] = str(event_string)
    return event_strings

def main():
    return

if __name__ == "__main__":
    main()
