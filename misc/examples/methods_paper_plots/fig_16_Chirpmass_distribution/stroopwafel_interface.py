#!/usr/bin/env python
import h5py as h5
import os, sys
import pandas as pd
import shutil
import time
import numpy as np
sys.path.append('/n/home04/lvanson/Programs/stroopwafel/') # Specific for Lieke running on cannon :D
from stroopwafel import sw, classes, prior, sampler, distributions, constants, utils
import argparse

# TODO fix issues with adaptive sampling
# TODO add in functionality for alternative pythonSubmit names and locations

#######################################################
### 
### For User Instructions, see 'docs/sampling.md'
### 
#######################################################


### Include options from local pythonSubmit file      
usePythonSubmit = True #If false, use stroopwafel defaults

### Set default stroopwafel inputs - these are overwritten by any command-line arguments

compas_executable = os.path.join(os.environ.get('COMPAS_ROOT_DIR'), 'src/COMPAS')   # Location of the executable      # Note: overrides pythonSubmit value
num_systems = int(1e7)                    # Number of binary systems to evolve                                              # Note: overrides pythonSubmit value
output_folder = '../../output'           # Location of output folder (relative to cwd)                                     # Note: overrides pythonSubmit value
random_seed_base = 0                # The initial random seed to increment from                                       # Note: overrides pythonSubmit value

num_cores = 50                        # Number of cores to parallelize over 
num_per_core = int(1e4)               # Number of binaries per batch
mc_only = False                      # Exclude adaptive importance sampling (currently not implemented, leave set to True)
run_on_hpc = True                    # Run on slurm based cluster HPC

output_filename = 'samples.csv'     # output filename for the stroopwafel samples
debug = True                        # show COMPAS output/errors
hdf5 = True


### Default options for interesting systems when using AIS: ['BBH', 'DNS', 'BHNS']
sys_int = 'BBH'

def create_dimensions():
    """
    This Function that will create all the dimensions for stroopwafel, a dimension is basically one of the variables you want to sample
    Invoke the Dimension class to create objects for each variable. Look at the Dimension class definition in classes.py for more.
    It takes the name of the dimension, its max and min value. 
    The Sampler class will tell how to sample this dimension. Similarly, prior tells it how it calculates the prior. You can find more of these in their respective modules
    OUT:
        As Output, this should return a list containing all the instances of Dimension class.
    """
    m1 = classes.Dimension('--initial-mass-1', 10, 150, sampler.kroupa, prior.kroupa)
    q = classes.Dimension('q', 0.01, 1, sampler.uniform, prior.uniform, should_print = False)
    a = classes.Dimension('--semi-major-axis', .01, 1000, sampler.flat_in_log, prior.flat_in_log)
    #kick_magnitude_1 = classes.Dimension('--kick-magnitude-1', 0, 500, sampler.uniform, prior.uniform)
    #kick_theta_1 = classes.Dimension('--kick-theta-1', -np.pi / 2, np.pi / 2, sampler.uniform_in_cosine, prior.uniform_in_cosine)
    #kick_phi_1 = classes.Dimension('--kick-phi-1', 0, 2 * np.pi, sampler.uniform, prior.uniform)
    #kick_magnitude_2 = classes.Dimension('--kick-magnitude-2', 0, 500, sampler.uniform, prior.uniform)
    #kick_theta_2 = classes.Dimension('--kick-theta-2', -np.pi / 2, np.pi / 2, sampler.uniform_in_cosine, prior.uniform_in_cosine)
    #kick_phi_2 = classes.Dimension('--kick-phi-2', 0, 2 * np.pi, sampler.uniform, prior.uniform)
    return [m1, q, a]#, kick_magnitude_1, kick_theta_1, kick_phi_1, kick_magnitude_2, kick_theta_2, kick_phi_2]

def update_properties(locations, dimensions):
    """
    This function is not mandatory, it is required only if you have some dependent variable. 
    For example, if you want to sample Mass_1 and q, then Mass_2 is a dependent variable which is product of the two.
    Similarly, you can assume that Metallicity_2 will always be equal to Metallicity_1
    IN:
        locations (list(Location)) : A list containing objects of Location class in classes.py. 
        You can play with them and update whatever fields you like or add more in the property (which is a dictionary)
    OUT: Not Required
    """
    m1 = dimensions[0]
    q = dimensions[1]
    for location in locations:
        location.properties['--initial-mass-2'] = location.dimensions[m1] * location.dimensions[q]
        location.properties['--metallicity'] =10**(np.random.uniform(-4, np.log10(0.03)))
        location.properties['--eccentricity'] = 0
        #location.properties['Kick_Mean_Anomaly_1'] = np.random.uniform(0, 2 * np.pi, 1)[0]
        #location.properties['Kick_Mean_Anomaly_2'] = np.random.uniform(0, 2 * np.pi, 1)[0]





#################################################################################
#################################################################################
###                                                                           ###
###         USER SHOULD NOT SET ANYTHING BELOW THIS LINE                      ###
###                                                                           ###
#################################################################################
#################################################################################





def configure_code_run(batch):
    """
    This function tells stroopwafel what program to run, along with its arguments.
    IN:
        batch(dict): This is a dictionary which stores some information about one of the runs. It has an number key which stores the unique id of the run
            It also has a subprocess which will run under the key process. Rest, it depends on the user. User is free to store any information they might need later 
            for each batch run in this dictionary. For example, here I have stored the 'output_container' and 'grid_filename' so that I can read them during discovery of interesting systems below
    OUT:
        compas_args (list(String)) : This defines what will run. It should point to the executable file along with the arguments.
        Additionally one must also store the grid_filename in the batch so that the grid file is created
    """
    batch_num = batch['number']
    grid_filename = os.path.join(output_folder, 'grid_' + str(batch_num) + '.csv')
    output_container = 'batch_' + str(batch_num)
    random_seed = random_seed_base + batch_num*NUM_SYSTEMS_PER_RUN  # ensure that random numbers are not reused across batches
    compas_args = [compas_executable, '--grid', grid_filename, '--output-container', output_container, '--random-seed' , random_seed]
    [compas_args.extend([key, val]) for key, val in commandOptions.items()]
    for params in extra_params:
        compas_args.extend(params.split("="))
    batch['grid_filename'] = grid_filename
    batch['output_container'] = output_container
    return compas_args

def interesting_systems(batch):
    """
    This is a mandatory function, it tells stroopwafel what an interesting system is. User is free to define whatever looks interesting to them.
    IN:
        batch (dict): As input you will be given the current batch which just finished its execution. You can take in all the keys you defined in the configure_code_run method above
    OUT:
        Number of interesting systems
        In the below example, I define all the NSs as interesting, so I read the files, get the SEED from the system_params file and define the key is_hit in the end for all interesting systems 
    """
    try:
        folder = os.path.join(output_folder, batch['output_container'])
        #shutil.move(batch['grid_filename'], folder + '/grid_' + str(batch['number']) + '.csv')
        # First try to find the csv file
        if not hdf5:
            system_parameters = pd.read_csv(folder + '/BSE_System_Parameters.csv', skiprows = 2)
            system_parameters.rename(columns = lambda x: x.strip(), inplace = True)
            seeds = system_parameters['SEED']
            double_compact_objects = pd.read_csv(folder + '/BSE_Double_Compact_Objects.csv', skiprows = 2)
            double_compact_objects.rename(columns = lambda x: x.strip(), inplace = True)
        else:
            # Check for hdf5 file
            sfile = h5.File(folder + '/batch_'+ str(batch['number']) +'.h5' ,'r')
            seeds = sfile['BSE_System_Parameters']['SEED'][:]
        for index, sample in enumerate(batch['samples']):
            seed = seeds[index]
            sample.properties['SEED'] = seed
            sample.properties['is_hit'] = 0
            sample.properties['batch'] = batch['number']

        if hdf5:
            double_compact_objects = sfile['BSE_Double_Compact_Objects']

        st1 = double_compact_objects['Stellar_Type(1)'][:]
        st2 = double_compact_objects['Stellar_Type(2)'][:]
        merger_flag = double_compact_objects['Merges_Hubble_Time'][:]    
        dco_seeds = double_compact_objects['SEED'][:]

        #Generally, this is the line you would want to change.
        if sys_int == 'BBH':
            dco_mask = np.logical_and(st1 == 14, st2 == 14)
        if sys_int == 'DNS':
            dco_mask = np.logical_and(st1 == 13, st2 == 13)
        if sys_int == 'BHNS':
            dco_mask = np.logical_and(st1 == 14, st2 == 13) | np.logical_and(st1 == 13, st2 == 14)
        merge_mask = merger_flag == 1
        dns_mask = np.logical_and(merge_mask, dco_mask)

        # Lieke: select systems of interest
        interesting_systems_seeds = set(dco_seeds[dns_mask])
        for index, sample in enumerate(batch['samples']):
            if sample.properties['SEED'] in interesting_systems_seeds:
                sample.properties['is_hit'] = 1

        # If you were working with an hdf5 file, make sure to close it again
        if hdf5:
            sfile.close()

        return sum(dns_mask) #len(dns)

    # You probably had no DCO's in your batch
    except IOError as error:
        print('You ran into an error, ', error)
        return 0

def selection_effects(sw):
    """
    This is not a mandatory function, it was written to support selection effects
    Fills in selection effects for each of the distributions
    IN:
        sw (Stroopwafel) : Stroopwafel object
    """
    if hasattr(sw, 'adapted_distributions'):
        biased_masses = []
        rows = []
        for distribution in sw.adapted_distributions:
            folder = os.path.join(output_folder, 'batch_' + str(int(distribution.mean.properties['batch'])))
            try:
                dco_file = pd.read_csv(folder + '/BSE_Double_Compact_Objects.csv', skiprows = 2)
                dco_file.rename(columns = lambda x: x.strip(), inplace = True)
            except:
                sfile = h5.File(folder + '/batch_'+ str(batch['number']) +'.h5' ,'r')
                dco_file = sfile['BSE_Double_Compact_Objects']
                sfile.close()

            row = dco_file.loc[dco_file['SEED'][:] == distribution.mean.properties['SEED']]
            rows.append([row.iloc[0]['Mass(1)'], row.iloc[0]['Mass(2)']])
            biased_masses.append(np.power(max(rows[-1]), 2.2))
        # update the weights
        mean = np.mean(biased_masses)
        for index, distribution in enumerate(sw.adapted_distributions):
            distribution.biased_weight = np.power(max(rows[index]), 2.2) / mean

def rejected_systems(locations, dimensions):
    """
    This method takes a list of locations and marks the systems which can be
    rejected by the prior distribution
    IN:
        locations (List(Location)): list of location to inspect and mark rejected
    OUT:
        num_rejected (int): number of systems which can be rejected
    """
    m1 = dimensions[0]
    q = dimensions[1]
    a = dimensions[2]
    mass_1 = [location.dimensions[m1] for location in locations]
    mass_2 = [location.properties['--initial-mass-2'] for location in locations]
    try:
        metallicity = [location.properties['--metallicity'] for location in locations]
    except:
        Z = dimensions[3]
        metallicity = [location.dimensions[Z] for location in locations]
    eccentricity = [location.properties['--eccentricity'] for location in locations]
    num_rejected = 0
    for index, location in enumerate(locations):
        radius_1 = utils.get_zams_radius(mass_1[index], metallicity[index])
        radius_2 = utils.get_zams_radius(mass_2[index], metallicity[index])
        star_to_roche_lobe_radius_ratio_1 = radius_1 / (location.dimensions[a] * (1 - eccentricity[index]) * utils.calculate_roche_lobe_radius(mass_1[index], mass_2[index]))
        star_to_roche_lobe_radius_ratio_2 = radius_2 / (location.dimensions[a] * (1 - eccentricity[index]) * utils.calculate_roche_lobe_radius(mass_2[index], mass_1[index]))
        location.properties['is_rejected'] = 0
        if (mass_2[index] < constants.MINIMUM_SECONDARY_MASS) or (location.dimensions[a] <= (radius_1 + radius_2)) \
        or star_to_roche_lobe_radius_ratio_1 > 1 or star_to_roche_lobe_radius_ratio_2 > 1:
            location.properties['is_rejected'] = 1
            num_rejected += 1
    return num_rejected


if __name__ == '__main__':

    # STEP 1 : Import and assign input parameters for stroopwafel 
    parser=argparse.ArgumentParser()
    parser.add_argument('--num_systems', help = 'Total number of systems', type = int, default = num_systems)  
    parser.add_argument('--num_cores', help = 'Number of cores to run in parallel', type = int, default = num_cores)
    parser.add_argument('--num_per_core', help = 'Number of systems to generate in one core', type = int, default = num_per_core)
    parser.add_argument('--debug', help = 'If debug of COMPAS is to be printed', type = bool, default = debug)
    parser.add_argument('--mc_only', help = 'If run in MC simulation mode only', type = bool, default = mc_only)
    parser.add_argument('--run_on_hpc', help = 'If we are running on a (slurm-based) HPC', type = bool, default = run_on_hpc)
    parser.add_argument('--output_filename', help = 'Output filename', default = output_filename)
    parser.add_argument('--output_folder', help = 'Output folder name', default = output_folder)
    namespace, extra_params = parser.parse_known_args()

    start_time = time.time()
    #Define the parameters to the constructor of stroopwafel
    TOTAL_NUM_SYSTEMS = namespace.num_systems #total number of systems you want in the end
    NUM_CPU_CORES = namespace.num_cores #Number of cpu cores you want to run in parallel
    NUM_SYSTEMS_PER_RUN = namespace.num_per_core #Number of systems generated by each of run on each cpu core
    debug = namespace.debug #If True, will print the logs given by the external program (like COMPAS)
    run_on_hpc = namespace.run_on_hpc #If True, it will run on a clustered system helios, rather than your pc
    mc_only = namespace.mc_only # If you dont want to do the refinement phase and just do random mc exploration
    output_filename = namespace.output_filename #The name of the output file
    output_folder = os.path.join(os.getcwd(), namespace.output_folder)

    # Set commandOptions defaults - these are Compas option arguments
    commandOptions = dict()
    commandOptions.update({'--output-path' : output_folder}) 

    # Over-ride with pythonSubmit parameters, if desired
    if usePythonSubmit:
        try:
            from pythonSubmit import pythonProgramOptions
            programOptions = pythonProgramOptions()
            pySubOptions = programOptions.generateCommandLineOptionsDict()

            # Remove extraneous options
            pySubOptions.pop('compas_executable', None)
            pySubOptions.pop('--grid', None)
            pySubOptions.pop('--output-container', None)
            pySubOptions.pop('--number-of-binaries', None)
            pySubOptions.pop('--output-path', None)
            pySubOptions.pop('--random-seed', None)

            commandOptions.update(pySubOptions)

        except:
            print("Invalid pythonSubmit file, using default stroopwafel options")
            usePythonSubmit = False
    

    print("Output folder is: ", output_folder)
    if os.path.exists(output_folder):
        #command = input ("The output folder already exists. If you continue, I will remove all its content. Press (Y/N)\n")
        #if (command == 'Y'):
        shutil.rmtree(output_folder)
        #else:
        #    exit()
    os.makedirs(output_folder)


    # STEP 2 : Create an instance of the Stroopwafel class
    sw_object = sw.Stroopwafel(TOTAL_NUM_SYSTEMS, NUM_CPU_CORES, NUM_SYSTEMS_PER_RUN, output_folder, output_filename, debug = debug, run_on_helios = run_on_hpc, mc_only = mc_only)


    # STEP 3: Initialize the stroopwafel object with the user defined functions and create dimensions and initial distribution
    dimensions = create_dimensions()
    sw_object.initialize(dimensions, interesting_systems, configure_code_run, rejected_systems, update_properties_method = update_properties)
    #sw_object.initialize(dimensions, None, configure_code_run, None, update_properties_method = update_properties)


    intial_pdf = distributions.InitialDistribution(dimensions)
    # STEP 4: Run the 4 phases of stroopwafel
    sw_object.explore(intial_pdf) #Pass in the initial distribution for exploration phase
    sw_object.adapt(n_dimensional_distribution_type = distributions.Gaussian) #Adaptation phase, tell stroopwafel what kind of distribution you would like to create instrumental distributions
    ## Do selection effects
    #selection_effects(sw)
    sw_object.refine() #Stroopwafel will draw samples from the adapted distributions
    sw_object.postprocess(distributions.Gaussian, only_hits = False) #Run it to create weights, if you want only hits in the output, then make only_hits = True

    end_time = time.time()
    print ("Total running time = %d seconds" %(end_time - start_time))

