import numpy as np
import subprocess
import sys
import os
import pickle
import itertools 
from subprocess import call
import yaml

# Check if we are using python 3
python_version = sys.version_info[0]
print("python_version =", python_version)

class pythonProgramOptions:
    """
    A class to store and access COMPAS program options in python
    """

    def __init__(self, config_file='compasConfigDefault.yaml', grid_filename=None, random_seed_filename='randomSeed.txt', output_directory=None):
        # Do './COMPAS --help' to see all options
        #-- Define variables

        # environment variable COMPAS_EXECUTABLE_PATH is used for docker runs
        # if COMPAS_EXECUTABLE_PATH is not set (== None) we assume this is an
        # interactive run with python3
        # if COMPAS_EXECUTABLE_PATH is set (!= None) we assume this is a run
        # inside a docker container - we have different directories inside a 
        # docker container (src, obj, bin), and the COMPAS executable resides
        # in the bin directory (rather than the src directory)
        
        # Load yaml file with options
        with open(config_file) as file:
            # The FullLoader parameter handles the conversion from YAML
            # scalar values to Python the dictionary format
            config = yaml.load(file, Loader=yaml.FullLoader)
        
        self.booleanChoices = config['booleanChoices']
        self.numericalChoices = config['numericalChoices']
        self.stringChoices = config['stringChoices']
        self.listChoices = config['listChoices']
                
        if self.booleanChoices   is None: self.booleanChoices   = {}
        if self.numericalChoices is None: self.numericalChoices = {}
        if self.stringChoices    is None: self.stringChoices    = {}
        if self.listChoices      is None: self.listChoices      = {}

        compas_executable_override = os.environ.get('COMPAS_EXECUTABLE_PATH')
        print('compas_executable_override', compas_executable_override)        

        if (compas_executable_override is None):
            # Make sure you have set and exported the COMPAS_ROOT_DIR environment variable
            git_directory = os.environ.get('COMPAS_ROOT_DIR')
            self.compas_executable = os.path.join(git_directory, 'src/COMPAS')
        else:
            self.compas_executable = compas_executable_override

        # If random_seed_filename is specified, overwrite the random seed from the yaml file
        if os.path.isfile(random_seed_filename): 
            self.numericalChoices['--random-seed'] = int(np.loadtxt(random_seed_filename))

        # If grid is specified in pythonProgramOptions(), ignore the values from yaml file
        if grid_filename:
            self.grid_filename = grid_filename
            self.stringChoices['--grid']= self.grid_filename
        else:
            if not self.stringChoices:
                self.grid_filename = None
            else:
                self.grid_filename = None
                if '--grid' in self.stringChoices:
                    self.grid_filename = self.stringChoices['--grid']

        print('grid_filename', self.grid_filename)

        # environment variable COMPAS_LOGS_OUTPUT_DIR_PATH is used primarily for docker runs
        # if COMPAS_LOGS_OUTPUT_DIR_PATH is set (!= None) it is used as the value for the
        # --output-path option
        # if COMPAS_LOGS_OUTPUT_DIR_PATH is not set (== None) the current working directory
        # is used as the value for the --output-path option
        compas_logs_output_override = os.environ.get('COMPAS_LOGS_OUTPUT_DIR_PATH')
        
        if (compas_logs_output_override is None):
            self.stringChoices['--output-path'] = os.getcwd()
            self.stringChoices['--output-container'] = output_directory   # names the directory to be created and in which log files are created.  Default in COMPAS is "COMPAS_Output"
        else:
            self.stringChoices['--output-path'] = compas_logs_output_override
            self.stringChoices['--output-container'] = output_directory 
        
        
        # environment variable COMPAS_INPUT_DIR_PATH is used primarily for docker runs
        # if COMPAS_INPUT_DIR_PATH is set (!= None) it is prepended to input filenames
        # (such as grid_filename and logfile_definitions)
        # if COMPAS_INPUT_DIR_PATH is not set (== None) the current working directory
        # is prepended to input filenames
        compas_input_path_override = os.environ.get('COMPAS_INPUT_DIR_PATH')
        
        if self.grid_filename != None:
            if compas_input_path_override == None:
                self.grid_filename = os.getcwd() + '/' + self.grid_filename
            else:
                self.grid_filename = compas_input_path_override + '/' + self.grid_filename

        if '--logfile-definitions' in self.stringChoices:
            self.logfile_definitions = self.stringChoices['--logfile-definitions']  # logfile record definitions file name (e.g. 'logdefs.txt')
        else:
            self.logfile_definitions = None

        if self.logfile_definitions != None:
            if compas_input_path_override == None:
                self.logfile_definitions = os.getcwd() + '/' + self.logfile_definitions
            else:
                self.logfile_definitions = compas_input_path_override + '/' + self.logfile_definitions

        self.makeCommandString()
                        
    def makeCommandString(self):
        """
        This function generates a dictionary mapping COMPAS options to their specified 
        values (or empty strings for boolean options). These are then combined into a string
        that can be run directly as a terminal command, or passed to the stroopwafel interface
        where some of them may be overwritten. Options not to be included in the command 
        line should be set to pythons None (except booleans, which should be set to False)
        """

        ### Collect all options into a dictionary mapping option name to option value
        self.command = {'compas_executable' : self.compas_executable} 

        for boolKey, boolVal in self.booleanChoices.items():
            if boolVal is True:
                self.command.update({boolKey: ''}) 
            else:
                self.command.update({boolKey: 'FALSE'}) 
                
        for numKey, numVal in self.numericalChoices.items():
            if not numVal == None:
                self.command.update({numKey : str(numVal)})    
                
        for strKey, strVal in self.stringChoices.items():
            if not strVal == None:
                self.command.update({strKey: strVal})
                
        for listKey, listVal in self.listChoices.items():
            if listVal:
                self.command.update({listKey : ' '.join(map(str, listVal))})
        
        # Ensure the Compas executable is first, and not repeated.
        # Options are non-ordered.
        self.shellCommand = self.command['compas_executable']
        del self.command['compas_executable'] 
        for key, val in self.command.items():
            self.shellCommand += ' ' + key + ' ' + val

        return

if __name__ == "__main__":
    #-- Get the program options
    myoptions = pythonProgramOptions(config_file='compasConfigDemo.yaml')
    print(myoptions.shellCommand)
    with open('runSubOptions.txt', 'w') as fwrite:
        fwrite.write(myoptions.shellCommand)
    #-- Run execute COMPAS shell string
    call(myoptions.shellCommand,shell=True)
