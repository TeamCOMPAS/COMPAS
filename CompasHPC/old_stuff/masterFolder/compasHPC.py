#!/usr/bin/env python

import os

#-- Change directory to where COMPAS was run from
path = os.path.dirname(os.path.abspath(__file__))
os.chdir(path)

import pythonSubmit as ps 

#-- Get the program options
programOptions = ps.pythonProgramOptions()

ps.runCompas(programOptions)