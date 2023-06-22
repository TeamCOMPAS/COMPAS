#!/usr/bin/env python3
import numpy as np
import h5py  as h5
import os

import scipy.special
from   scipy.optimize import newton
from   astropy.cosmology import FlatLambdaCDM
from   astropy.cosmology import WMAP9 #as cosmo
import astropy.units as u

from . import ClassCOMPAS
from . import ClassMSSFR
import selection_effects
import totalMassEvolvedPerZ as MPZ
class CosmicIntegrator(object):
    """

    The cosmological integrator calculates the rate 
    an object given the SFR that went into it at birth

    The entire class consists of several subclasses and their instances
    are done in this order (since they are dependent of each other)


    -Class cosmo: Astropy class from python library to set the cosmological relations
                  between luminosity distance, redshift, age univere etc. We define it once
                  here so we can consistently pass around in other classes.


    -Class MergerTimes: Creates the numpy arrays for the shells at which we calculate
                        The mergers. In turn it also creates the lookback times for the systems
                        of the Data class

    -Class MSSFR: Class that calculates the metallicity Specific StarFormation rate
                  for each individual object at each individual merger time of interest

    Each class lives in a separate file (except cosmo which is inbuilt library astropy)

    
    """

    def __init__(self, pathCOMPAS=None, Cosmology='WMAP',hubbleConstant = 67.8,\
                omegaMatter=0.308,redshiftFirstSFR=10., \
                minRedshift=0.0,   maxRedshift=2., nrRedshiftBins=20,\
                RedshiftTabulated =True, RedshiftTabulatedResolution=100000,
                GWdetector_sensitivity='O1', GWdetector_snrThreshold=8, verbose = False):

        #################################################
        #                                               #
        #     initialize universe and shell integral    #
        #                                               #
        #################################################

        #Define topology universe using astropy
        self.verbose                  = verbose
        self.pathCOMPAS               = pathCOMPAS

        if Cosmology == 'WMAP':
            self.cosmology            = WMAP9
        if Cosmology == 'Custom Flat':
            self.cosmology            = FlatLambdaCDM(H0=hubbleConstant *\
                                        u.km / u.s / u.Mpc, Om0=omegaMatter)
        self.redshiftFirstSFR         = redshiftFirstSFR
        self.ageFirstSFR              = self.cosmology.age(self.redshiftFirstSFR).value
        
        #These are the redshifts shells we integrate over
        self.minRedshift              = minRedshift
        self.maxRedshift              = maxRedshift
        self.nrRedshiftBins           = nrRedshiftBins
        #These are set by function
        self.Shell_centerRedshift     = None
        self.Shell_volume             = None
        self.Shell_dz                 = None
        self.Shell_luminosityDistance = None
        self.createConcentricRedshiftShells()


        ################################################################
        #     initialize classCOMPAS popsynth                          #
        #     uses dummy variables which dont work such that           #
        #     user HAS TO define it and check their input              #
        ################################################################
        print("reminders of what to set in the following order:")
        print()

        if self.verbose:
            print("Creating instance COMPAS class User has to still set DCO and Data")
        #setting Mlower/Mupper etc to None to force warning for user
        self.COMPAS  = ClassCOMPAS.COMPASData(path=self.pathCOMPAS, lazyData=True,\
                                          Mlower=None, Mupper=None, \
                                          binaryFraction=None)
        #if pathCOMPAS is None:
        #    # ClassCOMPAS will assume default pathCOMPAS "COMPAS_Output.h5"
        #    self.COMPAS  = ClassCOMPAS.COMPASData(path=self.pathCOMPAS, lazyData=True,\
        #                                      Mlower=None, Mupper=None, \
        #                                      binaryFraction=None)
        #else:
        #    self.COMPAS  = ClassCOMPAS.COMPASData(path=self.pathCOMPAS, lazyData=True,\
        #                                      Mlower=None, Mupper=None, \
        #                                     binaryFraction=None)

        #####################################################
        #     set the MSSFR class                           #
        #####################################################
        if self.verbose:
            print("Creating instance MSSFR class User has to still set grid")
        self.MSSFR = ClassMSSFR.MSSFR(metallicityGrid=None, cosmo=self.cosmology)       

        #additionally needed for sensitivity
        self.GWdetector_sensitivity     = GWdetector_sensitivity
        self.GWdetector_snrThreshold    = GWdetector_snrThreshold

        print("ClassCosmicIntegrator: Remember to setBirthTimesAnd2Darrays()")
        print("                        to prepare calculation/results")
        print("                        do you have the GW-detector you want?")

        #################################################
        #                                               #
        #     the eventual results of interest          #
        #     if we already have data we can set it     #
        #     else we have to manually call tis function#
        #################################################
        #Calculating birth age is simple subtraction age-merger -delaytime
        self.PerSystemPerRedshift_ageBirth       = None
        #Cannot do this in redshift so convert table above to a redshift table
        self.PerSystemPerRedshift_redshiftBirth  = None
        # dN Gpc-3 per year at z=redshift.     #each row is redshift merger
        self.PerSystemPerRedshift_ratesIntrinsic = None
        # dN per year in detector from shell 
        self.PerSystemPerRedshift_ratesObserved  = None

        self.RedshiftTabulated                   = RedshiftTabulated
        self.RedshiftTabulatedResolution         = RedshiftTabulatedResolution
        self.redshiftAgeTable                    = None 
        

    def createConcentricRedshiftShells(self):
        if self.verbose:
            print("Creating redshift shells for integral")

        #Thanx Jim Barrett for cleaning up this part a bit :D
        #Flat universe with Hubble constant of 70 and OmegaM of 0.3
        redshiftEdges = np.linspace(self.minRedshift,\
                                    self.maxRedshift,\
                                    self.nrRedshiftBins+1) #The bin edges in redshift
        #Central value of each redshift bin, this is the value used in cosmic Int
        self.Shell_dz             = np.diff(redshiftEdges)
        self.Shell_centerRedshift = 0.5*(redshiftEdges[:-1] + redshiftEdges[1:])      
        #Corresponding luminosity Distances [Mpc]
        self.Shell_luminosityDistance = \
        [x.value for x in self.cosmology.luminosity_distance(self.Shell_centerRedshift)] 

        #cosmology.comoving volume is in units of Mpc^-3
        #The line below translates this to Gpc^-3 (times 1e-9) 
        #and gives me the spehrical volume of each redshift bin.
        comovEdges    = [x.value*1e-9 for x in self.cosmology.comoving_volume(redshiftEdges)]
        #The difference between the spherical volumes gives me 
        #the shell volume of each redshift we are looking at
        self.Shell_volume   = np.diff(comovEdges)



    def setAgeBirthSystems(self):
        #each row is a redshift each column a system
        self.PerSystemPerRedshift_ageBirth       = np.zeros(shape=(int(self.nrRedshiftBins),\
                                                            len(self.COMPAS.delayTimes)))

        #the age of the universe at the birth of the system is the 
        #age of the universe at merger minus the delay time.
        #for the age of the universe in a redshift shell 
        #we use the center redshift value
        for nrz, redshift in enumerate(self.Shell_centerRedshift):
            ageUniverseAtMergerGyr = self.cosmology.age(redshift)
            delayTimeGyr           = self.COMPAS.delayTimes / 1000.
            ageBirth               = ageUniverseAtMergerGyr.value - delayTimeGyr
            maskUnreal             = ageBirth < self.ageFirstSFR
            self.PerSystemPerRedshift_ageBirth[nrz] = ageBirth
            self.PerSystemPerRedshift_ageBirth[nrz][maskUnreal] = -1


    def setRedshiftBirthSystems(self):
        self.PerSystemPerRedshift_redshiftBirth = np.zeros(shape=(int(self.nrRedshiftBins),\
                                                            len(self.COMPAS.delayTimes)))

        #calculating redshift from age is cheap, inverse not
        #either we precalulate a grid and look up closest value
        #cheap and quick but technically less precise
        #or we use an inverse function, which is slower
        if (self.RedshiftTabulated) and (self.redshiftAgeTable is None):
            redshifts  =  np.linspace(0.0000001,\
                                      self.redshiftFirstSFR,\
                                      self.RedshiftTabulatedResolution)
            ages       = self.cosmology.age(redshifts).value
            self.redshiftAgeTable = np.column_stack((redshifts,ages))

        for nrR, row in enumerate(self.PerSystemPerRedshift_ageBirth):
            maskUnphysical = (row == -1) #born before first SFR
            maskPhysical   = np.logical_not(maskUnphysical)
            
            if self.RedshiftTabulated:
                #find indices in precalculated table 
                inds      = np.digitize(row[maskPhysical], self.redshiftAgeTable[:,1])
                redshifts = self.redshiftAgeTable[inds,0]
            else:
                """
                All credits go to
                https://mail.scipy.org/pipermail/astropy/2013-December/002950.html
                The problem is that astropy seems to only have redshift->time
                not time->redshift
                """
                redshifts = np.zeros(np.sum(maskPhysical))
                for nra, age in enumerate(maskPhysical):
                     redshifts[nra] = newton(lambda x: \
                                    self.cosmology.age(x).value-age, 0)
            #fill in values
            self.PerSystemPerRedshift_redshiftBirth[nrR][maskUnphysical] = -1
            self.PerSystemPerRedshift_redshiftBirth[nrR][maskPhysical]   = redshifts

    def setBirthTimesAnd2Darrays(self):            

        if self.COMPAS.delayTimes is not None:
            if self.verbose:
                print("creating 2D arrays with birth ages and redshift")
               
            
            #from the COMPAS data we need the birth ages of each system
            #in our integral, only if we actually have data. 
            #Introduced this because you might want an empty class (without data)
            #for testing and explaining code in notebooks :D

            #Calculating birth age is simple subtraction age-merger -delaytime
            self.setAgeBirthSystems()
            #Cannot do this in redshift so convert table above to a redshift table
            self.setRedshiftBirthSystems()
            if self.verbose:
                print("creating 2D array with intrinsic and observed rate to be calculated")
             # dN Gpc-3 per year at z=redshift.     #each row is redshift merger
            self.PerSystemPerRedshift_ratesIntrinsic = np.zeros(shape=(int(self.nrRedshiftBins),\
                                                                len(self.COMPAS.delayTimes)))
            # dN per year in detector from shell 
            self.PerSystemPerRedshift_ratesObserved  = np.zeros(shape=(int(self.nrRedshiftBins),\
                                                                len(self.COMPAS.delayTimes)))    
        else:
            print()
            print("cannot set 2D-array of rates")
            print("COMPAS data is empty (COMPAS.setCOMPASData) " )

    def cosmologicalIntegration(self):
        if self.verbose:
            print("Doing the actual cosmic integration")
            print("Filling in the 2D arrays with rate per system per redshift")
        #For each row in 2D array which corresponds to a merger redshift
        # Get birth Age in Gyr , redshifts birth, and metallicities from COMPAS
        #Calculate MSSFR for that row and fill in the answer
        for nr in range(int(self.nrRedshiftBins)):
            for nrZ, Z in enumerate(self.COMPAS.metallicityGrid):
                maskZ    = self.COMPAS.metallicitySystems == Z
                MSSFR  = self.MSSFR.returnMSSFR(metallicity=Z,\
                                                  agesBirth=self.PerSystemPerRedshift_ageBirth[nr][maskZ],
                                                  redshiftBirth=self.PerSystemPerRedshift_redshiftBirth[nr][maskZ])
                RatesZ   = np.divide(MSSFR, self.COMPAS.totalMassEvolvedPerZ[nrZ])
                self.PerSystemPerRedshift_ratesIntrinsic[nr][maskZ] = RatesZ   # intrinric dN Gpc-3 per year at z=redshift.
                probObservingZ     = selection_effects.detection_probability(\
                                    self.COMPAS.mass1[maskZ],self.COMPAS.mass2[maskZ], self.Shell_centerRedshift[nr],\
                                    self.Shell_luminosityDistance[nr], self.GWdetector_snrThreshold,\
                                    sensitivity=self.GWdetector_sensitivity)
                NrMergersInShell  = np.multiply(RatesZ,self.Shell_volume[nr])                    # intrinsic rate per system per shell
                NrMergersInShell  = NrMergersInShell * (1./(1.+self.Shell_centerRedshift[nr])) # rate observer frame per system per shell
                self.PerSystemPerRedshift_ratesObserved[nr][maskZ]  = np.multiply(NrMergersInShell,probObservingZ)      # observed  dN per year prob between 0-1 
        if np.sum(self.PerSystemPerRedshift_ratesObserved[-1]) != 0 :
            print("The detected rate of the outermost redshift shell is nonzero, did we integrate far enough?")

