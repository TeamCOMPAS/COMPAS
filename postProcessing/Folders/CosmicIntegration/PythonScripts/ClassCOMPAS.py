#!/usr/bin/env python3
import numpy as np
import h5py  as h5
import os
import totalMassEvolvedPerZ as MPZ

class COMPASData(object):


    def __init__(self, path=None, fileName='COMPAS_output.h5',\
                 lazyData=True, Mlower=None, \
                 Mupper=None, binaryFraction=None):
        self.path                = path
        self.fileName            = fileName
        if (self.path is None):
            print("Just to double check you create instance of ClassCOMPAS without path/Data")
        elif not  os.path.isfile(path+fileName):
            raise ValueError("h5 file not found. Wrong path given? %s"\
                            %(path+fileName))
        elif os.path.isfile(path+fileName):
            pass



                              


        #Crucial values to be able to calculate MSSFR
        self.metallicityGrid     = None
        self.metallicitySystems  = None
        self.delayTimes          = None   #Myr
        #Crucial values I need for selection effects
        self.mass1               = None   #Msun
        self.mass2               = None   #Msun
        self.DCOmask             = None

        #Additional arrays that might be nice to store
        #to more quickly make some plots.
        #If you need more memory might help a tiny bit to not do
        self.lazyData            = lazyData
        self.mChirp              = None    #Msun
        self.q                   = None
        self.Hubble              = None

        #Needed to recover true solar mass evolved
        self.Mlower              = Mlower  #Msun
        self.Mupper              = Mupper  #Msun
        self.binaryFraction      = binaryFraction
        self.totalMassEvolvedPerZ = None   #Msun
        

        print("ClassCOMPAS: Remember to self.setGridAndMassEvolved()")
        print("                   then  self.setCOMPASDCOmask()")
        print("                   then  self.setCOMPASData()")
    
    def setCOMPASDCOmask(self, types='BBH', withinHubbleTime=True, optimistic=False):

            Data    = h5.File(self.path+self.fileName)
            fDCO    = Data['DoubleCompactObjects']
            if types == 'BBH':
                maskTypes = (fDCO['Stellar_Type_1'][()] == 14) &\
                            (fDCO['Stellar_Type_2'][()] == 14)
            elif types == 'BNS':
                maskTypes = (fDCO['Stellar_Type_1'][()] == 13) &\
                            (fDCO['Stellar_Type_2'][()] == 13)
            elif types == 'BHNS':
                maskTypes = ((fDCO['Stellar_Type_1'][()] == 14) &\
                            (fDCO['Stellar_Type_2'][()] == 13)) |\
                            ((fDCO['Stellar_Type_1'][()] == 13) &\
                            (fDCO['Stellar_Type_2'][()] == 14))
            else:
                raise ValueError('types=%s not in BBH, BNS, BHNS' %(types))
                    
            if withinHubbleTime == True:
                maskHubble = (fDCO['Merges_Hubble_Time'][()]==True)
            else:
                #Array where all are true
                maskHubble = np.ones(len(fDCO['Merges_Hubble_Time'][()]), dtype=bool)

            try:
                #we never want in first timestep after CEE, because 
                #we define it as a system that should not have survived the CEE
                RLOF_CEE = fDCO['Rlof_Secondary_Post_Common_Eenvelope'][()]
                maskNoRLOFafterCEE =  (RLOF_CEE==False)
            except KeyError:
                print("Warning no RLOF_SECONDARY_POST_COMMON_ENVELOPE column")
                print("I will not mask the data for this")
                #specifically use different column name as dummy 
                maskNoRLOFafterCEE =  np.ones(len(fDCO['Stellar_Type_1'][()]), dtype=bool)

            try:
                if optimistic == True:
                    #we do not care about the optimistic flag (both False and True allowed)
                    #Array where all are true, specifically use different column name
                    #in case the data does not have this column
                    maskOptimistic = np.ones(len(fDCO['Optimistic_Common_Envelope'][()]),\
                                     dtype=bool)
                else:
                    #optimistic scenario not allowed (pessimistic) hence the flag must be false
                    #This removes systems with CEE from HG donors (no core envelope separation)
                    maskOptimistic = fDCO['Optimistic_Common_Envelope'][()] == False
            except KeyError:
                print("Warning no Optimistic_Common_Envelope column")
                print("I will not mask the data for this")
                #specifically use different column name as dummy 
                maskOptimistic =  np.ones(len(fDCO['Stellar_Type_1'][()]), dtype=bool)

                              

            self.DCOmask = maskTypes & maskHubble & maskOptimistic #& maskNoRLOFafterCEE

            Data.close()

    def setGridAndMassEvolved(self):
    
        #The COMPAS simulation does not evolve all stars 
        #give me the correction factor for the total mass evolved
        #I assume each metallicity has the same limits, and does correction
        #factor, but the total mass evolved might be different.
        #This does not change when we change types and other masks this is 
        #general to the entire simulation so calculate once
        _, self.totalMassEvolvedPerZ =\
        MPZ.totalMassEvolvedPerZ(path=self.path, fileName=self.fileName, Mlower=self.Mlower, \
                                 Mupper=self.Mupper, binaryFraction=self.binaryFraction)
        #Want to recover entire metallicity grid, assume that every metallicity
        #evolved shows in all systems again should not change within same run
        #so dont redo if we reset the data
        Data = h5.File(self.path+self.fileName)
        metallicities =Data['SystemParameters']['Metallicity@ZAMS_1'][()]
        self.metallicityGrid     = np.unique(metallicities)
        Data.close()

    def setCOMPASData(self):
        Data    = h5.File(self.path+self.fileName)
        fDCO    = Data['DoubleCompactObjects']

        #Stuff I need for cosmological integral

        #sorry not the prettiest line is a boolean slice of seeds
        #this only works because seeds in systems file and DCO file are printed
        #in same order

        #Recover metallicity from initial parameters
        self.seedsDCO            = fDCO['SEED'][()][self.DCOmask]
        initialSeeds             = Data['SystemParameters']['SEED'][()]
        initialZ                 = Data['SystemParameters']['Metallicity@ZAMS_1'][()]
        maskMetallicity          = np.in1d(initialSeeds, self.seedsDCO)
        self.metallicitySystems  = initialZ[maskMetallicity]


        self.delayTimes          = np.add(fDCO['Time'][()][self.DCOmask] , \
                                          fDCO['Coalescence_Time'][()][self.DCOmask])
        self.mass1               = fDCO['Mass_1'][()][self.DCOmask]
        self.mass2               = fDCO['Mass_2'][()][self.DCOmask]

        #Stuff of data I dont need for integral
        #but I might be to laze to read in myself
        #and often use. Might turn it of for memory efficiency
        if self.lazyData:
            self.q                   = np.divide(self.mass2, self.mass1)
            boolq                    = self.mass2 > self.mass1
            self.q[boolq]            = np.divide(self.mass1[boolq], self.mass2[boolq])
            self.mChirp = np.divide((np.multiply(self.mass2, self.mass1)**(3./5.) ),\
                                           (np.add(self.mass2, self.mass1)**(1./5.)))
            self.Hubble              = fDCO['Merges_Hubble_Time'][...].squeeze()[self.DCOmask]

        Data.close()
    def recalculateTrueSolarMassEvolved(self, Mlower, Mupper, binaryFraction):
        #Possibility to test assumptions of True solar mass evolved
        self.Mlower              = Mlower
        self.Mupper              = Mupper
        self.binaryFraction      = binaryFraction
        _, self.totalMassEvolvedPerZ =\
        MPZ.totalMassEvolvedPerZ(pathCOMPASh5=self.path , Mlower=self.Mlower, \
                                 Mupper=self.Mupper, binaryFraction=self.binaryFraction)
