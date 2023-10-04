import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy import units as u
from astropy import constants as c



class calculate_spin(object):
    """
    This class calculates the Black Hole (BH) or Neutron Star (NS) spin based on a given spin prescription
    It returns the spin of the compact object M1, and compact object M2    
    """

    
    def __init__(self, data_path=None, spin_model_name='uniform'):
    
        self.path                = data_path
        
        if (self.path is None):
            print("Warning: your hdf5 path is set to 'None'")
        elif not  os.path.isfile(data_path):
            raise ValueError("h5 file not found. Wrong path given?", "path given = %s"%data_path)
        elif os.path.isfile(data_path):
            self.h5file           = h5.File(data_path)
            
            
        self.spin_model_name = spin_model_name
        
        
    def random_uniform_spin(self, low, high):   
        """
        returns black hole spins with a random spin value uniformly sampled between 'low' and 'high'
        """
        
        import random
        
        sample_size = len(self.h5file['BSE_Double_Compact_Objects']['Mass(1)'][...].squeeze())
        
        self.spinM1 = np.random.uniform(low=low, high=high, size=sample_size)
        self.spinM2 = np.random.uniform(low=low, high=high, size=sample_size)
        
        
        return self.spinM1, self.spinM2
    
    
    def fixed_value_spin(self, spin_value):   
        """
        returns black hole spins with a fixed spin value 'spin_value'
        """
        M1samples = self.h5file['BSE_Double_Compact_Objects']['Mass(1)'][...].squeeze()
        
        self.spinM1 = np.ones_like(a=M1samples) * spin_value
        self.spinM2 = np.ones_like(a=M1samples) * spin_value
        
        return self.spinM1, self.spinM2
    
        
    def setCOMPASData(self):
        """ reads in some of the COMPAS parameters needed from hdf5 file """
        

        # """ reads in some of the COMPAS parameters needed from hdf5 file """

        fDCO      = self.h5file['BSE_Double_Compact_Objects'] # hdf5 file with the DCO information
        fSN       = self.h5file['BSE_Supernovae']  # hdf5 file with the SN information
        fini      = self.h5file['BSE_System_Parameters']  # hdf5 file with initial parameters
        # #
        self.SEED_dco = fDCO['SEED'][...].squeeze()


        
        mask_dco_SN  = np.in1d(fSN['SEED'][...].squeeze(), self.SEED_dco) # mask in the SNe files the SNe that correspond to our DCO               
        # first, remove simulteneous_SN
        whichSN = fSN['Supernova_State'][...].squeeze()[mask_dco_SN] 
        mask_no_simulteneous_SN = (whichSN[0::2]!=3) & (whichSN[1::2]!=3)
        
        # remove BNS
        self.st1 = fDCO['Stellar_Type(1)'][...].squeeze()
        self.st2 = fDCO['Stellar_Type(2)'][...].squeeze()
        
        mask_contains_BH = ((self.st1==14) | (self.st2==14)) & (mask_no_simulteneous_SN==1)
        
        
        mask_dco_ini  = np.in1d(fini['SEED'][...].squeeze(), self.SEED_dco[mask_contains_BH]) #  update mask removing simulteneous SNe
        mask_dco_SN   = np.in1d(fSN['SEED'][...].squeeze() , self.SEED_dco[mask_contains_BH]) #  update mask removing simulteneous SNe      

        self.st1 =  self.st1[mask_contains_BH]
        self.st2 =  self.st2[mask_contains_BH]
        self.M1 = fDCO['Mass(1)'][...].squeeze()[mask_contains_BH]   # Compact object mass [Msun] of the initially more massive star
        self.M2 = fDCO['Mass(2)'][...].squeeze()[mask_contains_BH]  # Compact object mass [Msun] of the initially less massive star

        
        
        mask_dco_ini  = np.in1d(fini['SEED'][...].squeeze(), self.SEED_dco[mask_contains_BH]) # mask in the initial parameter files the SNe that correspond to our DCO that are BHs
        mask_dco_SN  = np.in1d(fSN['SEED'][...].squeeze(), self.SEED_dco[mask_contains_BH]) # mask in the SNe files the SNe that correspond to our DCO that are BHs        
        
        
        self.metallicitySystems  = fini['Metallicity@ZAMS(1)'][...].squeeze()[mask_dco_ini] 

        whichSN = fSN['Supernova_State'][...].squeeze()[mask_dco_SN]   # this is 1 if the initially primary star goes SN and 2 if the secondary goes supernova
        whichSN2 = whichSN[1::2] # get whichStar for the first SN   (there are 2 SNe for all DCOs)       
        mask_SN_is_M1 = (whichSN==1)
        mask_SN_is_M2 = (whichSN==2) 
        
        self.separationPreSN2= fSN['SemiMajorAxis<SN'][...].squeeze()[mask_dco_SN][1::2] # the separation just before each SN  in [Rsun], we need only the separation for the second SN to occur, so the [1::2]  

        self.convert_a_to_P_circular(separation=self.separationPreSN2*u.Rsun, M1=self.M1*u.Msun, M2=self.M2*u.Msun)  # obtain the Period before the SNe
        self.PeriodPreSN2 = self.PeriodPreSN2.to(u.d).value
        self.MassCOM1CoreSN = fSN['Mass_CO_Core@CO(SN)'][...].squeeze()[mask_dco_SN][mask_SN_is_M1]   # obtain the CO core mass before the SNe when M2 goes SN
        self.MassCOM2CoreSN = fSN['Mass_CO_Core@CO(SN)'][...].squeeze()[mask_dco_SN][mask_SN_is_M2]   # obtain the CO core mass before the SNe when M1 goes SN

        self.spinM1, self.spinM2 = np.zeros_like(self.M1), np.zeros_like(self.M2)  # start by giving all primaries zero spin and all secondaries zero spin 
        # # did M1 form in the first SN?
        self.M1formedSecond =  (whichSN2==1) # mask that is 1 if the  compact object M1 formed first in the DCO
        # # did M2 form in the first SN?
        self.M2formedSecond =  (whichSN2==2)  # mask that is 1 if the compact object M2 formed first in the DCO
        # Wolf-Rayet mass 
        self.mWR =  fSN['Mass_Total@CO(SN)'][...].squeeze()[mask_dco_SN][1::2]   # obtain the WR mass before the SNe

    

        
    def function_f_Bavera21(self, c1, c2, c3):
        """
        calculates the f() function in Bavera et al (2021) based on given coefficients 
        m_WR with units using astropy
        """

        top = -c1
        bottom = c2 + np.exp(-c3*self.mWR)

        f = top/bottom


        return f  
    
   
    def convert_a_to_P_circular(self, separation, M1, M2):
        """calculate Period from separation
        separation is separation (needs to be given in astropy units)
        M1 and M2 are masses of the binary

        """
        G = c.G # [gr cm s^2]


        mu = G*(M1+M2)
        self.PeriodPreSN2 = 2*np.pi * np.sqrt(separation**3/mu)

          
    def calculate_alpha_beta_Bavera21(self):
        """ returns alpha and beta for the Bavera et al. (2021) spin calculation """
        # numerical coefficients form text below Eq 2
        # we use the values at helium depletion, since we later on use the C/O core mass. 
        c1_alpha, c2_alpha, c3_alpha =  0.059305, 0.035552, 0.270245
        c1_beta,  c2_beta, c3_beta   =  0.026960, 0.011001, 0.420739

        alpha = self.function_f_Bavera21(c1_alpha, c2_alpha, c3_alpha)
        beta  = self.function_f_Bavera21(c1_beta,  c2_beta,  c3_beta)

        return alpha, beta



        
        
    def Bavera21(self):
        """
        Returns spinM1 and spinM2, the spins of the compact objects formed from
        the initial most massive star (M1) and initial least massive star (M2), respectively. 
        
        In this approximation only a BH that is formed second can be tidally spun up, if its 
        pre-SN separation is tight enough. 

        based on Eq 1 and 2 from https://arxiv.org/pdf/2105.09077.pdf
    
    
        """
        
        # the following function reads in the COMPAS DCO parameters such as WR mass and sets self.spinM1, self.spinM2 
        # to arrays of zeros
        self.setCOMPASData() 
        
        ## Calculate spin for the ones that tidally spin up ##
        
        alpha, beta = self.calculate_alpha_beta_Bavera21()  #obtain alpha and beta 
        
        # define mask of BHs that can have spin through tidal spin up
        # if BH (self.st==14) & formed second (self.M1formedFirst==0) & has pre-SN period (WR-BH period) of < 1 day, calculate spin with the approximation
        maskGiveSpin1 = ((self.st1==14) & (self.M1formedSecond==1)) & (np.log10(self.PeriodPreSN2) < 1)
        maskGiveSpin2 = ((self.st2==14) & (self.M2formedSecond==1)) & (np.log10(self.PeriodPreSN2) < 1)
        
        # Bavera approximation: 
        first_term = (alpha* (np.log10(self.PeriodPreSN2)**2)) 
        second_term =  ( beta * np.log10(self.PeriodPreSN2))  
        self.spinM1[maskGiveSpin1]  =  first_term[maskGiveSpin1]  + second_term[maskGiveSpin1]  
        self.spinM2[maskGiveSpin2]  =  first_term[maskGiveSpin2]  + second_term[maskGiveSpin2] 
        
        #for some binaries the properties fall outside of the fit, we set these artificial to spin 0 
        n_negative_spin1 = np.sum(self.spinM1<0)
        self.spinM1[self.spinM1<0] = np.zeros(n_negative_spin1)
        n_negative_spin2 = np.sum(self.spinM2<0)
        self.spinM2[self.spinM2<0] = np.zeros(n_negative_spin2)

        
        return self.spinM1, self.spinM2
    
    
    
    
    def Qin21(self):
        """
        Returns spinM1 and spinM2, the spins of the compact objects formed from
        the initial most massive star (M1) and initial least massive star (M2), respectively. 
        
        In this approximation only a BH that is formed second can be tidally spun up, if its 
        pre-SN separation is tight enough. 
        
        see Qin+18, approximation originally given in https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.3682C 
        (and Equation 5 in https://arxiv.org/pdf/2103.02608.pdf)
        
        """
        
        # the following function reads in the COMPAS DCO parameters such as WR mass and sets self.spinM1, self.spinM2 
        # to arrays of zeros
        self.setCOMPASData() 
        
        m_, c_ = -5./3, 0.5 # from Qin + 2018 

        # if BH & formed second and tight enough at WR-BH phase, calculate spin with Qin+18 approximation for systems that have spin 1 
        maskGiveSpin1 = ((self.st1==14) & (self.M1formedSecond==1) &  (np.log10(self.PeriodPreSN2) < -0.3))
        maskGiveSpin2 = ((self.st2==14) & (self.M2formedSecond==1) &  (np.log10(self.PeriodPreSN2) < -0.3))
        # the systems above have maximum spin
        self.spinM1[maskGiveSpin1] = np.ones(np.sum(maskGiveSpin1)) # fill with ones 
        self.spinM2[maskGiveSpin2] = np.ones(np.sum(maskGiveSpin2)) # fill with ones 
  
        
        # now assign the spin for systems that lie in between the 0 and 1 spin using the fitting formulae
        maskGiveSpin1 = ((self.st1==14) & (self.M1formedSecond==1) &  (np.log10(self.PeriodPreSN2) >= -0.3) &  (np.log10(self.PeriodPreSN2) < 0.3))
        maskGiveSpin2 = ((self.st2==14) & (self.M2formedSecond==1) &  (np.log10(self.PeriodPreSN2) >= -0.3) &  (np.log10(self.PeriodPreSN2) < 0.3))

        self.spinM1[maskGiveSpin1] =  m_ * np.log10(self.PeriodPreSN2[maskGiveSpin1])  + c_   
        self.spinM2[maskGiveSpin2] =  m_ * np.log10(self.PeriodPreSN2[maskGiveSpin2])  + c_   
              
    
        return self.spinM1, self.spinM2

    
    
    def Z_dependent_spin_param(self, metallicity_range, prescription_name):
        """ 
        returns metallicity dependent spin parameter coefficients of fitting formula for the Geneva or MESA prescription 
        by Belczynski et al. (2020); https://www.aanda.org/10.1051/0004-6361/201936528. Prescription_name should be in ['Geneva', 'MESA']
        it either returns (Geneva:) b, m_1, m_2, a_low,  or (MESA:) a_1, b_1, m_1, a_2 and b_2
        """
        
        # return coefficients from paragraphs below equation 3 and 4 in Belczysnki et al. (2020) 
        if metallicity_range >= 0.010:
            return [2.258, 16.0, 24.2, 0.13] if prescription_name=='Geneva' else  [-0.0016, 0.115, np.inf, np.inf, np.inf]

        elif (metallicity_range >= 0.004) & (metallicity_range < 0.010):
            return [3.578, 31.0, 37.8, 0.25] if prescription_name=='Geneva' else  [-0.0006, 0.105, np.inf, np.inf, np.inf]

        elif (metallicity_range >= 0.0012) & (metallicity_range < 0.004):
            return [2.434, 18.0, 27.7, 0.0 ] if prescription_name=='Geneva' else  [0.0076, 0.050, -0.0019, 0.165, 12.09]


        elif (metallicity_range < 0.0012):
            return [3.666, 32.0, 38.8, 0.25] if prescription_name=='Geneva' else  [-0.0010, 0.125, np.inf, np.inf, np.inf]


    
    def Geneva(self):
        '''
        returns spins of the primary (spinM1) and secondary (spinM2) 
        based on the Geneva BH natal spin model as described in Belczsynski+20 
        See Equation 3 and Fig 1 in https://www.aanda.org/10.1051/0004-6361/201936528 
        this model assumes moderate BH birth spins, with a spin magnitude dependent on metallicity and CO core mass
        Warning: this function currently assumes the two stars have identical metallicity (add a Z mask to change this)
        ''' 

        # apply spin model on metallicity ranges that linearly span Z values by Belczysnki (cf. http://arxiv.org/abs/2301.01312) 
        Z_ranges = [[0.0, 0.0012], [0.0012, 0.004], [0.004, 0.010], [0.010, 1.0]]  # Tab. 1 conditions in http://arxiv.org/abs/2301.01312 

        self.setCOMPASData()  # read in required COMPAS data (CO mass, and metallicity)
        
        # loop over four metallicity ranges of prescription
        for _, Z_range in enumerate(Z_ranges):
            
            mask_inside_Z = (self.metallicitySystems >= Z_range[0]) & (self.metallicitySystems < Z_range[1])  # mask systems inside metallicity range
            b_, m1_, m2_, a_low = self.Z_dependent_spin_param(metallicity_range=Z_range[0], prescription_name='Geneva')  # obtains metallcity dependent model coefficients 
            
            ### apply Equation 3 in Belczynski et al. (2020) ###
            # set spins to max spin value of 0.85 for stars below lower CO core mass (m_1) condition:            
            self.spinM1[(self.MassCOM1CoreSN <= m1_)&(mask_inside_Z==1)], self.spinM2[(self.MassCOM2CoreSN <= m1_)&(mask_inside_Z==1)] = 0.85, 0.85
            
            # set spins to linear interpolation for CO masses in between m_1 and m_2 (a=-0.088 in eq. 3):
            mask_spin1 = (self.MassCOM1CoreSN > m1_) &  (self.MassCOM1CoreSN < m2_) & (mask_inside_Z==1)
            mask_spin2 = (self.MassCOM2CoreSN > m1_) &  (self.MassCOM2CoreSN < m2_) & (mask_inside_Z==1)
            self.spinM1[mask_spin1], self.spinM2[mask_spin2] = (-0.088*self.MassCOM1CoreSN[mask_spin1]) + b_, (-0.088*self.MassCOM2CoreSN[mask_spin2]) + b_
            
            # set spins to min spin value 'a_low' for stars above lower CO core mass (m_2) condition:
            self.spinM1[(self.MassCOM1CoreSN >= m2_)&(mask_inside_Z==1)], self.spinM2[(self.MassCOM2CoreSN >= m2_)&(mask_inside_Z==1)] = a_low, a_low


        return self.spinM1, self.spinM2
    

    
    def MESA(self):
        '''
        returns spins of the primary (spinM1) and secondary (spinM2) 
        based on the MESA BH natal spin model as described in Belczsynski et al. (2020) 
        See Equation 4 and Fig 2 in https://www.aanda.org/10.1051/0004-6361/201936528 
        this model assumes moderate BH birth spins, with a spin magnitude dependent on metallicity and CO core mass
        Warning: this function currently assumes the two stars have identical metallicity (add a Z mask to change this)
        ''' 
        # MESA spin prescription see from Table 1 in http://arxiv.org/abs/2301.01312 
        
        # apply spin model on metallicity ranges that linearly span Z values by Belczysnki (cf. http://arxiv.org/abs/2301.01312) 
        Z_ranges = [[0.0, 0.0012], [0.0012, 0.004], [0.004, 0.010], [0.010, 1.0]] # Tab. 1 conditions in  in http://arxiv.org/abs/2301.01312 

        self.setCOMPASData()  # read in required COMPAS data (CO mass, and metallicity)
        
        # loop over four metallicity ranges of prescription
        for _, Z_range in enumerate(Z_ranges):
            mask_inside_Z = (self.metallicitySystems >= Z_range[0]) & (self.metallicitySystems < Z_range[1])  # mask systems inside metallicity range
            a1_, b1_, a2_, b2_, m1_ = self.Z_dependent_spin_param(metallicity_range=Z_range[0], prescription_name='MESA')  # obtains metallcity dependent model coefficients 
            

            ### apply Equation 4 in Belczynski et al. (2020) ###
            # set spins to linear interpolation for CO masses in below m_1 condition:
            mask_spin1 = (self.MassCOM1CoreSN <= m1_)&(mask_inside_Z==1)
            mask_spin2 = (self.MassCOM2CoreSN <= m1_)&(mask_inside_Z==1)
            self.spinM1[mask_spin1], self.spinM2[mask_spin2] = (a1_*self.MassCOM1CoreSN[mask_spin1]) + b1_, (a1_*self.MassCOM2CoreSN[mask_spin2]) + b1_
             
            
            
            if np.isfinite(m1_):
                # set spins to second linear interpolation equation if we are in the Z range where m_1 is not infinite. 
                mask_spin1 = (self.MassCOM1CoreSN > m1_)&(mask_inside_Z==1)
                mask_spin2 = (self.MassCOM2CoreSN > m1_)&(mask_inside_Z==1)
                self.spinM1[mask_spin1], self.spinM2[mask_spin2] = (a2_*self.MassCOM1CoreSN[mask_spin1]) + b2_, (a2_*self.MassCOM2CoreSN[mask_spin2]) + b2_


        return self.spinM1, self.spinM2
    





class calculate_spin_olderCOMPASdata(object):
    """
    This class calculates the Black Hole (BH) or Neutron Star (NS) spin based on a given spin prescription
    It returns the spin of the compact object M1, and compact object M2    
    """

    
    def __init__(self, data_path=None, spin_model_name='uniform'):
    
        self.path                = data_path
        
        if (self.path is None):
            print("Warning: your hdf5 path is set to 'None'")
        elif not  os.path.isfile(data_path):
            raise ValueError("h5 file not found. Wrong path given?", "path given = %s"%data_path)
        elif os.path.isfile(data_path):
            self.h5file           = h5.File(data_path)
            
            
        self.spin_model_name = spin_model_name
        
        
    def random_uniform_spin(self, low, high):   
        """
        returns black hole spins with a random spin value uniformly sampled between 'low' and 'high'
        """
        
        import random
        
        sample_size = len(self.h5file['doubleCompactObjects']['M1'][...].squeeze())
        
        self.spinM1 = np.random.uniform(low=low, high=high, size=sample_size)
        self.spinM2 = np.random.uniform(low=low, high=high, size=sample_size)
        
        
        return self.spinM1, self.spinM2
    
    
    def fixed_value_spin(self, spin_value):   
        """
        returns black hole spins with a fixed spin value 'spin_value'
        """
        M1samples = self.h5file['doubleCompactObjects']['M1'][...].squeeze()
        
        self.spinM1 = np.ones_like(a=M1samples) * spin_value
        self.spinM2 = np.ones_like(a=M1samples) * spin_value
        
        
        return self.spinM1, self.spinM2
    
        
    def setCOMPASData(self):
        """ reads in some of the COMPAS parameters needed from hdf5 file """
        
        fDCO      = self.h5file['doubleCompactObjects'] # hdf5 file with the DCO information
        fSN       = self.h5file['supernovae']  # hdf5 file with the SN information
        #
        self.M1 = fDCO['M1'][...].squeeze()   # Compact object mass [Msun] of the initially more massive star
        self.M2 = fDCO['M2'][...].squeeze()  # Compact object mass [Msun] of the initially less massive star
        self.metallicitySystems  = fDCO['Metallicity1'][...].squeeze() 
        
        self.seedsDCO = fDCO['seed'][...].squeeze()  # get the seeds in the DCO file 
        self.seedsSN = fSN['randomSeed'][...].squeeze()    # get the seeds in the SN file 
        maskSNdco = np.in1d(self.seedsSN,  self.seedsDCO) # mask in the SNe files the SNe that correspond to our DCO
        whichSN = fSN['whichStar'][...].squeeze()[maskSNdco]   # this is 1 if the initially primary star goes SN and 2 if the secondary goes supernova
        whichSN2 = whichSN[1::2] # get whichStar for the first SN   (there are 2 SNe for all DCOs)       
        
        self.separationPreSN2= fSN['separationBefore'][...].squeeze()[maskSNdco][1::2] # the separation just before each SN  in [Rsun], we need only the separation for the second SN to occur, so the [1::2]  

        self.convert_a_to_P_circular(separation=self.separationPreSN2*u.Rsun, M1=self.M1*u.Msun, M2=self.M2*u.Msun)  # obtain the Period before the SNe
        self.PeriodPreSN2 = self.PeriodPreSN2.to(u.d).value
        self.MassCOM2CoreSN = fSN['MassCOCoreSN'][...].squeeze()[maskSNdco][1::2]   # obtain the CO core mass before the SNe
        self.MassCOM1CoreSN = fSN['MassCOCoreSN'][...].squeeze()[maskSNdco][0::2]
        
        self.st1 = fDCO['stellarType1'][...].squeeze()   # obtain the final stellar type of the Primary 
        self.st2 = fDCO['stellarType2'][...].squeeze()   # obtain the final stellar type of the Secondary
        
        self.spinM1, self.spinM2 = np.zeros_like(self.M1), np.zeros_like(self.M2)  # start by giving all primaries zero spin and all secondaries zero spin 
        # did M1 form in the first SN?
        self.M1formedSecond =  (whichSN2==1) # mask that is 1 if the  compact object M1 formed first in the DCO
        # did M2 form in the first SN?
        self.M2formedSecond =  (whichSN2==2)  # mask that is 1 if the compact object M2 formed first in the DCO
        mask_SN1not1or2 = (whichSN2!=1) & (whichSN2!=2)
        
        
        self.mWR =  fSN['MassStarSN'][...].squeeze()[maskSNdco][1::2]   # obtain the CO core mass before the SNe

        
    def function_f_Bavera21(self, c1, c2, c3):
        """
        calculates the f() function in Bavera et al (2021) based on given coefficients 
        m_WR with units using astropy
        """

        top = -c1
        bottom = c2 + np.exp(-c3*self.mWR)

        f = top/bottom


        return f  
    
   
    def convert_a_to_P_circular(self, separation, M1, M2):
        """calculate Period from separation
        separation is separation (needs to be given in astropy units)
        M1 and M2 are masses of the binary

        """
        G = c.G # [gr cm s^2]


        mu = G*(M1+M2)
        self.PeriodPreSN2 = 2*np.pi * np.sqrt(separation**3/mu)

          
    def calculate_alpha_beta_Bavera21(self):
        """ returns alpha and beta for the Bavera et al. (2021) spin calculation """
        # numerical coefficients form text below Eq 2
        # we use the values at helium depletion, since we later on use the C/O core mass. 
        c1_alpha, c2_alpha, c3_alpha =  0.059305, 0.035552, 0.270245
        c1_beta,  c2_beta, c3_beta   =  0.026960, 0.011001, 0.420739

        alpha = self.function_f_Bavera21(c1_alpha, c2_alpha, c3_alpha)
        beta  = self.function_f_Bavera21(c1_beta,  c2_beta,  c3_beta)

        return alpha, beta



        
        
    def Bavera21(self):
        """
        Returns spinM1 and spinM2, the spins of the compact objects formed from
        the initial most massive star (M1) and initial least massive star (M2), respectively. 
        
        In this approximation only a BH that is formed second can be tidally spun up, if its 
        pre-SN separation is tight enough. 

        based on Eq 1 and 2 from https://arxiv.org/pdf/2105.09077.pdf
    
    
        """
        
        # the following function reads in the COMPAS DCO parameters such as WR mass and sets self.spinM1, self.spinM2 
        # to arrays of zeros
        self.setCOMPASData() 
        
        ## Calculate spin for the ones that tidally spin up ##
        
        alpha, beta = self.calculate_alpha_beta_Bavera21()  #obtain alpha and beta 
        
        # define mask of BHs that can have spin through tidal spin up
        # if BH (self.st==14) & formed second (self.M1formedFirst==0) & has pre-SN period (WR-BH period) of < 1 day, calculate spin with the approximation
        maskGiveSpin1 = ((self.st1==14) & (self.M1formedSecond==1)) & (np.log10(self.PeriodPreSN2) < 1)
        maskGiveSpin2 = ((self.st2==14) & (self.M2formedSecond==1)) & (np.log10(self.PeriodPreSN2) < 1)
        
        # Bavera approximation: 
        first_term = (alpha* (np.log10(self.PeriodPreSN2)**2)) 
        second_term =  ( beta * np.log10(self.PeriodPreSN2))  
        self.spinM1[maskGiveSpin1]  =  first_term[maskGiveSpin1]  + second_term[maskGiveSpin1]  
        self.spinM2[maskGiveSpin2]  =  first_term[maskGiveSpin2]  + second_term[maskGiveSpin2] 
        
        #for some binaries the properties fall outside of the fit, we set these artificial to spin 0 
        n_negative_spin1 = np.sum(self.spinM1<0)
        self.spinM1[self.spinM1<0] = np.zeros(n_negative_spin1)
        n_negative_spin2 = np.sum(self.spinM2<0)
        self.spinM2[self.spinM2<0] = np.zeros(n_negative_spin2)
        
        
        return self.spinM1, self.spinM2
    
    
    
    
    def Qin21(self):
        """
        Returns spinM1 and spinM2, the spins of the compact objects formed from
        the initial most massive star (M1) and initial least massive star (M2), respectively. 
        
        In this approximation only a BH that is formed second can be tidally spun up, if its 
        pre-SN separation is tight enough. 
        
        see Qin+18, approximation originally given in https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.3682C 
        (and Equation 5 in https://arxiv.org/pdf/2103.02608.pdf)
        
        """
        
        # the following function reads in the COMPAS DCO parameters such as WR mass and sets self.spinM1, self.spinM2 
        # to arrays of zeros
        self.setCOMPASData() 
        
        m_, c_ = -5./3, 0.5 # from Qin + 2018 

        # if BH & formed second and tight enough at WR-BH phase, calculate spin with Qin+18 approximation for systems that have spin 1 
        maskGiveSpin1 = ((self.st1==14) & (self.M1formedSecond==1) &  (np.log10(self.PeriodPreSN2) < -0.3))
        maskGiveSpin2 = ((self.st2==14) & (self.M2formedSecond==1) &  (np.log10(self.PeriodPreSN2) < -0.3))
        # the systems above have maximum spin
        self.spinM1[maskGiveSpin1] = np.ones(np.sum(maskGiveSpin1)) # fill with ones 
        self.spinM2[maskGiveSpin2] = np.ones(np.sum(maskGiveSpin2)) # fill with ones 
  
        
        # now assign the spin for systems that lie in between the 0 and 1 spin using the fitting formulae
        maskGiveSpin1 = ((self.st1==14) & (self.M1formedSecond==1) &  (np.log10(self.PeriodPreSN2) >= -0.3) &  (np.log10(self.PeriodPreSN2) < 0.3))
        maskGiveSpin2 = ((self.st2==14) & (self.M2formedSecond==1) &  (np.log10(self.PeriodPreSN2) >= -0.3) &  (np.log10(self.PeriodPreSN2) < 0.3))

        self.spinM1[maskGiveSpin1] =  m_ * np.log10(self.PeriodPreSN2[maskGiveSpin1])  + c_   
        self.spinM2[maskGiveSpin2] =  m_ * np.log10(self.PeriodPreSN2[maskGiveSpin2])  + c_   
              
    
        return self.spinM1, self.spinM2

    
    
    def Z_dependent_spin_param(self, metallicity_range, prescription_name):
        """ 
        returns metallicity dependent spin parameter coefficients of fitting formula for the Geneva or MESA prescription 
        by Belczynski et al. (2020); https://www.aanda.org/10.1051/0004-6361/201936528. Prescription_name should be in ['Geneva', 'MESA']
        it either returns (Geneva:) b, m_1, m_2, a_low,  or (MESA:) a_1, b_1, m_1, a_2 and b_2
        """
        
        # return coefficients from paragraphs below equation 3 and 4 in Belczysnki et al. (2020) 
        if metallicity_range >= 0.010:
            return [2.258, 16.0, 24.2, 0.13] if prescription_name=='Geneva' else  [-0.0016, 0.115, np.inf, np.inf, np.inf]

        elif (metallicity_range >= 0.004) & (metallicity_range < 0.010):
            return [3.578, 31.0, 37.8, 0.25] if prescription_name=='Geneva' else  [-0.0006, 0.105, np.inf, np.inf, np.inf]

        elif (metallicity_range >= 0.0012) & (metallicity_range < 0.004):
            return [2.434, 18.0, 27.7, 0.0 ] if prescription_name=='Geneva' else  [0.0076, 0.050, -0.0019, 0.165, 12.09]


        elif (metallicity_range < 0.0012):
            return [3.666, 32.0, 38.8, 0.25] if prescription_name=='Geneva' else  [-0.0010, 0.125, np.inf, np.inf, np.inf]


    
    def Geneva(self):
        '''
        returns spins of the primary (spinM1) and secondary (spinM2) 
        based on the Geneva BH natal spin model as described in Belczsynski+20 
        See Equation 3 and Fig 1 in https://www.aanda.org/10.1051/0004-6361/201936528 
        this model assumes moderate BH birth spins, with a spin magnitude dependent on metallicity and CO core mass
        Warning: this function currently assumes the two stars have identical metallicity (add a Z mask to change this)
        ''' 

        # apply spin model on metallicity ranges that linearly span Z values by Belczysnki (cf. http://arxiv.org/abs/2301.01312) 
        Z_ranges = [[0.0, 0.0012], [0.0012, 0.004], [0.004, 0.010], [0.010, 1.0]]  # Tab. 1 conditions in http://arxiv.org/abs/2301.01312 

        self.setCOMPASData()  # read in required COMPAS data (CO mass, and metallicity)
        
        # loop over four metallicity ranges of prescription
        for _, Z_range in enumerate(Z_ranges):
            
            mask_inside_Z = (self.metallicitySystems >= Z_range[0]) & (self.metallicitySystems < Z_range[1])  # mask systems inside metallicity range
            b_, m1_, m2_, a_low = self.Z_dependent_spin_param(metallicity_range=Z_range[0], prescription_name='Geneva')  # obtains metallcity dependent model coefficients 
            
            ### apply Equation 3 in Belczynski et al. (2020) ###
            # set spins to max spin value of 0.85 for stars below lower CO core mass (m_1) condition:
            self.spinM1[(self.MassCOM1CoreSN <= m1_)&(mask_inside_Z==1)], self.spinM2[(self.MassCOM2CoreSN <= m1_)&(mask_inside_Z==1)] = 0.85, 0.85
            
            # set spins to linear interpolation for CO masses in between m_1 and m_2 (a=-0.088 in eq. 3):
            mask_spin1 = (self.MassCOM1CoreSN > m1_) &  (self.MassCOM1CoreSN < m2_) & (mask_inside_Z==1)
            mask_spin2 = (self.MassCOM2CoreSN > m1_) &  (self.MassCOM2CoreSN < m2_) & (mask_inside_Z==1)
            self.spinM1[mask_spin1], self.spinM2[mask_spin2] = (-0.088*self.MassCOM1CoreSN[mask_spin1]) + b_, (-0.088*self.MassCOM2CoreSN[mask_spin2]) + b_
            
            # set spins to min spin value 'a_low' for stars above lower CO core mass (m_2) condition:
            self.spinM1[(self.MassCOM1CoreSN >= m2_)&(mask_inside_Z==1)], self.spinM2[(self.MassCOM2CoreSN >= m2_)&(mask_inside_Z==1)] = a_low, a_low


        return self.spinM1, self.spinM2
    

    
    def MESA(self):
        '''
        returns spins of the primary (spinM1) and secondary (spinM2) 
        based on the MESA BH natal spin model as described in Belczsynski et al. (2020) 
        See Equation 4 and Fig 2 in https://www.aanda.org/10.1051/0004-6361/201936528 
        this model assumes moderate BH birth spins, with a spin magnitude dependent on metallicity and CO core mass
        Warning: this function currently assumes the two stars have identical metallicity (add a Z mask to change this)
        ''' 
        # MESA spin prescription see from Table 1 in http://arxiv.org/abs/2301.01312 
        
        # apply spin model on metallicity ranges that linearly span Z values by Belczysnki (cf. http://arxiv.org/abs/2301.01312) 
        Z_ranges = [[0.0, 0.0012], [0.0012, 0.004], [0.004, 0.010], [0.010, 1.0]] # Tab. 1 conditions in  in http://arxiv.org/abs/2301.01312 

        self.setCOMPASData()  # read in required COMPAS data (CO mass, and metallicity)
        
        # loop over four metallicity ranges of prescription
        for _, Z_range in enumerate(Z_ranges):
            mask_inside_Z = (self.metallicitySystems >= Z_range[0]) & (self.metallicitySystems < Z_range[1])  # mask systems inside metallicity range
            a1_, b1_, a2_, b2_, m1_ = self.Z_dependent_spin_param(metallicity_range=Z_range[0], prescription_name='MESA')  # obtains metallcity dependent model coefficients 
            

            ### apply Equation 4 in Belczynski et al. (2020) ###
            # set spins to linear interpolation for CO masses in below m_1 condition:
            mask_spin1 = (self.MassCOM1CoreSN <= m1_)&(mask_inside_Z==1)
            mask_spin2 = (self.MassCOM2CoreSN <= m1_)&(mask_inside_Z==1)
            self.spinM1[mask_spin1], self.spinM2[mask_spin2] = (a1_*self.MassCOM1CoreSN[mask_spin1]) + b1_, (a1_*self.MassCOM2CoreSN[mask_spin2]) + b1_
             
            
            
            if np.isfinite(m1_):
                # set spins to second linear interpolation equation if we are in the Z range where m_1 is not infinite. 
                mask_spin1 = (self.MassCOM1CoreSN > m1_)&(mask_inside_Z==1)
                mask_spin2 = (self.MassCOM2CoreSN > m1_)&(mask_inside_Z==1)
                self.spinM1[mask_spin1], self.spinM2[mask_spin2] = (a2_*self.MassCOM1CoreSN[mask_spin1]) + b2_, (a2_*self.MassCOM2CoreSN[mask_spin2]) + b2_


        return self.spinM1, self.spinM2
    
        
    

        
    
