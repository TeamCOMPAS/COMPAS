#!/usr/bin/env python3
import numpy as np
import scipy.special
from   scipy.optimize    import newton
from   astropy.cosmology import FlatLambdaCDM
from   astropy.cosmology import WMAP9 #as cosmo
import astropy.units as u




class MSSFR(object):
    """
    This class is to calculate the metallicity specific star formation
    rate in a specific metallicity bin.

    It combines a
    -  MZ -relation     : galaxy stellar mass - metallicity relatiion
                          This translates a metallicity to a 
                          a galaxy stellar mass which has the same 
                          average metallicity.
    -  SFR-prescription : Star formation rate prescription.
                          Amount of solar mass that goes into forming
                          stars.
    -  GSMF function    : Galaxy stellar mass function
                          A density function of the distribution of galaxy stellar masses.


    The entire GSMF is the mass in all stars in all galaxies. Hence the fraction of the GSMF
    is the fraction of the stellar mass in all galaxies at that time. 
    We assume that the average cosmic SFR spreads evenly among all the galaxy stellar mass.
    Meaning a galaxy with a stellar mass twice that of another galaxy, has twice
    the SFR. Hence this method does not include local star burts (SMC/LMC for example)
    

    The overall idea. Given a metallicity bin with upper and lower bound.
    Translate this to a galaxy stellar mass bin with upper and lower bound using
    a MZ-relation. Then find the fraction of the GSMF density function that is taken by 
    this glaxy stellar mass bin. We assumed that the fraction of the stellar mass is 
    proportional to the fraction of the SFR. Hence we multiply this fraction by the SFR.
    Note that because we take "cosmological averages" this method does not take into 
    account galaxy specific metallicity distributions.


    Alternatively we can also directly use a metallicity density/probability function
    and calculate the fraction occupied by the metallicity bin in this distribution. This
    is the logNormal distribution, which is still multiplied by the SFR and a 
    'cosmological average'.

    """

    def __init__(self, verbose=False,\
                 metallicityGrid=None,  binInLogSpace=True,\
                 metallicityLowerLimit=1E-90, metallicityUpperLimit=1.,\
                 logOH12sun=8.69,     solarMetallicity=0.0142, cosmo=None):

        #if True we print statements
        self.verbose               = verbose

        #With regards to metallicity binning the integral
        self.metallicityGrid       = metallicityGrid
        self.metallicityBinEdges   = None  
        self.binInLogSpace         = binInLogSpace 
        #lower limit is not zero to avoid divide by zero error
        self.metallicityLowerLimit = metallicityLowerLimit
        self.metallicityUpperLimit = metallicityUpperLimit

        print("ClassMSSFR: Remember to set self.SFRprescription() + subparts")
        print("                            self.Zprescription()   +sub-parts")

        if self.metallicityGrid is not None:
            self.calculateMetallicityBinEdges()
        else:
            print("                            set self.calculateMetallicityBinEdges()")


        #import constants default from is Asplund
        self.logOH12sun            = logOH12sun        #fraction number density solar
        self.solarMetallicity      = solarMetallicity  #fraction mass in metals solar

        #if you use ZM_GSMF or LogNormal metallicity distribution
        self.Zprescription         = None
        #current Options  = logNormal' or 'MZ_GSMF

        #If you use MZ and GSMF relations
        self.GSMFprescription      = None
        # Current Options = Panter et al. (2004) Single, 
        #                   Furlong et al. (2015) Single,
        #                   Furlong et al. (2015) Double
        self.ZMprescription        = None
        # Current Options = Ma et al. (2015),
        #                   Langer et al. (2006)
        #                   Langer et al. +offset (2006)


        #if you use a logNormal    
        self.logNormalPrescription = None
        self.customLogNormal       = [None, None, None] #[Z0, alpha, sigma]

        #Current Options  = Phenomenological

        # SFR prescriptions
        self.SFRprescription       = None 
        self.customSFR             = [None, None, None, None] #[a,b,c,d] see function customSFR
        #Current Options  = Madau et al. (2014)
        #                   Madau et al. (2017)
        #                   Strolger et al (2004)


        self.cosmology             = None
        #by allowing to pass the cosmology class
        #we only have to set it once in the project reducing mistakes
        if cosmo == 'oldDefault':
            self.cosmology = FlatLambdaCDM(H0=67.8 *\
                                        u.km / u.s / u.Mpc, Om0=0.308)
        elif cosmo is not None:
            self.cosmology = cosmo
        else:
            #Have not defined a cosmology jet use a default
            self.cosmology = WMAP9 #as cosmo

        
        
    def calculateMetallicityBinEdges(self):
        """Calculates the bins used in the Riemann sum over metallicities

        """

        if self.binInLogSpace:
            logMetallicities = np.log10(self.metallicityGrid)
            b=  logMetallicities[:-1] + (logMetallicities[1:] - logMetallicities[:-1])/2.
            b = 10.**b #the boundaries for integration are not in log space so
                   #convert to "normal" numbers.
        else:
            b=  (self.metallicityGrid[1:] - self.metallicityGrid[:-1])/2. \
                + self.metallicityGrid[:-1] 

        self.metallicityBinEdges = np.zeros(len(b)+2)

        #the lowest/highest metallicity bin edge are set in options
        #the calculated b edges are all in between

        self.metallicityBinEdges[0]    = self.metallicityLowerLimit
        self.metallicityBinEdges[-1]   = self.metallicityUpperLimit
        self.metallicityBinEdges[1:-1] = b



    







    def LogOH12vsLogZZsun(self, value, inValue='logOH12'):
        """
        This function translates between 

        logOH12  = number density of oxygen to hyrdogen
        logZZsun = metallicity mass fraction in solar units

        asssumptions: logOH12sun
        The solar number density from oxygen to hydrogen.
        This assumption is non-trivial and can affect results
        """
        if (inValue == 'logZZsun'):
            #translate from logZZsun to logOH12
            logOH12 = value + self.logOH12sun
            return logOH12
        elif (inValue == 'logOH12'):
            #translate from logOH12 to logZZsun
            logZZsun = (value-self.logOH12sun)
            return logZZsun


    ########################################################
    #                                                      #
    #                                                      #
    #       Everything related to GSMF ZM relations        #
    #                                                      #
    #                                                      #
    ########################################################
    def returnFractionMZ_GSMF(self, Zlower, Zupper, redshift):

        """
        Returns the fraction of the SFR.
        It combines a
        -  MZ -relation     : galaxy stellar mass - metallicity relatiion
                              This translates a metallicity to a 
                              a galaxy stellar mass which has the same 
                              average metallicity.
        -  GSMF function    : Galaxy stellar mass function
                              A density function of the distribution of galaxy stellar masses.

        The entire GSMF is the mass in all stars in all galaxies.
        Hence the fraction of the GSMF is the fraction of the stellar
        mass in all galaxies at that time. We assume that the 
        average cosmic SFR spreads evenly among all the galaxy stellar mass.
        Meaning a galaxy with a stellar mass twice that of another galaxy, has twice
        the SFR. Hence this method does not include local star burts (SMC/LMC for example)
        
        The overall idea. Given a metallicity bin with upper and lower bound.
        Translate this to a galaxy stellar mass bin with upper and lower bound using
        a MZ-relation. Then find the fraction of the GSMF density function that is taken by 
        this glaxy stellar mass bin. We assumed that the fraction of the stellar mass is 
        proportional to the fraction of the SFR. Hence we multiply this fraction by the SFR.
        Note that because we take "cosmological averages" this method does not take into 
        account galaxy specific metallicity distributions.

        """
        z = np.copy(redshift) #found that the mask z>4 overwrites input

        #these prescriptions work in metallicities of unit solar
        Zupper_Zsun = Zupper/ self.solarMetallicity 
        Zlower_Zsun = Zlower/ self.solarMetallicity 

        #Not compressed coding but easier to read


        if self.ZMprescription == 'Ma et al. (2016)':
            Mupper = 10**self.Ma2015ZM(Zupper_Zsun, z)
            Mlower = 10**self.Ma2015ZM(Zlower_Zsun, z)
        elif self.ZMprescription == 'Langer et al. (2006)':
            Mupper = 10**self.Langer2005ZM(Zupper_Zsun, z)
            Mlower = 10**self.Langer2005ZM(Zlower_Zsun, z)
        elif self.ZMprescription == 'Langer et al. +offset (2006)':
            Mupper = 10**self.Langer2005OffsetZM(Zupper_Zsun, z)
            Mlower = 10**self.Langer2005OffsetZM(Zlower_Zsun, z)          
        elif self.ZMprescription == 'Savaglio2005Bisector':
            print("Danger Savaglio2005Bisector not fully working not tested")
            Mupper = 10**self.Sava2005BisectorZM(Zupper_Zsun, z)
            Mlower = 10**self.Sava2005BisectorZM(Zlower_Zsun, z)
        elif self.ZMprescription == 'Savaglio2005Bisector':
            print("Danger Savaglio2005Bisector not fully working turnover in function")
            Mupper = 10**self.Savaglio2005ZM(Zupper_Zsun, z)
            Mlower = 10**self.Savaglio2005ZM(Zlower_Zsun, z)
        else:     
             raise ValueError("The mass-metallicity prescription is not recognised.\n"+\
                              "Current (trusted) options are:\n"+\
                              "Ma et al. (2016), Langer et al. (2006),"+\
                              "Langer et al. +offset (2006)\n"+\
                              "Apologies but I will now break")
        
        #To integrate the schechter function we need to know
        #- alpha which governs the power law slope
        #- phi for the normalization (in case of the double schechter)
        #- Mc which determines the location of the turnover
        if   self.GSMFprescription == 'Panter et al. (2004) Single': #REDSHIFT INDEPENDENT!
            logMc, phi, alpha = self.PanterSingleRedshiftIndependent(z)
        elif self.GSMFprescription == 'Furlong et al. (2015) Single':
            logMc, phi, alpha = self.lineairFitSingleSchechterFurlong(z)
        elif self.GSMFprescription == 'Furlong et al. (2015) Double':
            logMc,  phi1, a1, phi2, a2 = self.lineairFitDoubleSchechterFurlong(z)
        else:
            raise ValueError( "This GSMF prescription is not recognised.\n"+\
                              "Current options are:\n"+\
                              "Panter et al. (2004) Single, Furlong et al. (2015) Single, Furlong et al. (2015) Double \n"+\
                              "Apologies but I will now break")

        Mc = 10**logMc

        #The schechter can either be a single or double schechter function
        if self.GSMFprescription in ['Panter et al. (2004) Single','Furlong et al. (2015) Single']:
            schechterUpperIncomplete = scipy.special.gammainc(alpha+2,  Mupper/Mc)
            schechterLowerIncomplete = scipy.special.gammainc(alpha+2,  Mlower/Mc)
            fraction                 = schechterUpperIncomplete - schechterLowerIncomplete

        elif self.GSMFprescription in ['Furlong et al. (2015) Double']:
            #again could do trick to make smaller but this is more
            #robust when passing either integers/floats/arrays as redshift
            schechter1Complete   = scipy.special.gamma(a1+2.0)
            schechter1Incomplete = scipy.special.gammainc(a1+2.0, Mupper/(Mc)) * schechter1Complete
            schechter2Complete   = scipy.special.gamma(a2+2.0)
            schechter2Incomplete = scipy.special.gammainc(a2+2.0, Mupper/(Mc)) * schechter2Complete
            integralUpperBound   = (phi1*schechter1Incomplete + phi2*schechter2Incomplete) /\
                                   (phi1*schechter1Complete   + phi2*schechter2Complete)
            #lower bound
            schechter1Complete   = scipy.special.gamma(a1+2.0)
            schechter1Incomplete = scipy.special.gammainc(a1+2.0, Mlower/(Mc)) * schechter1Complete
            schechter2Complete   = scipy.special.gamma(a2+2.0)
            schechter2Incomplete = scipy.special.gammainc(a2+2.0, Mlower/(Mc)) * schechter2Complete
            integralLowerBound   = (phi1*schechter1Incomplete + phi2*schechter2Incomplete) /\
                                   (phi1*schechter1Complete   + phi2*schechter2Complete)
            fraction             = np.subtract(integralUpperBound, integralLowerBound)
        return fraction

    #####################################################
    #                                                   #
    #           GSMF Schechter functions                #
    #                                                   #
    #####################################################  
    def lineairFitSingleSchechterFurlong(self, z):

        try:
            nr = len(z)
        except:
            z  = np.array([z])
            nr = 1
            
        fitlogMc        = np.array([11.14, 11.11, 11.06, 10.91, 10.78, 10.60])
        fitphi1          = np.array([0.84, 0.84, 0.74, 0.45, 0.22, 0.12])*10**(-3)
        fita1            = np.array([-1.43, -1.45, -1.48, -1.57, -1.66, -1.74])

        
        fitredshifts = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 4.0])
        
        fitvalues = [fitlogMc, fitphi1, fita1]
        thresholds = [0.0, 0.5, 1.0, 2.0, 3.0, 10000000]
        r          = []
        for nrv, values in enumerate(fitvalues):
            dz     = np.diff(fitredshifts)
            dydz   = np.divide(np.diff(values), dz)
            yvalues= np.zeros(nr)
            for nrz, redshift in enumerate(thresholds[:-1]):
                mask = (z>= thresholds[nrz]) & (z<=thresholds[nrz+1])
                if nrz == 0 :
                    #interpolate from z0.5 down
                    dz   = 0.5 - z[mask]
                    interpolatedValue  = values[nrz+1] - np.multiply(dz, dydz[nrz])
                else:
                    #interpolate up
                    dz   = z[mask] - redshift
                    interpolatedValue  = values[nrz] + np.multiply(dz, dydz[nrz])
                yvalues[mask] =  interpolatedValue
            r.append(yvalues)
            
        logMc, phi1, a1 = r
        #Ignoring phi fit since we do not use the implicit SFR values
        phi1 = np.ones(len(a1))
        #This fudge makes the schechter function solvable with gamma functions
        #This happens at redshifts above 7. Else a <-2 makes Gamma(a+2)=Nan
        #In reality this means that we alter the fraction
        a1[a1<-1.99]=-1.99
        return logMc, phi1, a1


    def lineairFitDoubleSchechterFurlong(self, z):

        try:
            nr = len(z)
        except:
            z  = np.array([z])
            nr = 1

        
        fitredshifts = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 4.0])
        fitlogMc     = np.array([10.95, 10.88, 10.74, 10.44, 10.19, 10.00])
        fitphi1      = np.array([1.45, 1.61, 1.51, 1.06, 0.63, 0.24])*10**-3
        fita1        = np.array([-1.31, -1.24, -0.98, -0.25, 0.23, 0.43])
        fitphi2      = np.array([0.0,  0.08,   0.48, 0.8, 0.61, 0.43])*10**-3
        fita2        = np.array([-2.22, -1.79, -1.62, -1.58, -1.64, -1.69])
        
        ratioPhi1Phi2 = np.divide(fitphi2, fitphi1) #Do it in this order because phi2 is zero
        
        fitvalues  = [fitlogMc, fitphi1, fita1, fitphi2, fita2, ratioPhi1Phi2]
        thresholds = [0.0, 0.5, 1.0, 2.0, 3.0, 10000000]
        r          = []
        for nrv, values in enumerate(fitvalues):
            dz     = np.diff(fitredshifts)
            dydz   = np.divide(np.diff(values), dz)
            yvalues= np.zeros(nr)
            for nrz, redshift in enumerate(thresholds[:-1]):
                mask = (z>= thresholds[nrz]) & (z<=thresholds[nrz+1])
                if nrz == 0 :
                    #interpolate from z0.5 down
                    dz   = 0.5 - z[mask]
                    interpolatedValue  = values[nrz+1] - np.multiply(dz, dydz[nrz])
                else:
                    #interpolate up
                    dz   = z[mask] - redshift
                    interpolatedValue  = values[nrz] + np.multiply(dz, dydz[nrz])
                yvalues[mask] =  interpolatedValue
            r.append(yvalues)


        logMc, phi1, a1, phi2, a2, ratio = r
        ratio[ratio<0.]=0.
        phi1 = np.ones(len(ratio))
        phi2 = ratio
        #When phi == 0 the integral doesn matter, but is still done
        #Hence when the interpolation gives an a<-2. even if we do not use
        #it it will brake the code  Gammafunction(a+2) = inf/nan for a<-2.
        #So I pass dummy value
        a2[phi2 == 0]  =0.

        #dont want metallicity to increase again at z<0.5
        mask = (z<0.5) & (a2<-1.79)
        a2[mask]=-1.79
        return logMc, phi1, a1, phi2, a2




    def PanterSingleRedshiftIndependent(self,z):
        #GSMF from panter et al 2004
        #https://academic.oup.com/mnras/article/355/3/764/952806
        #Used in Norman & Langer  2006
        #http://iopscience.iop.org/article/10.1086/500363/meta
        #Note that this is redshift independent!
        logMc = np.log10(7.64*(10**10))
        phi   = 7.8* (10**(-3)) 
        alpha = -1.16
        return np.ones(len(z))*logMc, np.ones(len(z))*phi, np.ones(len(z))*alpha

    #####################################################
    #                                                   #
    #           ZM - relations (inverse of papers)      #
    #           MZ-relations outside class bottom file  #
    #                                                   #
    #####################################################  
    def Ma2015ZM(self, ZZSun, z):
        logM = ((np.log10(ZZSun)  -7.95 + self.logOH12sun -0.93*\
               (np.exp(-0.43*z)) )/0.35  +10.)
        return logM

    def Langer2005ZM(self, ZZsun, z):
        Mstar  =  7.64*10**10
        logM   =  np.log10((((ZZsun*(10**(0.3*z)))**2)*Mstar)) 
        return logM

    def Langer2005OffsetZM(self, ZZsun, z):
        ZZsun = ZZsun * 10**(0.3*(z-0.7))
        M= (ZZsun**2.09) * (10**(-8.49+2.09*self.logOH12sun))
        return np.log10(M)

    def Sava2005BisectorZM(self, ZZsun, z):
        #used for z=0.7
        #Gemini Deep survey
        logZZsun =  np.log10(ZZsun)
        logOH12  =  self.LogOH12vsLogZZsun(logZZsun, inValue='logZZsun')
        logM     =  (logOH12-4.062)/(0.478)
        return logM

    def Savaglio2005ZM(self, ZZsun, z):
        logTh   = np.log10(self.cosmology.age(z).value)
        logZZsun= np.log10(ZZsun)
        logOH12 =  self.LogOH12vsLogZZsun(logZZsun, inValue='logZZsun')
        a       = -0.09649
        b       = -0.4030*logTh + 2.5315
        c       = -7.5903 -logOH12  + 5.1733*logTh - 0.3944*logTh*logTh 
        d       = b*b - 4*a*c
        mask1   = d<0.
        mask2   = np.logical_not(mask1)
        logM    = np.zeros(len(ZZsun))
        logM[mask2] = ((-b) + np.sqrt(d[mask2])) / (2*a) #The ones with a non-negative root
                                                            #are the solutions of interest
        logM[mask1]= None  #Set the ones with complex solution to -1
                           #Need to account for these later on
        return logM




    #################################################
    #                                               #
    #           SFR -  Prescriptions                #
    #                                               #
    #################################################
    def returnSFR(self, redshifts, age):
        if (self.SFRprescription   == 'Madau et al. (2014)'):
            SFR = self.SFR_Madau(redshifts)
        elif (self.SFRprescription   == 'Madau et al. (2017)'):
            SFR = self.SFR_Madau2(redshifts)
        elif (self.SFRprescription == 'Strolger et al. (2004)'):
            SFR = self.SFR_Strolger(age)
        elif (self.SFRprescription == 'Neijssel et al. (2019)'):
            SFR = self.SFR_Neijssel(redshifts)
        elif (self.SFRprescription == 'Custom SFR'):
            SFR = self.SFR_Custom(redshifts)
        else:
            raise ValueError( "This SFR prescription is not recognised.\n"+\
                              "Current options are:\n"+\
                              "Madau et al. (2014),Madau et al. (2017),Neijssel et al. (2019), Strolger et al. (2004) \n"+\
                              "Apologies but I will now break")
        return SFR

    def SFR_Madau(self, z):
        """
        https://arxiv.org/pdf/1403.0007v3.pdf pg 48 eq(15)

        """
        SFR = 0.015* ((1+z)**2.7) / ( 1 + ((1+z)/2.9)**5.6) * 1e9 #[1e9 for GPc-3]
        return SFR # [Msun yr-1 Gpc-3] in comoving volume

    def SFR_Madau2(self, z):
        """
        https://arxiv.org/pdf/1606.07887.pdf
        """
        SFR = 0.01 * ((1+z)**2.6) / (1 + ((1+z)/3.2)**6.2) * 1e9 #[1e9 for GPc-3]
        return SFR # [Msun yr-1 Gpc-3] in comoving volume

    def SFR_Strolger(self, tGyrs):
        """
        https://arxiv.org/pdf/1308.1546.pdf pg 2 eq(1)
        """
        #extinction corrected model
        a,b,c,d, t0=0.182, 1.26, 1.865, 0.071, 13.47


        SFR = 10**9 *a * (tGyrs**b * np.exp(-tGyrs/c) + d*np.exp(d*(tGyrs-t0)/c))
        if np.isnan(SFR.any()):
            raise ValueError("Nan in SFR calculation for %s" %(self.SFRprescription))
        return SFR #Msun yr-1 Gpc-3 in comoving volume

    def SFR_Neijssel(self, z):
        """
        fingerspitzengefuhl Cosmic integration paper
        """
        SFR = 0.01 * ((1+z)**2.77) / (1 + ((1+z)/2.9)**4.7) * 1e9 #[1e9 for GPc-3]
        return SFR # [Msun yr-1 Gpc-3] in comoving volume        

    def SFR_Custom(self, z):
        """
        Custom SFR same functional form as Madau et al
        """
        a = self.customSFR[0]            
        b = self.customSFR[1] 
        c = self.customSFR[2] 
        d = self.customSFR[3] 
        SFR = a * ((1+z)**b) / (1 + ((1+z)/c)**d) * 1e9 #[1e9 for GPc-3]
        return SFR # [Msun yr-1 Gpc-3] in comoving volume
    #################################################
    #                                               #
    #           MZ - relations                      #
    #                                               #
    #################################################   

    def Tremmonti2004MZ(self, logM, z):
        #z=0.1  validity 8.5 < logM  < 11.5
        #SDSS Sloan digital sky survey
        #R23 for metallicity, with Monte-Carlo draw for starFormation History
        #Compared against survey 53.000 star-forming galaxies at z=0.1
        logOH12  = -1.492 + 1.847*logM - 0.08026*logM*logM
        logZZsun = self.LogOH12vsLogZZsun(logOH12)
        return logZZsun


    def Savaglio2005MZ(self, logM, z):
        #used galaxies from 0.4 < z < 1.0 plots go  from 7 < logM < 11
        #Gemini Deep survey
        logTh = np.log10(self.cosmology.age(z).value)
        logOH12 = -7.5903 + 2.5315*logM + \
                  -0.09649*logM*logM  + 5.1733*logTh - 0.3944*logTh*logTh + \
                  -0.4030*logTh*logM
        logZZsun = self.LogOH12vsLogZZsun(logOH12)
        return logZZsun

    def Ma2015MZ(self, logM, z):
        #simulations at z=0 validity expected z0-6  4 < logMstellar < 11
        #implicitly assumes logOH12sun = 9.0 I made it variable so we can use our own logOH sun
        logZZsun = 0.35*(logM - 10) + 0.93*np.exp(-0.43*z) +7.95 - self.logOH12sun
        return logZZsun

    def Langer2005MZ(self, logM, z):
        M_star  = (7.64*10**10)
        Z_Zsun = np.sqrt((10**logM)/(M_star))* (10**(-0.3*z))
        return np.log10(Z_Zsun)

    def Sava2005BisectorMZ(self, logM, z):
        #Furthermore the ZM relations of Savaglio goes from z 0.4-1.9
        #only valid at z=0.7
        #Gemini Deep survey
        logOH12 = 0.478*logM+4.062
        logZZsun = self.LogOH12vsLogZZsun(logOH12)
        return logZZsun


    def returnFractionLogNormal(self, Zlower, Zupper, redshift):
        """
        Instead of combining a GSMF and MZ relation
        you can also directly construct a redshift dependent metallicity distribution

        In this case we focus on a log-normal distribution

        log10(Z) is normally distributed with:

        meanM  = log(Zmean)
                 This is redshift dependent and Zmean scales as
                 Zmean = Z0 * 10^(alpha * redshift)
                 Z0 is the mean metallicity at redshift zero
                 alpha is the scaling of the redshift

        sigma  = the standard deviation of log(Z) distribution
                 we assume for now that this is redshift independent

        """
        if self.logNormalPrescription == 'Neijssel Phenomenological':
            "Based on norman Langer"
            Z0       = 0.035
            alpha    = -0.23
            sigma    = 0.39
        elif self.logNormalPrescription == 'Custom Phenomenological':

            if None in self.customLogNormal:
                raise ValueError("customLogNormal not set .customLogNormal =[Z0, alpha, sigma]")
            Z0       = self.customLogNormal[0]
            alpha    = self.customLogNormal[1]
            sigma    = self.customLogNormal[2]

        else:
            raise ValueError( "This logNormal prescription is not recognised.\n"+\
                              "Current options are:\n"+\
                              "Phenomenological, Custom Phenomenological' \n"+\
                              "Apologies but I will now break")
        Zmean    = Z0 * (10**(alpha*redshift))
        mu       = np.log(Zmean)  - (sigma * sigma) /2.
        Xupper   = (np.log(Zupper) - mu)/float(np.sqrt(2)*sigma)
        CDFUpper = 0.5 + 0.5 * scipy.special.erf(Xupper)
        Xlower   = (np.log(Zlower) - mu)/float(np.sqrt(2)*sigma)
        CDFLower =  0.5 + 0.5 * scipy.special.erf(Xlower)
        fraction = CDFUpper - CDFLower 
        return fraction


    def returnMSSFR(self, metallicity=None, agesBirth=None, redshiftBirth=None):
        """
        This is the main function called in the cosmic integration routine.

        For a given binNr in the metallicity bins it calculates the MSSFR

        """
        # find the bin number
        lowerBinNr      = np.where(self.metallicityGrid == metallicity)[0][0]

        #we will only calculate the systems of which ages
        # and redshifts are not -1 (flagged as born too soon)
        mask = agesBirth != -1
        #independent of metallicity so we can do all at once
        SFR        = np.zeros(len(agesBirth))
        SFR[mask]  = self.returnSFR(redshiftBirth[mask], agesBirth[mask])
        #The calculation of the fraction is done such that it does 
        #a single metallicity bin at the time
        fractions = np.zeros(len(agesBirth))

        Zlower = self.metallicityBinEdges[lowerBinNr]
        Zupper = self.metallicityBinEdges[lowerBinNr+1]
        if self.Zprescription == 'MZ_GSMF':
            fractions[mask] = self.returnFractionMZ_GSMF(Zlower, Zupper, redshiftBirth[mask])
        elif self.Zprescription == 'logNormal':
            fractions[mask] = self.returnFractionLogNormal(Zlower, Zupper, redshiftBirth[mask])
        else:
            raise ValueError("self.Zprescription not recognized, MZ_GSMF or logNormal?")
        #Always raise an error if the fraction is unphysical
        if ((fractions>1).any()) or ((fractions<0).any()):
            raise ValueError("unphysical SFR fraction check metallicity binning limits")
        #i.e. if system is born before first star formation (z=-1)
        #we set the rate to zero
        MSSFR  = np.multiply(SFR, fractions)
        return MSSFR



    def printSFRoptions(self):

        print("""
        Default instance:
        self.SFRprescription       = None 

        pass string
        Current Options  = 'Madau et al. (2014)'
                           'Madau et al. (2017)'
                           'Strolger et al. (2004)'
                           'Neijssel et al. (2019)'
                           'Custom SFR'

        If you use Custom SFR also set the constants
        you want to use

        self.customSFR = [a,b,c,d] 
        
        which are used in a Madau et al like formula (see source code)

        """)

    def printZMoptions(self):


        print("""
        If you use MZ and GSMF relations
        i.e. self.Zprescription    = 'MZ_GSMF'

        Default instance
        self.ZMprescription        = None
        Current Options = Ma et al. (2016),
                      Langer et al. (2006)
                      Langer et al. +offset (2006)
        """)



    def printGSMFoptions(self):


        print("""
        If you use MZ and GSMF relations
        i.e. self.Zprescription    = 'MZ_GSMF'

        Default instance
        self.GSMFprescription      = None
        Current Options = Panter et al. (2004) Single, 
                          Furlong et al. (2015) Single,
                          Furlong et al. (2015) Double
        """)


    def printLogNormaloptions(self):


        print("""
        If you use Log normal 
        i.e. self.Zprescription    = 'logNormal'

        Default instance
        self.logNormalPrescription      = None
        Current Options = 'Neijssel Phenomenological'
                          'Custom Phenomenological'

        If you use Custom Phenomenological remember to 
        additionally set 

        self.customLogNormal = [Z0, alpha, sigma]
        """)

