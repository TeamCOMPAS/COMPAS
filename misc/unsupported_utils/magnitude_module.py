"""
These are a set of functions to calculate the magnitudes of stars given luminosities, temperatures, spectral distributions, passband transmissivities, and wavelengths. These functions were used mainly to calculate the magnitudes of systems from COMPAS runs. 

If you would like to use a spectral distribution that is not the planck distribution you can go to https://pysynphot.readthedocs.io/en/latest/index.html. They provide various atmospheric models and explain how to access them. Once you have access to one of these atmospheric model spectral distributions you can use the functions in this file.

Note: PySynPhot does have functionality to create passbands and convolve them with any spectral distribution in PySynPhot, so you can try to just use PySynPhot instead of some functions in this file 

Note: Information for different passbands can be found at http://svo2.cab.inta-csic.es/theory/fps/

Note: These functions were originally used in a python notebook, so it might be best to use them in that way.
"""

import numpy as np

KB = 1.380649e-23 # Boltzmann's Constant in SI units
C = 2.99792458e8 # Speed of Light in SI units
H = 6.62607015e-34 # Planck's constant in SI units
SIGMA = 5.670374419e-8 # Stefan-Boltzmann Constant in SI units
SUNRADIUS = 6.957e8 # Radius of Sun in meters 
SUNLUM = 3.828e26 # in watts
SUNTEMP = 5778 # in Kelvin
STARDIST = 3.086e17 # 10 parsecs in meters

def getEnergyCoeff(spectrum, transms, wavelengths):
    '''
    Returns the average flux per wavelength of a star using a specified spectral distribution, for a given passband.
    
        Parameters:
            spectrum (1D array): The spectral distribution used for your source (units: watts/(m**3))
            transms (1D array): The transmissivities of the passband (units: None)
            wavelengths (1D array): The appropriate wavelengths for passband (units: meters)
            
        Returns:
            energyCoeff (number): The average flux per wavelength (units: watts/(m**2*nm))
    
    '''
    stepSize = getStepSize(wavelengths) # Units of meters
    
    denominatorIntegrand = np.multiply(transms, wavelengths)
    numeratorIntegrand = np.multiply(1e-9*spectrum, denominatorIntegrand)

    
    numerator = integrate(numeratorIntegrand, stepSize)
    denominator = integrate(denominatorIntegrand, stepSize)
    
    energyCoeff = np.divide(numerator, denominator)
    
    return energyCoeff


def getEnergyCoeffs(spectrums, transms, wavelengths):
    '''
    Returns an array of average flux per wavelength for an array of spectral distributions, for a given passband
    
        Parameters:
            spectrums (1D array): The spectral distribution used for your sources (units: watts/(m**3))
            transms (1D array): The transmissivities of the passband (units: None)
            wavelengths (1D array): The appropriate wavelengths for passband (units: meters)
            
        Returns:
            energyCoeffs (1D array): Array of average flux per wavelength (units: watts/(m**2*nm)
    '''
    energyCoeffs = np.array([getEnergyCoeff(spectrum, transms, wavelengths) for spectrum in spectrums])
    
    return energyCoeffs


def getAvgEnergy(lum, temp, energyCoeff, distToStar):
    '''
    Returns the average flux per wavelength that detector passband sees
    
        Parameters: 
            lum (number): Luminosity of star (units: solar luminosity)
            temp (number): Temperature of star (units: kelvin)
            energyCoeff (number): Average flux per wavelength that star emits at given temperature (units: watts/(m**2*nm))
            distToStar (number): The distance between the detector and the star (units: meters)
            
        Returns:
            avgEnergy (number): Average flux per wavelength that detector passband sees (units: watts/(m**2*nm))
            
    '''
    starRadius = getStarRadius(lum*SUNLUM, temp)
    avgEnergy = energyCoeff*(starRadius/distToStar)**2
    
    return avgEnergy


def getStarMagnitude(lum, temp, energyCoeff, distToStar, zeroPoint):
    '''
    Returns the magnitude of a star for a passband
    
        Parameters:
            lum (number): Luminosity of star (units: solar luminosity)
            temp (number): Temperature of star (units: kelvin)
            energyCoeff (number): Average flux per wavelength that star emits at given temperature (units: watts/(m**2*nm))
            distToStar (number): The distance between the detector and the star (units: meters)
            zeroPoint (number): zero point of passband (units: none)
            
        Returns:
            magnitude (number): Magnitude of a star for a passband (units: mag)
    '''
    avgEnergy = getAvgEnergy(lum, temp, energyCoeff, distToStar)
    #magnitude = -2.5*np.log10(avgEnergy)+zeroPoint
    magnitude = -2.5*np.log10(avgEnergy/zeroPoint)

    
    return magnitude


def getStarsMagnitudes(lums, temps, energyCoeffs, distToStar, zeroPoint):
    '''
    Returns magnitudes of stars for a passband
    
        Parameters:
            lums (1D array): Luminosities of stars (units: solar luminosity)
            temps (1D array): Temperatures of stars (units: kelvin)
            energyCoeffs (1D array): Average flux per wavelength that stars emit at given temperatures (units: watts/(m**2*nm))
            distToStar (number): The distance between the detector and the star (units: meters)
            zeroPoint (number): zero point of passband (units: none)
            
        Returns:
            magnitudes (1D array): Magnitudes of stars for a passband (units: mag)
    '''
    numValues = len(lums)
    
    magnitudes = np.array([getStarMagnitude(lums[i], temps[i], energyCoeffs[i], distToStar, zeroPoint) for i in range(numValues)])
    
    return magnitudes


def getStarMagnitudesOverTime(lums, temps, energyCoeffs, distToStar, zeroPoint):
    '''
    Returns magnitudes of stars for a passband over time (Different rows of 2D arrays correspond to different instances of time)
    
        Parameters:
            lums (2D array): Luminosities of stars over time (units: solar luminosity)
            temps (2D array): Temperatures of stars over time (units: kelvin)
            energyCoeffs (2D array): Average flux per wavelength that stars emit at given temperatures over time (units: watts/(m**2*nm))
            distToStar (number): The distance between the detector and the star (units: meters)
            zeroPoint (number): zero point of passband (units: none)
            
        Returns:
            magnitudes (2D array): Magnitudes of stars for a passband over time (units: mag)
    '''
    numTimes = len(lums)
    
    magnitudesOverTime = np.array([np.array(getStarsMagnitudes(lums[t], temps[t], energyCoeffs[t], distToStar, zeroPoint)) for t in range(numTimes)])
    
    return magnitudesOverTime


def getStarMagnitudeMap(lums, temps, energyCoeffs, distToStar, zeroPoint):
    '''
    Returns 2D map of magnitudes for an array of luminosities and temperatures, for a passband
    
        Parameters:
            lums (1D array): Range of luminosities
            temps (1D array): Range of temperatures
            energyCoeffs (1D array): Average flux per wavelength that stars emit at given temperatures (units: watts/(m**2*nm))
            distToStar (number): The distance between the detector and the star (units: meters)
            zeroPoint (number): zero point of passband (units: none)
            
        Returns:
            magGrid (2D array): A 2D array of magnitudes for each point on the luminosity-temperature coordinate grid (units: mag)
    '''
    numTemps = len(temps)
    numLums = len(lums)
    
    magGrid = np.zeros([numTemps, numLums])
    
    for j in range(numLums):
        for i in range(numTemps):
            magGrid[i][j] = getStarMagnitude(lums[i], temps[j], energyCoeffs[j], distToStar, zeroPoint)
                       
    return magGrid


def getBinaryMagnitude(lum1, lum2, temp1, temp2, energyCoeff1, energyCoeff2, distToStar, zeroPoint):
    '''
    Returns the magnitude of a binary for a passband
            
        Parameters:
            lum1 (number): Luminosity of star 1 (units: solar luminosity)
            lum2 (number): Luminosity of star 2 (units: solar luminosity)
            
            temp1 (number): Temperature of star 1 (units: kelvin)
            temp2 (number): Temperature of star 2 (units: kelvin)
            
            energyCoeff1 (number): Average flux per wavelength that star 1 emits at given temperature (units: watts/(m**2*nm))
            energyCoeff2 (number): Average flux per wavelength that star 2 emits at given temperature (units: watts/(m**2*nm))
            
            distToStar (number): The distance between the detector and the star (units: meters)
            zeroPoint (number): zero point of passband (units: none)
            
        Returns:
            binaryMagnitude (number): Magnitude of a binary for a passband (units: mag)
    
    '''
    starRadius1 = getStarRadius(lum1*SUNLUM, temp1)
    starRadius2 = getStarRadius(lum2*SUNLUM, temp2)
    
    avgEnergy1 = energyCoeff1*(starRadius1/distToStar)**2
    avgEnergy2 = energyCoeff2*(starRadius2/distToStar)**2
    
    avgEnergy = avgEnergy1 + avgEnergy2
    
    binaryMagnitude = -2.5*np.log10(avgEnergy/zeroPoint)
    
    return binaryMagnitude


def getPlanckDistribution(temp, wavelengths):
    '''
    Returns the Planck Distribution for a given temperature and set of wavelengths
        
        Parameters:
            temp (number): Temperature at which distribution is to be calculated (units: Kelvin)
            wavelengths (1D array): Wavelengths at which distribution is to be calculated (units: meters)
        
        Returns:
            planck_dist (1D array): 1D array of values of the planck distribution at the given temperature for the range of given wavelengths (units: watts/(m**3))
        
    '''
    exponential = np.exp(H*C/(wavelengths*KB*temp))
    planck_dist = (2*np.pi*H*C**2)/(wavelengths**5*(exponential-1))
    
    return planck_dist


def getStarRadius(lum, temp):
    '''
    Returns the radius of a star, given its temperature and luminosity
    
        Parameters:
            lum (number): Luminosity of star (units: watts)
            temp (number): Temperature of star (units: kelvin)
            
        Returns:
            starRadius (number): Radius of star with the given luminosity and temperature (units: meters)
    '''
    starRadius = np.sqrt(lum/(4*np.pi*SIGMA*temp**4))
    
    return starRadius


def integrate(array, stepSize):
    '''
    Returns the left-hand riemann sum of a function in an interval (approximate definite integral), given a step size
    
        Parameters:
            array (1D array): An array with the values of your function at a certain interval (units: anything)
            stepSize (number): The step size used to integrate your function (units: anything really)
            
        Returns:
            sum (number): The left-hand riemann sum of the function in the specified interval (units: anything)
    '''
    sum = 0
    
    for val in array:
        sum += val*stepSize
    
    return sum 


def getStepSize(array):
    '''
    Returns the step size of an equally spaced array 
    
        Parameters:
        array (1D array): Equally spaced array (units: anything)
        
        Returns:
        stepSize (number): step size of array (units: anything)
    '''
    stepSize = array[1] - array[0]
    return stepSize