import numpy as np 
import h5py as h5 #for reading in data
import os

def threePartBrokenPowerLaw(x, x1=0.01, x2=0.08, x3=0.5, x4=200, a1=-0.3, \
                            a2=-1.3, a3=-2.3, C1=1):
    #Not that everything outside the range x1<x4 is set to zero
    yvalues = np.zeros(len(x))
    
    #calculate values of the x values that are x1<=x<x2
    mask1 = (x>=x1) & (x<x2)
    yvalues[mask1] = C1 * (x[mask1]**a1)
    
    #calculate values of the x values that are x2<=x<x3
    mask2 = (x>=x2) & (x<x3)
    C2    = C1 * (x2**(a1-a2))
    yvalues[mask2] = C2 * (x[mask2]**a2)
    
    #calculate values of the x values that are x3<=x<=x4
    mask3 = (x>=x3) & (x<=x4)
    C3    = C1 * (x2**(a1-a2)) * (x3**(a2-a3))
    yvalues[mask3] = C3 * (x[mask3]**a3)
    
    return yvalues

def CDFbrokenPowerLaw(x, x1, x2, x3, x4, a1, a2, a3, C1):
    yvalues = np.zeros(len(x))
    
    C2    = float(C1 * (x2**(a1-a2)))
    C3    = float(C2 * (x3**(a2-a3)))
    
    N1 = float(((1./(a1+1)) * C1 * (x2**(a1+1))) - ((1./(a1+1)) * C1 * (x1**(a1+1))))
    N2 = float(((1./(a2+1)) * C2 * (x3**(a2+1))) - ((1./(a2+1)) * C2 * (x2**(a2+1))))
    N3 = float(((1./(a3+1)) * C3 * (x4**(a3+1))) - ((1./(a3+1)) * C3 * (x3**(a3+1))))
    
    bottom = N1+N2+N3
    
    mask1 = (x>=x1) & (x<x2)
    top1 = (((1./(a1+1)) * C1 * (x[mask1]**(a1+1)) - (1./(a1+1)) * C1 * (x1**(a1+1))))
    yvalues[mask1] = top1/bottom
    
    #calculate values of the x values that are x2<=x<x3
    mask2 = (x>=x2) & (x<x3)
    top2 =  N1 + (((1./(a2+1)) * C2 * (x[mask2]**(a2+1)) - (1./(a2+1)) * C2 * (x2**(a2+1))))
    yvalues[mask2] = top2/bottom
    
    #calculate values of the x values that are x3<=x<=x4
    mask3 = (x>=x3) & (x<=x4)
    top3 =  N1 + N2 + (((1./(a3+1)) * C3 * (x[mask3]**(a3+1)) - (1./(a3+1)) * C3 * (x3**(a3+1))))
    yvalues[mask3] = top3/bottom
    return yvalues


def invertCDFbrokenPowerLaw(CDF, x1, x2, x3, x4, a1, a2, a3, C1):
    #I specifically do floats against python rounding when dividing
    
    #The constants needed
    C2    = float(C1 * (x2**(a1-a2)))
    C3    = float(C2 * (x3**(a2-a3)))
    
    N1 = float(((1./(a1+1)) * C1 * (x2**(a1+1))) - ((1./(a1+1)) * C1 * (x1**(a1+1))))
    N2 = float(((1./(a2+1)) * C2 * (x3**(a2+1))) - ((1./(a2+1)) * C2 * (x2**(a2+1))))
    N3 = float(((1./(a3+1)) * C3 * (x4**(a3+1))) - ((1./(a3+1)) * C3 * (x3**(a3+1))))
    
    bottom = N1+N2+N3
    
    CDFx2 = CDFbrokenPowerLaw(np.array([x2,x2]), x1, x2, x3, x4, a1, a2, a3, C1)[0]
    CDFx3 = CDFbrokenPowerLaw(np.array([x3,x3]), x1, x2, x3, x4, a1, a2, a3, C1)[0]

    
    xvalues = np.zeros(len(CDF))
    
    mask1 = (CDF < CDFx2)
    xvalues[mask1] =  (((CDF[mask1]*(N1+N2+N3))  + \
                      ( (1./(a1+1))*C1*(x1**(a1+1))))/((1./(a1+1))*C1))**(1./(a1+1))
    
    mask2 = (CDFx2<= CDF) & (CDF < CDFx3)
    xvalues[mask2] = ((((CDF[mask2]*(N1+N2+N3))-(N1))  + \
                      ( (1./(a2+1))*C2*(x2**(a2+1))))/((1./(a2+1))*C2))**(1./(a2+1))
    
    mask3 = (CDFx3<= CDF) 
    xvalues[mask3] = ((((CDF[mask3]*(N1+N2+N3))-(N1+N2))  + \
                      ((1./(a3+1))*C3*(x3**(a3+1))))/((1./(a3+1))*C3))**(1./(a3+1))
    
    return xvalues





def createSampleUniverse(binaryFraction=1., x1=0.01, x2=0.08, x3=0.5, x4=200, a1=-0.3, \
                            a2=-1.3, a3=-2.3, C1=1, sampleSize=5000000, Mmin=0.01, Mmax=200):
    
    binaryFraction = binaryFraction

    #Given the defined three-part broken powerlaw,
    #We can sample a subset by using Mmin, Mmax.
    #We convert it to the respective values between
    #0-1 and only sample uniformly between those (woohoo :D)
    #Mmin and Mmax have to be between x1 and x4
    CDFmin = CDFbrokenPowerLaw(np.array([Mmin]), x1, x2, x3, x4, a1, a2, a3, C1)
    CDFmax = CDFbrokenPowerLaw(np.array([Mmax]), x1, x2, x3, x4, a1, a2, a3, C1)
    
    #All the random drawing that we need
    drawM1         = np.random.uniform(CDFmin,CDFmax,sampleSize)
    drawBinary     = np.random.uniform(0,1,sampleSize)
    drawM2         = np.random.uniform(0,1,sampleSize)

    #All the arrays we want to fill
    M1 = np.zeros(sampleSize)
    M2 = np.zeros(sampleSize)

    #Define the IMF broken powerlaw and calculate masses from inverted CDF


    M1 = invertCDFbrokenPowerLaw(drawM1, x1, x2, x3, x4, a1, a2, a3, C1)

    #Binary fraction is easier, since we draw between 0-1, every draw with 
    #value above binary fraction = single star and every value below is binary
    #for a single star we set M2=0 Msun.
    #Note that we assume that the binary Fraction is mass independent
    #Future work to implenet Max Moe ps and qs options
    maskBinary = drawBinary < binaryFraction  #booleans

    #again for the secondary we assume the mass ratio distribution to be flat across
    #the whole parameter range so then the drawM2 (if it is in a binary) 
    #just becomes the mass fraction.

    M2[maskBinary] = np.multiply(drawM2[maskBinary],M1[maskBinary])
    #all the ones outside the mask remain zero
    return M1, M2


def inverseCDF(C, CDF, index, xmin, xmax):
    #CDF sincle powerlaw
    a =  (1./(index + 1)) * C * CDF**(index+1)
    b =  (1./(index + 1)) * C * xmin**(index+1)
    c =  (1./(index + 1)) * C * xmax**(index+1)
    top    = ((CDF * (c-b))+b)*(index + 1)
    bottom = C
    return (top/bottom)**(1./(index+1))



def retrieveMassEvolvedPerZ(path):
    f = h5.File(path, 'r') # open in read-only

    allSystems = f['BSE_System_Parameters']
    metals = (allSystems['Metallicity@ZAMS(1)'])[()]
    m1s = (allSystems['Mass@ZAMS(1)'])[()]
    m2s = (allSystems['Mass@ZAMS(2)'])[()]
    total = []
    for Z in np.unique(metals):
        mask = metals == Z
        total.append(np.sum(m1s[mask]) + np.sum(m2s[mask]))
    f.close()
    return np.array(total)



def totalMassEvolvedPerZ(path=None, Mlower=None, Mupper=None, binaryFraction=0.7, \
                         x1=0.01, x2=0.08, x3=0.5, x4=200., a1=-0.3, a2=-1.3, a3=-2.3, C1=1., Mmax=200):

    #the default values assume a Kroupa IMF for M1
    if path is None:
        raise TypeError("\n Need to give path to data")
    if Mlower is None:
        raise TypeError("\n Need to give lower limit M1 of pythonSubmit")
    if Mupper is None:
        raise TypeError("\n Need to give upper limit M1 of pythonSubmit")
    

    M1, M2 = createSampleUniverse(binaryFraction=binaryFraction, x1=x1, x2=x2, x3=x3, x4=x4, \
                                  a1=a1, a2=a2, a3=a3, C1=C1, Mmax=Mmax)

    totalMassInStarFormation = np.sum(M1) + np.sum(M2)

    #Now mask M1 and M2 to see what lies in the range of COMPAS
    maskM1 = (M1>=Mlower) & (M1<=Mupper)
    maskBinaries = (M2!=0)
    mask = maskM1 & maskBinaries

    totalMassEvolvedCOMPAS = np.sum(M1[mask]) + np.sum(M2[mask])

    #multiplication fraction
    fraction = totalMassEvolvedCOMPAS/float(totalMassInStarFormation)

    #so we need to muliplu the mass evolved per metallicity times (1/fraction)
    #to know the total mass evolved per metallicity
    MassEvolvedPerZ = retrieveMassEvolvedPerZ(path)

    multiplicationFactor = 1./fraction

    totalMassEvolvedPerMetallicity = (MassEvolvedPerZ)/(fraction)
    return multiplicationFactor, totalMassEvolvedPerMetallicity



