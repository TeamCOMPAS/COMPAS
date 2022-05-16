"""
Sampler for the correlated parameter distributions of 
primary mass, mass ratio, period, and eccentricity 
(M1, q, P, e) of isolated singles, binaries, and higher 
order systems on the Main Sequence, following Moe & DiStefano 2017.

This script is modified from the original IDL script 
by Max Moe, converted into python by Mads Sorensen, and
adapted for COMPAS by Reinhold Willcox.
See: https://arxiv.org/abs/1808.06488

The output of this script is a grid file suitable for 
COMPAS input of parameters for single stars, isolated
binaries, and the inner binary of higher order systems.
The properties of the outer components of higher order 
systems are not calculated, but they can be understood
to be stable since they come directly from the observations,
which must necessarily only include stable systems 
(and not, e.g, chaotic triples).

Due to the correlated nature of these parameters, it is not 
straightforward to supply alternative constraints, e.g
including only short period systems. Thus the only free
parameters are the maximum and minimum total mass, 
and the number of systems sampled. To change these,
set them within the `if __name__ == "__main__":`
function at the bottom of this file.

See the original comments below for clarity.
"""


# -*- coding: utf-8 -*-

import numpy as np
import scipy.integrate

def idl_tabulate(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in range(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret

#-----------------------------------
# Belows is the adapted version of Maxwell Moe's IDL code that generates a population
# of single and binary stars based on the paper Mind your P's and Q's
# By Maxwell Moe and Rosanne Di Stefano
#
# The python code has been adopted by Mads Sørensen
#-----------------------------------
# Version history:
# V. 0.1 2017/02/03
# By Mads Sørensen
# - This is a pure adaption from IDL to Python.
# - The function idl_tabulate is similar to
# the IDL function int_tabulated except, this function seems to be slightly
# more exact in its solution.
# Therefore, relative to the IDL code, there are small numerical differences.
# - The output is written to the console.
#-----------------------------------
#
# Comments below beginning with  is the original nodes by Maxwell Moe.
# Please read these careful for understanding the script.
# NOTE - This version produces only the statistical distributions of
#       single stars, binaries, and inner binaries in hierarchical triples.
#       Outer tertiaries in hierarchical triples are NOT generated.
#       Moreover, given a set of companions, all with period P to
#       primary mass M1, this version currently uses an approximation to
#       determine the fraction of those companions that are inner binaries
#       vs. outer triples. Nevertheless, this approximation reproduces
#       the overall multiplicity statistics.
#
#
# Step 1 - Tabulate probably density functions of periods,
#         mass ratios, and eccentricities based on
#         analytic fits to corrected binary star populations.
# Step 2 - Implement Monte Carlo method to generate stellar
#         population from those density functions.
#


def createParameterDistributionsAndSampler():

    """
    #
    # Step 1
    #
    # Vector of primary masses M1 (Msun), logarithmic orbital period P (days),
    # mass ratios q = Mcomp/M1, and eccentricities e
    #
    """

    # 0.8 < M1 < 40 (where we have statistics corrected for selection effects)
    M1v = 10**(np.linspace(0,1,101)*(np.log10(40.)-np.log10(0.8))+np.log10(0.8))
    
    # 0.15 < log P < 8.0
    logPv = np.linspace(0,157,158)*0.05+0.15
    
    # 0.10 < q < 1.00
    qv = np.linspace(0.1,1,91)
    
    # 0.0001 < e < 0.9901
    # set minimum to non-zero value to avoid numerical errors
    ev = np.linspace(0,0.99,100)+0.0001
    
    # Note that companions outside this parameter space (e.g., q < 0.1,
    # log P (days) > 8.0) are not constrained in M+D16 and therefore
    # not considered.
    
    
    # Number of elements within each respective vector
    numM1=101
    numlogP=158
    numq=91
    nume=100
    
    
    # Distribution functions - define here, but evaluate within for loops.
    
    # Frequency of companions with q > 0.1 per decade of orbital period.
    # Bottom panel in Fig. 36 of M+D16
    flogP_sq=np.zeros([numlogP,numM1])
    
    
    # Given M1 and P, the cumulative distribution of mass ratios q
    cumqdist = np.zeros([numq,numlogP,numM1])
    
    # Given M1 and P, the cumulative distribution of eccentricities e
    cumedist = np.zeros([nume,numlogP,numM1])
    
    # Given M1 and P, the probability that the companion
    # is a member of the inner binary (currently an approximation).
    # 100% for log P < 1.5, decreases with increasing P
    probbin = np.zeros([numlogP,numM1])
    
    
    # Given M1, the cumulative period distribution of the inner binary
    # Normalized so that max(cumPbindist) = total binary frac. (NOT unity)
    cumPbindist = np.zeros([numlogP,numM1])
    
    
    # Slope alpha of period distribution across intermediate periods
    # 2.7 - DlogP < log P < 2.7 + DlogP, see Section 9.3 and Eqn. 23.
    # Slightly updated from version 1.
    alpha=0.018
    DlogP=0.7
    
    
    # Heaviside function for twins with 0.95 < q < 1.00
    H=qv*0.
    ind=np.where(qv >= 0.95)
    H[ind]=1.0
    H=H/idl_tabulate(qv,H) #normalize so that integral is unity
    
    
    # Relevant indices with respect to mass ratio
    indlq=np.where(qv >= 0.3)
    indsq=np.where(qv < 0.3)
    indq0p3=np.min(indlq)
    
    
    # Loop through primary mass
    for i in range(0, numM1):
        myM1=M1v[i]
        # Twin fraction parameters that are dependent on M1 only section 9.1
        FtwinlogPle1=0.3-0.15*np.log10(myM1)   # Eqn. 6
        logPtwin=8.-myM1                       # Eqn. 7a
        if (myM1 >= 6.5):
            logPtwin=1.5                       # Eqn. 7b
    
        # Frequency of companions with q > 0.3 at different orbital periods
        # and dependent on M1 only section 9.3 (slightly modified since v1)
        flogPle1   = 0.020+0.04*np.log10(myM1)+0.07*(np.log10(myM1))**2.   # Eqn. 20
        flogPeq2p7 = 0.039+0.07*np.log10(myM1)+0.01*(np.log10(myM1))**2.   # Eqn. 21
        flogPeq5p5 = 0.078-0.05*np.log10(myM1)+0.04*(np.log10(myM1))**2.   # Eqn. 22
    
        # Loop through orbital period P
        for j in range(0, numlogP):
            mylogP = logPv[j]
    
    
            # Given M1 and P, set excess twin fraction section 9.1 and Eqn. 5
            if(mylogP <= 1.):
                Ftwin=FtwinlogPle1
            else:
                Ftwin=FtwinlogPle1*(1.- (mylogP-1)/(logPtwin-1.))
            if(mylogP >= logPtwin):
                Ftwin=0.
    
    
            # Power-law slope gamma_largeq for M1 < 1.2 Msun and various P Eqn. 9
            if(mylogP <= 5.0):
                gl_1p2 = -0.5
            if(mylogP > 5.0):
                gl_1p2 = -0.5-0.3*(mylogP-5.0)
    
            # Power-law slope gamma_largeq for M1 = 3.5 Msun and various P Eqn. 10
            if(mylogP <= 1.0):
                gl_3p5 = -0.5
            if((mylogP > 1.0)and(mylogP <= 4.5)):
                gl_3p5 = -0.5-0.2*(mylogP-1.0)
            if((mylogP > 4.5)and(mylogP <= 6.5)):
                gl_3p5 = -1.2-0.4*(mylogP-4.5)
            if(mylogP > 6.5):
                gl_3p5 = -2.0
    
            # Power-law slope gamma_largeq for M1 > 6 Msun and various P Eqn. 11
            if(mylogP <= 1.0):
                gl_6 = -0.5
            if((mylogP > 1.0)and(mylogP <= 2.0)):
                gl_6 = -0.5-0.9*(mylogP-1.)
            if((mylogP > 2.0)and(mylogP <= 4.0)):
                gl_6 = -1.4-0.3*(mylogP-2.)
            if(mylogP > 4.0):
                gl_6 = -2.0
    
            # Given P, interpolate gamma_largeq w/ respect to M1 at myM1
            if(myM1 <= 1.2):
                gl=gl_1p2
            if((myM1 > 1.2)and(myM1 <= 3.5)):
                #gl=interpol([gl_1p2,gl_3p5],np.log10([1.2,3.5]),np.log10(myM1))
                gl=np.interp(np.log10(myM1), np.log10([1.2,3.5]), [gl_1p2,gl_3p5])
            if((myM1 > 3.5)and(myM1 <= 6.0)):
                #gl=interpol([gl_3p5,gl_6],alog10([3.5,6.0]),alog10(myM1))
                gl=np.interp(np.log10(myM1), np.log10([3.5,6.0]), [gl_3p5,gl_6])
            if(myM1 > 6.0):
                gl=gl_6
    
    
            # Power-law slope gamma_smallq for M1 < 1.2 Msun and all P Eqn. 13
            gs_1p2=0.3
    
            # Power-law slope gamma_smallq for M1 = 3.5 Msun and various P Eqn. 14
            if(mylogP <= 2.5):
                gs_3p5 = 0.2
            if((mylogP > 2.5)and(mylogP <= 5.5)):
                gs_3p5 = 0.2-0.3*(mylogP-2.5)
            if(mylogP > 5.5):
                gs_3p5 =-0.7-0.2*(mylogP-5.5)
    
            # Power-law slope gamma_smallq for M1 > 6 Msun and various P Eqn. 15
            if(mylogP <= 1.0):
                gs_6 = 0.1
            if((mylogP > 1.0)and(mylogP <= 3.0)):
                gs_6 = 0.1-0.15*(mylogP-1.)
            if((mylogP > 3.0)and(mylogP <= 5.6)):
                gs_6 =-0.2-0.50*(mylogP-3.)
            if(mylogP > 5.6):
                gs_6 =-1.5
    
            # Given P, interpolate gamma_smallq w/ respect to M1 at myM1
            if(myM1 <= 1.2):
                gs=gs_1p2
            if((myM1 > 1.2)and(myM1 <= 3.5)):
                gs = np.interp(np.log10(myM1),np.log10([1.2,3.5]),[gs_1p2,gs_3p5])
            if((myM1 > 3.5)and(myM1 <= 6.0)):
               gs = np.interp(np.log10(myM1),np.log10([3.5,6.0]),[gs_3p5,gs_6])
            if(myM1 > 6.0):
                gs=gs_6
    
    
            # Given Ftwin, gamma_smallq, and gamma_largeq at the specified M1 & P,
            # tabulate the cumulative mass ratio distribution across 0.1 < q < 1.0
            fq=qv**gl                                   # slope across 0.3 < q < 1.0
            fq=fq/idl_tabulate(qv[indlq],fq[indlq])     # normalize to 0.3 < q < 1.0
            fq=fq*(1.-Ftwin)+H*Ftwin                    # add twins
            fq[indsq]=fq[indq0p3]*(qv[indsq]/0.3)**gs   # slope across 0.1 < q < 0.3
            cumfq=np.cumsum(fq)-fq[0]                   # cumulative distribution
            cumfq=cumfq/np.max(cumfq)                   # normalize cumfq(q=1.0) = 1
            cumqdist[:,j,i]=cumfq                       # save to grid
    
    
            # Given M1 and P, q_factor is the ratio of all binaries 0.1 < q < 1.0
            # to those with 0.3 < q < 1.0
            q_factor = idl_tabulate(qv,fq)
    
    
            # Given M1 & P, calculate power-law slope eta of eccentricity dist.
            if(mylogP >= 0.7):
                # For log P > 0.7 use fits in Section 9.2.
                # Power-law slope eta for M1 < 3 Msun and log P > 0.7
                eta_3= 0.6-0.7/(mylogP-0.5)  # Eqn. 17
                # Power-law slope eta for M1 > 7 Msun and log P > 0.7
                eta_7= 0.9-0.2/(mylogP-0.5)  # Eqn. 18
            else:
                # For log P < 0.7, set eta to fitted values at log P = 0.7
                eta_3 = -2.9
                eta_7 = -0.1
    
    
            # Given P, interpolate eta with respect to M1 at myM1
            if(myM1 <= 3.):
                eta = eta_3
            if((myM1 > 3.)and(myM1 <= 7.)):
                eta = np.interp(np.log10(myM1),np.log10([3.,7.]),[eta_3, eta_7])
            if(myM1 > 7.):
                eta = eta_7
    
    
            # Given eta at the specified M1 and P, tabulate eccentricity distribution
            if(10**mylogP <= 2.):
                # For P < 2 days, assume all systems are close to circular
                # For adopted ev (spacing and minimum value), eta = -3.2 satisfies this
                fe=ev**(-3.2)
            else:
                fe=ev**eta
                e_max=1.-(10**mylogP/2.)**(-2./3.) # maximum eccentricity for given P
                ind=np.where(ev >= e_max)
                fe[ind]=0.                       # set dist. = 0 for e > e_max
                # Assume e dist. has power-law slope eta for 0.0 < e / e_max < 0.8 and
                # then linear turnover between 0.8 < e / e_max < 1.0 so that dist.
                # is continuous at e / e_max = 0.8 and zero at e = e_max
                ind=np.where((ev >= 0.8*e_max)&(ev <= 1.0*e_max))
                ind_cont=np.min(ind)-1
                fe[ind] = np.interp(ev[ind],[0.8*e_max,1.0*e_max],[fe[ind_cont],0.])
    
            cumfe=np.cumsum(fe)-fe[0]  # cumulative distribution
            cumfe=cumfe/np.max(cumfe)             # normalize cumfe(e=e_max) = 1
            cumedist[:,j,i]=cumfe              # save to grid
    
    
            # Given constants alpha and DlogP and
            # M1 dependent values flogPle1, flogPeq2p7, and flogPeq5p5,
            # calculate frequency flogP of companions with q > 0.3 per decade
            # of orbital period at given P (Section 9.3 and Eqn. 23)
            if(mylogP <= 1.):
                flogP=flogPle1
            if((mylogP > 1.)and(mylogP <= 2.7-DlogP)):
                flogP=flogPle1+(mylogP-1.)/(1.7-DlogP)*(flogPeq2p7-flogPle1-alpha*DlogP)
            if((mylogP > 2.7-DlogP)and(mylogP <= 2.7+DlogP)):
                flogP=flogPeq2p7+alpha*(mylogP-2.7)
            if((mylogP > 2.7+DlogP)and(mylogP <= 5.5)):
                flogP=flogPeq2p7+alpha*DlogP+(mylogP-2.7-DlogP)/(2.8-DlogP)*(flogPeq5p5-flogPeq2p7-alpha*DlogP)
            if(mylogP > 5.5):
                flogP=flogPeq5p5*np.exp(-0.3*(mylogP-5.5))
    
    
            # Convert frequency of companions with q > 0.3 to frequency of
            # companions with q > 0.1 according to q_factor save to grid
            flogP_sq[j,i]=flogP*q_factor
    
            # Calculate prob. that a companion to M1 with period P is the
            # inner binary.  Currently this is an approximation.
            # 100% for log P < 1.5
            # For log P > 1.5 adopt functional form that reproduces M1 dependent
            # multiplicity statistics in Section 9.4, including a
            # 41% binary star faction (59% single star fraction) for M1 = 1 Msun and
            # 96% binary star fraction (4% single star fraction) for M1 = 28 Msun
            if(mylogP <= 1.5):
                probbin[j,i]=1.0
            else:
                probbin[j,i]=1.0-0.11*(mylogP-1.5)**1.43*(myM1/10.)**0.56
            if(probbin[j,i] <= 0.):
                probbin[j,i]=0.
    
        # Given M1, calculate cumulative binary period distribution
        mycumPbindist=np.cumsum(flogP_sq[:,i]*probbin[:,i]) - flogP_sq[0,i]*probbin[0,i]
        # Normalize so that max(cumPbindist) = total binary star fraction (NOT 1)
        mycumPbindist=mycumPbindist/np.max(mycumPbindist)*idl_tabulate(logPv,flogP_sq[:,i]*probbin[:,i])
        cumPbindist[:,i]=mycumPbindist  # save to grid
    

    """
    #
    # Step 2
    #
    # Implement Monte Carlo method / random number generator to select
    # single stars and binaries from the grids of distributions
    #
    """
    
    
    # Create vector for PRIMARY mass function, which is the mass distribution
    # of single stars and primaries in binaries.
    # This is NOT the IMF, which is the mass distribution of single stars,
    # primaries in binaries, and secondaries in binaries.
    
    # Full primary mass vector across 0.08 < M1 < 150
    M1 = np.linspace(0,150000,150000)*0.001+0.08
    # Slope = -2.3 for M1 > 1 Msun
    fM1=M1**(-2.3)
    # Slope = -1.6 for M1 = 0.5 - 1.0 Msun
    ind=np.where(M1 <= 1.)
    fM1[ind]=M1[ind]**(-1.6)
    # Slope = -0.8 for M1 = 0.15 - 0.5 Msun
    ind=np.where(M1 <= 0.5)
    fM1[ind]=M1[ind]**(-0.8)/0.5**(1.6-0.8)
    # Cumulative primary mass distribution function
    cumfM1=np.cumsum(fM1)-fM1[0]
    cumfM1=cumfM1/np.max(cumfM1)
    

    # Define sampler function
    def parameterSampler(M1min=5.0):

        # Minimum primary mass to generate (must be >0.080 Msun)
        
        # Value of primary mass CDF where M1 = M1min
        cumf_M1min=np.interp(M1min,M1,cumfM1)
        
        # Select primary M1 > M1min from primary mass function
        myM1 = np.interp(cumf_M1min+(1.0-cumf_M1min)*np.random.rand(),cumfM1,M1)
    
        # Find index of M1v that is closest to myM1.
        #     For M1 = 40 - 150 Msun, adopt binary statistics of M1 = 40 Msun.
        #     For M1 = 0.08 - 0.8 Msun, adopt P and e dist of M1 = 0.8Msun,
        #     scale and interpolate the companion frequencies so that the
        #     binary star fraction of M1 = 0.08 Msun primaries is zero,
        #     and truncate the q distribution so that q > q_min = 0.08/M1
        indM1=np.where(abs(myM1-M1v) == min(abs(myM1-M1v)))
        indM1=indM1[0]
    
    
        # Given M1, determine cumulative binary period distribution
        mycumPbindist=(cumPbindist[:,indM1]).flatten
        # If M1 < 0.8 Msun, rescale to appropriate binary star fraction
        if(myM1 <= 0.8):
            mycumPbindist=mycumPbindist()*np.interp(np.log10(myM1),np.log10([0.08,0.8]),[0.0,1.0])
    
    
        # Given M1, determine the binary star fraction
        mybinfrac = np.max(mycumPbindist())
    
    
        # Generate random number myrand between 0 and 1
        myrand=np.random.rand()
    
    
        # If random number < binary star fraction, generate a binary
        if(myrand < mybinfrac):
            # Given myrand, select P and corresponding index in logPv
            mylogP=np.interp(myrand,mycumPbindist(),logPv)
            indlogP=np.where(abs(mylogP-logPv) == min(abs(mylogP-logPv)))
            indlogP=indlogP[0]
    
    
            # Given M1 & P, select e from eccentricity distribution
            mye=np.interp(np.random.rand(),cumedist[:,indlogP,indM1].flatten(),ev)
    
    
            # Given M1 & P, determine mass ratio distribution.
            # If M1 < 0.8 Msun, truncate q distribution and consider
            # only mass ratios q > q_min = 0.08 / M1
            mycumqdist= cumqdist[:,indlogP,indM1].flatten()
            if(myM1 < 0.8):
                q_min=0.08/myM1
                # Calculate cumulative probability at q = q_min
                cum_qmin = np.interp(q_min,qv,mycumqdist)
                # Rescale and renormalize cumulative distribution for q > q_min
                mycumqdist=mycumqdist-cum_qmin
                mycumqdist=mycumqdist/max(mycumqdist)
                # Set probability = 0 where q < q_min
                indq=np.where(qv <= q_min)
                mycumqdist[indq] = 0.0
    
            # Given M1 & P, select q from cumulative mass ratio distribution
            myq=np.interp(np.random.rand(),mycumqdist,qv)
    
    
            # Print M1, q, P & e to file
            myM2 = myq*myM1
            myP = np.power(10, mylogP)
    
        else:
            # If instead random number > binary star fraction, generate single star
    
            # Set dummy variables so that the single star is really run in a wide binary
            myM2 = 0.1
            myP = 1e8 # very wide
            mye = 0
    
        return myM1, myM2, myP, mye

    return parameterSampler



def createParameterGridfile(gridname, nSamples, m1Min=5, m1Max=100):

    # Define the distributions and extract the sampler
    parameterSampler = createParameterDistributionsAndSampler()

    with open(gridname, 'w') as fwrite:

        # Loop over the the number of desired samples, tossing any above the primary mass limit
        nEntries = 0
        while nEntries < nSamples:

            m1, m2, P, e = parameterSampler(m1Min)

            # Need to explicitly toss samples > m1Max
            if m1 > m1Max:
                continue # don't include in overall count
            nEntries += 1    

            fwrite.write('--initial-mass-1 {} --initial-mass-2 {} --orbital-period {} --eccentricity {}\n'.format(m1, m2, P, e))





if __name__ == "__main__":

    # USER SET PARAMETERS
    nSamples=1e5
    m1Min=5 
    m1Max=100

    # LEAVE THESE 
    gridname = 'grid_moedistefano_nSamples{}.txt'.format(int(nSamples))
    createParameterGridfile(gridname=gridname, nSamples=nSamples)
