import numpy as np               #for handling arrays
import h5py as h5                #for reading the COMPAS data
import matplotlib.pyplot as plt  #for plotting
import math

# Get some data to plot
pathToData = 'COMPAS_Output.h5'
Data  = h5.File(pathToData)

#print(list(Data.keys()))
#['CommonEnvelopes', 'DoubleCompactObjects', 'Supernovae', 'SystemParameters']

SPs = Data['SystemParameters']
# ['CE_Alpha', 'Eccentricity@ZAMS', 'Error', 'ID', 'LBV_Multiplier', 'Mass@ZAMS_1', 'Mass@ZAMS_2', 'Merger', 'Merger_At_Birth', 'Metallicity@ZAMS_1', 'Metallicity@ZAMS_2', 'Omega@ZAMS_1', 'Omega@ZAMS_2', 'SEED', 'SN_Kick_Magnitude_Random_Number_1', 'SN_Kick_Magnitude_Random_Number_2', 'SN_Kick_Mean_Anomaly_1', 'SN_Kick_Mean_Anomaly_2', 'SN_Kick_Phi_1', 'SN_Kick_Phi_2', 'SN_Kick_Theta_1', 'SN_Kick_Theta_2', 'Separation@ZAMS', 'Sigma_Kick_CCSN_BH', 'Sigma_Kick_CCSN_NS', 'Sigma_Kick_ECSN', 'Sigma_Kick_USSN', 'Stellar_Type@ZAMS_1', 'Stellar_Type@ZAMS_2', 'Stellar_Type_1', 'Stellar_Type_2', 'Unbound', 'WR_Multiplier']
CEs = Data['CommonEnvelopes']
# ['BE_Fixed_1', 'BE_Fixed_2', 'BE_Kruckow_1', 'BE_Kruckow_2', 'BE_Loveridge_1', 'BE_Loveridge_2', 'BE_Loveridge_Winds_1', 'BE_Loveridge_Winds_2', 'BE_Nanjing_1', 'BE_Nanjing_2', 'Binding_Energy<CE_1', 'Binding_Energy<CE_2', 'CE_Event_Count', 'Core_Mass_1', 'Core_Mass_2', 'Double_Core_CE', 'Eccentricity<CE', 'Eccentricity>CE', 'ID', 'Immediate_RLOF>CE', 'Kruckow_1', 'Kruckow_2', 'Lambda@CE_1', 'Lambda@CE_2', 'Lambda_Fixed_1', 'Lambda_Fixed_2', 'Lambda_Nanjing_1', 'Lambda_Nanjing_2', 'Loveridge_1', 'Loveridge_2', 'Loveridge_Winds_1', 'Loveridge_Winds_2', 'Luminosity<CE_1', 'Luminosity<CE_2', 'MT_History', 'Mass_1<CE', 'Mass_2<CE', 'Mass_Env_1', 'Mass_Env_2', 'Merger', 'Optimistic_CE', 'RLOF_1', 'RLOF_2', 'Radius_1<CE', 'Radius_1>CE', 'Radius_2<CE', 'Radius_2>CE', 'RocheLobe_1<CE', 'RocheLobe_1>CE', 'RocheLobe_2<CE', 'RocheLobe_2>CE', 'SEED', 'Separation<CE', 'Separation>CE', 'Simultaneous_RLOF', 'Stellar_Type_1', 'Stellar_Type_1<CE', 'Stellar_Type_2', 'Stellar_Type_2<CE', 'Tau_Circ', 'Tau_Dynamical<CE_1', 'Tau_Dynamical<CE_2', 'Tau_Nuclear<CE_1', 'Tau_Nuclear<CE_2', 'Tau_Radial<CE_1', 'Tau_Radial<CE_2', 'Tau_Sync', 'Tau_Thermal<CE_1', 'Tau_Thermal<CE_2', 'Teff<CE_1', 'Teff<CE_2', 'Time', 'Zeta_RLOF_Analytic', 'Zeta_Star_Compare']
SNe = Data['Supernovae']
# ['Applied_Kick_Velocity_SN', 'Drawn_Kick_Velocity_SN', 'Eccentricity', 'Eccentricity<2ndSN', 'Experienced_RLOF_SN', 'Fallback_Fraction_SN', 'Hydrogen_Poor_SN', 'Hydrogen_Rich_SN', 'ID', 'Kick_Velocity(uK)', 'Mass_CO_Core@CO_SN', 'Mass_CP', 'Mass_Core@CO_SN', 'Mass_He_Core@CO_SN', 'Mass_SN', 'Mass_Total@CO_SN', 'Orb_Velocity<2ndSN', 'Runaway_CP', 'SEED', 'SN_Kick_Phi_SN', 'SN_Kick_Theta_SN', 'SN_Type_SN', 'Separation', 'Separation<2ndSN', 'Stellar_Type_Prev_CP', 'Stellar_Type_Prev_SN', 'Stellar_Type_SN', 'Supernova_State', 'Systemic_Velocity', 'Time', 'True_Anomaly(psi)_SN', 'Unbound']
DCs = Data['DoubleCompactObjects']
# ['Coalescence_Time', 'Eccentricity@DCO', 'ID', 'MT_Case_1', 'MT_Case_2', 'Mass_1', 'Mass_2', 'Merges_Hubble_Time', 'Recycled_NS_1', 'Recycled_NS_2', 'SEED', 'Separation@DCO', 'Stellar_Type_1', 'Stellar_Type_2', 'Time']

##########

### Want to find relative rates of production channels of isolated pulsars. These are:
# A) 1st SN disruptive, exploder becomes isolated pulsar (only need 1st SN from SNe)
# B) 1st SN disruptive, companion later SNs into pulsar (need 2nd SN from DCs)
# C) 2nd SN disruptive, exploder becomes isolated pulsar (need 2nd SN from SNe)
# D) 2nd SN disruptive, companion was a pulsar now isolated (need 2nd SN from SNe)
# White dwarf accretion (can I even find this?)
# Would also like to know fractions from ECSN, USSN, and CCSN

spSeeds = SPs['SEED'][()]
ceSeeds = CEs['SEED'][()]
snSeeds = SNe['SEED'][()]
dcSeeds = DCs['SEED'][()]

uniqSNseeds = set(snSeeds)
repeatSNseeds = [ seed for seed in uniqSNseeds if np.count_nonzero(snSeeds == seed) > 1 ]

# Need to get seeds of 1st vs 2nd SN in SNe
# use companion types

snCompType = SNe['Stellar_Type_Prev_CP'][()]
snRemnType = SNe['Stellar_Type_SN'][()]
snIsUnbound = SNe['Unbound'][()]

snMask_NSRemn = snRemnType == 13
snMask_BHRemn = snRemnType == 14
snMask_NSComp = snCompType == 13
snMask_Unbound = snIsUnbound == True
snMask_1stSN = snCompType <= 12   # 13 = NS, 14 = BH, 15 = MR
snMask_2ndSN = snCompType  > 12   # 13 = NS, 14 = BH, 15 = MR
snMask_Fucked = snCompType == 15  # 13 = NS, 14 = BH, 15 = MR   # TODO

# Long logic for Case B pulsars
# To get the Case B systems is tricky - 
# These are systems which disrupted at the 1st SN, but it is the companion that later SNs into a pulsar. 
# Because they are disrutped at the first SN, they do not appear in SNe, but the seed will be in DCOs. 
# To find these systems, look at DCO systems which have at least 1 NS, and find their seed in SNe. 
# Isolate for only first SN disrupters in SNe. 
#   If there are 2 NSs at the end, then this system passes. 
#   If there's only 1, need to make sure that the 1st Supernova produces a BH. 

dcType1 = DCs['Stellar_Type_1'][()]
dcType2 = DCs['Stellar_Type_2'][()]

dcSeeds2NS = dcSeeds[(dcType1 == 13) & (dcType2 == 13)]
dcSeeds1NS = dcSeeds[dcType1 != dcType2]     # If not equal and only allowing BHs and NSs, then this must be a BHNS

snMask_dcSeeds2NS = [seed in dcSeeds2NS for seed in snSeeds]
snMask_dcSeeds1NS = [seed in dcSeeds1NS for seed in snSeeds]




###############################################
### Collect seeds for pulsars by Formation Case

# Case A) 1st SN disruptive, exploder becomes isolated pulsar (only need 1st SN from SNe)
seedsCaseA     = snSeeds[ snMask_Unbound & snMask_1stSN & snMask_NSRemn ]

# Case B) 1st SN disruptive, companion later SNs into pulsar (need 2nd SN from DCs)
seedsCaseB_2NS = snSeeds[ snMask_Unbound & snMask_1stSN & snMask_dcSeeds2NS ]
seedsCaseB_1NS = snSeeds[ snMask_Unbound & snMask_1stSN & snMask_dcSeeds1NS & snMask_BHRemn ]
seedsCaseB = np.concatenate((seedsCaseB_2NS, seedsCaseB_1NS))

# Case C) 2nd SN disruptive, exploder becomes isolated pulsar (need 2nd SN from SNe)
seedsCaseC     = snSeeds[ snMask_Unbound & snMask_2ndSN & snMask_NSRemn ]

# Case D) 2nd SN disruptive, companion was a pulsar now isolated (need 2nd SN from SNe)
seedsCaseD     = snSeeds[ snMask_Unbound & snMask_2ndSN & snMask_NSComp ]

print("\nCase A (1st SN disruptive, exploder becomes isolated pulsar) seeds: ")
print(seedsCaseA)
print("\nCase B (1st SN disruptive, companion later SNs into pulsar) seeds: ")
print(seedsCaseB)
print("\nCase C (2nd SN disruptive, exploder becomes isolated pulsar) seeds: ")
print(seedsCaseC)
print("\nCase D (2nd SN disruptive, companion was a pulsar now isolated) seeds: ")
print(seedsCaseD)

nCaseA = len(seedsCaseA)
nCaseB = len(seedsCaseB)
nCaseC = len(seedsCaseC)
nCaseD = len(seedsCaseD)
#nTotal = nCaseA + nCaseB + nCaseC + nCaseD     # this is a shit estimate
nBinaries = len(spSeeds) # total number of input systems
nTotal = 2*nBinaries


fCaseA = nCaseA / nTotal
fCaseB = nCaseB / nTotal
fCaseC = nCaseC / nTotal
fCaseD = nCaseD / nTotal

print()
print(nCaseA, "=> ", round(fCaseA*100, 2), "%")  
print(nCaseB, "=> ", round(fCaseB*100, 2), "%")
print(nCaseC, "=> ", round(fCaseC*100, 2), "%")
print(nCaseD, "=> ", round(fCaseD*100, 2), "%")
print(nTotal, " Total")
print(nBinaries, " Binaries")

