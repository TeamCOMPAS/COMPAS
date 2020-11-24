#!/usr/bin/env python3
import numpy as np
import h5py as h5
import os
import totalMassEvolvedPerZ as MPZ


class COMPASData(object):
    def __init__(
        self,
        path=None,
        fileName="COMPAS_output.h5",
        lazyData=True,
        Mlower=None,
        Mupper=None,
        binaryFraction=None,
    ):
        self.path = path
        self.fileName = fileName
        if self.path is None:
            print("Just to double check you create instance of ClassCOMPAS without path/Data")
        elif not os.path.isfile(path + fileName):
            raise ValueError(
                "h5 file not found. Wrong path given? %s" % (path + fileName)
            )
        elif os.path.isfile(path + fileName):
            pass

        # Crucial values to be able to calculate MSSFR
        self.metallicityGrid = None
        self.metallicitySystems = None
        self.delayTimes = None  # Myr
        # Crucial values I need for selection effects
        self.mass1 = None  # Msun
        self.mass2 = None  # Msun
        self.DCOmask = None
        self.allTypesMask = None
        self.BBHmask = None
        self.DNSmask = None
        self.BHNSmask = None
        self.initialZ = None

        # Additional arrays that might be nice to store
        # to more quickly make some plots.
        # If you need more memory might help a tiny bit to not do
        self.lazyData = lazyData
        self.mChirp = None  # Msun
        self.q = None
        self.optimisticmask = None

        # Needed to recover true solar mass evolved
        self.Mlower = Mlower  # Msun
        self.Mupper = Mupper  # Msun
        self.binaryFraction = binaryFraction
        self.totalMassEvolvedPerZ = None  # Msun

        print("ClassCOMPAS: Remember to self.setGridAndMassEvolved() [optional]")
        print("                   then  self.setCOMPASDCOmask()")
        print("                   then  self.setCOMPASData()")

    def setCOMPASDCOmask(
        self, types="BBH", withinHubbleTime=True, pessimistic=True, noRLOFafterCEE=True
    ):
        # By default, we mask for BBHs that merge within a Hubble time, assumming
        # the pessimistic CEE prescription (HG donors cannot survive a CEE) and
        # not allowing immediate RLOF post-CEE
        Data = h5.File(self.path + self.fileName, "r")
        fDCO = Data["DoubleCompactObjects"]
        fCEE = Data["CommonEnvelopes"]

        # Masks for DCO type
        maskBBH = (fDCO["Stellar_Type_1"][()] == 14) & (fDCO["Stellar_Type_2"][()] == 14)
        maskDNS = (fDCO["Stellar_Type_1"][()] == 13) & (fDCO["Stellar_Type_2"][()] == 13)
        maskBHNS = ((fDCO["Stellar_Type_1"][()] == 14) & (fDCO["Stellar_Type_2"][()] == 13)) \
                    |((fDCO["Stellar_Type_1"][()] == 13) & (fDCO["Stellar_Type_2"][()] == 14))
        maskAllTypes = maskBBH | maskDNS | maskBHNS

        if types == "BBH":
            maskTypes = maskBBH
        elif types == "BNS":
            maskTypes = maskDNS
        elif types == "BHNS":
            maskTypes = maskBHNS
        elif types == "all":
            maskTypes = maskAllTypes
        else:
            raise ValueError("type=%s not one of 'BBH', 'BNS', 'BHNS', 'all'" % (types))

        # Mask DCOs merging within Hubble time
        if withinHubbleTime:
            maskHubble = fDCO["Merges_Hubble_Time"][()] == True
        else:
            maskHubble = np.ones(len(fDCO["Merges_Hubble_Time"][()]), dtype=bool)

        # Masks related to CEEs
        CEEseeds = fCEE["SEED"][()]
        DCOseeds = fDCO["SEED"][()]

        if (pessimistic or noRLOFafterCEE):
            maskDCOCEEs = np.in1d(CEEseeds, DCOseeds) # Mask for CEEs involved in forming DCOs
            DCOCEEseeds = CEEseeds[maskDCOCEEs] # Seeds of CEEs involved in forming DCO

        # Mask for DCOs formed assumming pessimistic CEE
        if pessimistic:
            optimisticFlagForDCOCEEs = fCEE["Optimistic_CE"][()][maskDCOCEEs]
            # Seeds of DCOs that have optimistic CEE
            optimisticCEEseeds = np.unique(DCOCEEseeds[optimisticFlagForDCOCEEs])
            maskOptimistic = np.in1d(DCOseeds,optimisticCEEseeds)
            maskPessimistic = np.logical_not(maskOptimistic)
        else:
            maskPessimistic = np.ones(len(fDCO["ID"][()]), dtype=bool)

        # Mask for DCOs formed without RLOF immediately after CEE
        if noRLOFafterCEE:
            immediateRLOFflagforDCOCEEs = fCEE["Immediate_RLOF>CE"][()][maskDCOCEEs]
            # Seeds of DCOs that have immediate RLOF post-CEE
            immediateRLOFCEEseeds = np.unique(DCOCEEseeds[immediateRLOFflagforDCOCEEs])
            maskRLOFafterCEE = np.in1d(DCOseeds,immediateRLOFCEEseeds)
            maskNoRLOFafterCEE = np.logical_not(maskRLOFafterCEE)
        else:
            maskNoRLOFafterCEE = np.ones(len(fDCO["ID"][()]), dtype=bool)

        # Combine all the masks
        self.DCOmask = maskTypes & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.BBHmask = maskBBH & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.DNSmask = maskDNS & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.BHNSmask = maskBHNS & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.allTypesMask = maskAllTypes & maskHubble & maskPessimistic & maskNoRLOFafterCEE
        self.optimisticmask = maskPessimistic
        Data.close()

    def setGridAndMassEvolved(self):
        # The COMPAS simulation does not evolve all stars
        # give me the correction factor for the total mass evolved
        # I assume each metallicity has the same limits, and does correction
        # factor, but the total mass evolved might be different.
        # This does not change when we change types and other masks this is
        # general to the entire simulation so calculate once
        _, self.totalMassEvolvedPerZ = MPZ.totalMassEvolvedPerZ(
            path=self.path,
            fileName=self.fileName,
            Mlower=self.Mlower,
            Mupper=self.Mupper,
            binaryFraction=self.binaryFraction,
        )
        # Want to recover entire metallicity grid, assume that every metallicity
        # evolved shows in all systems again should not change within same run
        # so dont redo if we reset the data
        Data = h5.File(self.path + self.fileName, "r")
        if self.initialZ is None:
            self.initialZ = Data["SystemParameters"]["Metallicity@ZAMS_1"][()]
        self.metallicityGrid = np.unique(self.initialZ)
        Data.close()

    def setCOMPASData(self):
        Data = h5.File(self.path + self.fileName, "r")
        fDCO = Data["DoubleCompactObjects"]
        # sorry not the prettiest line is a boolean slice of seeds
        # this only works because seeds in systems file and DCO file are printed
        # in same order

        # Get metallicity grid of DCOs
        self.seedsDCO = fDCO["SEED"][()][self.DCOmask]
        initialSeeds = Data["SystemParameters"]["SEED"][()]
        if self.initialZ is None:
            self.initialZ = Data["SystemParameters"]["Metallicity@ZAMS_1"][()]
        maskMetallicity = np.in1d(initialSeeds, self.seedsDCO)
        self.metallicitySystems = self.initialZ[maskMetallicity]

        self.delayTimes = np.add(
            fDCO["Time"][()][self.DCOmask], fDCO["Coalescence_Time"][()][self.DCOmask]
        )
        self.mass1 = fDCO["Mass_1"][()][self.DCOmask]
        self.mass2 = fDCO["Mass_2"][()][self.DCOmask]

        # Stuff of data I dont need for integral
        # but I might be to laze to read in myself
        # and often use. Might turn it of for memory efficiency
        if self.lazyData:
            self.q = np.divide(self.mass2, self.mass1)
            boolq = self.mass2 > self.mass1
            self.q[boolq] = np.divide(self.mass1[boolq], self.mass2[boolq])
            self.mChirp = np.divide(
                (np.multiply(self.mass2, self.mass1) ** (3.0 / 5.0)),
                (np.add(self.mass2, self.mass1) ** (1.0 / 5.0)),
            )
            self.Hubble = fDCO["Merges_Hubble_Time"][...].squeeze()[self.DCOmask]

        Data.close()

    def recalculateTrueSolarMassEvolved(self, Mlower, Mupper, binaryFraction):
        # Possibility to test assumptions of True solar mass evolved
        self.Mlower = Mlower
        self.Mupper = Mupper
        self.binaryFraction = binaryFraction
        _, self.totalMassEvolvedPerZ = MPZ.totalMassEvolvedPerZ(
            pathCOMPASh5=self.path,
            Mlower=self.Mlower,
            Mupper=self.Mupper,
            binaryFraction=self.binaryFraction,
        )
