#ifndef __AIS_h__
#define __AIS_h__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>

#include <boost/algorithm/string.hpp>

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "Options.h"
#include "Errors.h"
#include "Rand.h"


class AIS {

public:

	AIS();


    // object identifiers - all classes have these
	OBJECT_ID     ObjectId()                            { return m_ObjectId; }
    OBJECT_TYPE   ObjectType()                          { return OBJECT_TYPE::AIS; }
    STELLAR_TYPE  StellarType()                         { return STELLAR_TYPE::NONE; }


    // getters - alphabetically
	double      CovLogA()  const                       { return m_CovLogA[m_RandomGaussianDraw]; }
	double      CovM1() const                          { return m_CovM1[m_RandomGaussianDraw]; }
	double      CovQ() const                           { return m_CovQ[m_RandomGaussianDraw]; }

	bool        DrawingFromAISDistributions() const    { return m_DrawingFromAISDistributions; }

	double      MuLogA() const                         { return m_MuLogA[m_RandomGaussianDraw]; }
	double      MuM1() const                           { return m_MuM1[m_RandomGaussianDraw]; }
	double      MuQ() const                            { return m_MuQ[m_RandomGaussianDraw]; }


    // member functions - alphabetically (sort of)

    /*
     * Calculate if we have a Double Compact Object of interest
     *
     * If the DCO is of interest, update AIS internal counter m_CounterDCOsAIS and return true,
     * otherwise return false without updating the internal counter
     *
     * This function uses a template parameter and so is defined in the header file so that it
     * can be declared before the BaseBinaryStar class (it is used inside the BaseBinaryStar class)
     *
     *
     * bool CalculateDCOHit(const T* const p_Binary)
     *
     * @param   [IN]    p_Binary                    The binary star to interrogate
     * @return                                      True if DCO is of interest, else false
     */
    template <class T>
    bool CalculateDCOHit(const T* const p_Binary) {

        int DCOhit = 0;                                                                                 // default is no hit

        if ((utils::Equals(OPTIONS->AIS_DCOTypeString(), "ALL" )                            ||          // check all DCOs
            (utils::Equals(OPTIONS->AIS_DCOTypeString(), "BBH" ) && p_Binary->IsBHandBH())  ||          // check only BH + BH
            (utils::Equals(OPTIONS->AIS_DCOTypeString(), "BNS" ) && p_Binary->IsNSandNS())  ||          // check only NS + NS
            (utils::Equals(OPTIONS->AIS_DCOTypeString(), "BHNS") && p_Binary->IsNSandBH()))) {          // check only NS + NS

            DCOhit = 1;                                                                                 // assume hit, unless...

            if ((OPTIONS->AIS_Hubble()      && !p_Binary->MergesInHubbleTime())     ||                  // exclude DCOs that do not merge in Hubble if user wants
                (OPTIONS->AIS_RLOF()        &&  p_Binary->RLOFSecondaryPostCEE())   ||                  // exclude DCOs that have RLOFafterSecondaryAfterCE if user wants
                (OPTIONS->AIS_Pessimistic() &&  p_Binary->OptimisticCommonEnvelope())) DCOhit = 0;      // exclude DCOs that have optimisticCEFlag if user wants
        }

        m_CounterDCOsAIS += DCOhit;                                                                     // update iternal counter

        return (DCOhit == 1);
    }

    void        DefineGaussians();
    void        Initialise();
    void        PrintExploratorySettings();
    bool        ShouldStopExploratoryPhase(const int p_PopulationSize);
    void        UpdateExploratoryPhaseFraction(const int p_PopulationSize);


protected:

    OBJECT_ID    m_ObjectId;                                                                            // instantiated object's unique object id


    // member variables - alphabetical
	DBL_VECTOR m_CovLogA;
	DBL_VECTOR m_CovM1;
	DBL_VECTOR m_CovQ;

	DBL_VECTOR m_MuLogA;
	DBL_VECTOR m_MuM1;
	DBL_VECTOR m_MuQ;


    // The following are AIS constants:
    //
    // IMF variables
    // Taken from Kroupa 2001 https://arxiv.org/pdf/astro-ph/0009005.pdf
    //
    // Some of these are just plain constants, but other need to be calculated
    // (the one that use pow()).  They are calculated once in the AIS constructor
    // so that we don't waste time calculating them every time we use them.
    //
    // Some compilers let us declare them all as constexpr, but some compilers don't
    // allow math functions (like pow()) in constexpr statements (the math functions
    // delivered with some compilers aren't declared as constexpr themselves), hence
    // the need to calculate these in Initialise().


    // Kroupa IMF is a broken power law with three slopes
    double kroupaPower1;
    double kroupaPower2;
    double kroupaPower3;

    // Often require the power law exponent plus one
    double kroupaPowerPlus1_1;
    double kroupaPowerPlus1_2;
    double kroupaPowerPlus1_3;

    double oneOverKroupaPower1Plus1;
    double oneOverKroupaPower2Plus1;
    double oneOverKroupaPower3Plus1;

    // There are two breaks in the Kroupa power law -- they occur here (in solar masses)
    double kroupaBreak1;
    double kroupaBreak2;

    // Some values that can be calculated at compile time
    // Works with gcc c++11 - may not work with other compilers
    double kroupaBreak1_Plus1_1;
    double kroupaBreak1_Plus1_2;
    double kroupaBreak1_Power1_2;

    double kroupaBreak2_Plus1_2;
    double kroupaBreak2_Plus1_3;
    double kroupaBreak2_Power2_3;


    // member functions - alphabetically
    int         m_CounterDCOsAIS;                                                                       // CounterrDCOsAIS counts the number of DCOs of interest ("hits") in the exploratory phase of Adaptive Importance Sampling
	bool        m_DrawingFromAISDistributions;                                                          // if true we sample from AIS distributions
    double      m_FractionExploratory;                                                                  // fraction of samples that should be spend on exploratory phase
    int         m_RandomGaussianDraw;                                                                   // random number that determines from which Gaussian we draw.
    double      m_FractionSampled;                                                                      // fraction of total samples that we have already sampled in exploratory phase

};

#endif  // __AIS_h__
