#include "AIS.h"


// Default constructor
AIS::AIS() {

    m_ObjectId = globalObjectId++;              // unique object id - remains for life of star (even through evolution to other phases)


    m_DrawingFromAISDistributions = false;      // if true we sample from AIS distributions
    m_RandomGaussianDraw          = 99999;
    m_FractionSampled             = 0.0;
    m_FractionExploratory         = 1.0;
    m_CounterDCOsAIS              = 0;          // counts the number of DCOs of interest ("hits") in the exploratory phase AIS
}


/*
 * Initialise the AIS object
 *
 *
 * void Initialise()
 */
void AIS::Initialise() {

    m_DrawingFromAISDistributions = true;                       //  we are drawing from AIS distributions in AIS phase 2 // floor move this to only do once True

    m_RandomGaussianDraw = RAND->RandomInt(m_CovLogA.size());   // the max int (i.e. index of gaussian) we want to sample. -1 since cov_loga also read in empty line  JR: todo: check I fixed the empty line...
}


/*
 * Print the selections of the 'hits' of the exploratory phase of the Adaptive Importance sampling
 *
 * Prints to stdout         JR: todo: do we want to print this to a file to record the setting?
 *
 *
 * void PrintExploratorySettings()
 */
void AIS::PrintExploratorySettings() {

	SAY("");
	SAY(" ------------------------------------------------------------- ");
	SAY(" Running the Adaptive Importance Sampling exploratory phase ");

    switch (_(utils::ToUpper((OPTIONS->AIS_DCOTypeString())).c_str())) {
        case _("ALL") : SAY(" - Selecting all DCOs");         break;
        case _("BBH") : SAY(" - Selecting all BBH mergers");  break;
        case _("BNS") : SAY(" - Selecting all BNS mergers");  break;
        case _("BHNS"): SAY(" - Selecting all BHNS mergers"); break;
        default       : SHOW_WARN(ERROR::INVALID_AIS_DCO_TYPE);     // show warning
    }

    if (OPTIONS->AIS_Hubble())      { SAY(" - Selecting only binaries that merge in Hubble time    "); }
    else                            { SAY(" - Selecting binaries that both do and do not merge in Hubble time   "); }

    if (OPTIONS->AIS_RLOF())        { SAY(" - Excluding binaries that have RLOFSecondaryAfterCEE == 1    "); }
    else                            { SAY(" - Not excluding binaries that have RLOFSecondaryAfterCEE == 1   "); }

    if (OPTIONS->AIS_Pessimistic()) { SAY(" - Selecting only Pessimistic CE binaries    "); }
    else                            { SAY(" - Selecting Optimistic CE binaries   "); }

    SAY(" ------------------------------------------------------------- ");
    SAY("");
}


/*
 * Reads gaussian means and covariances from AISgaussians.txt file and populates the member variable vectors:
 *
 *    m_CovLogA;
 *    m_CovM1;
 *    m_CovQ;
 *    m_MuLogA;
 *    m_MuM1;
 *    m_MuQ;
 *
 * The values stored in these vectors are used to draw random samples for Adaptive Importance Sampling (step 2)
 *
 *
 * void DefineGaussians()
 */
void AIS::DefineGaussians() {

    string filename = "AISgaussians.txt";                                           // JR: todo: put this in constants.h?  or maybe in AIS.h

	if (!utils::FileExists(filename)) {                                             // check gaussians file exists
        SHOW_WARN(ERROR::FILE_DOES_NOT_EXIST, filename);                            // does not exist - show warning
	}
	else {                                                                          // file exists
        std::ifstream infile;
    	infile.open(filename);                                                      // open the file
	    if (infile.fail()) {                                                        // open failed - show warning
            SHOW_WARN(ERROR::FILE_OPEN_ERROR, filename);
	    }
        else {
            int numGaussians = 0;                                                   // number of rows of Gaussians in file
            infile >> numGaussians;
            if (numGaussians > 0) {
                double muM1, muLoga, muQ, covM1, covLogA, covQ;                     // temp storage

                // JR: todo: probably should make this a bit more robust
                while (infile >> muM1) {                                            // read the file
                    infile >> muLoga >> muQ >> covM1 >> covLogA >> covQ;            // file contains means and covariances of Gaussians in 6 columns

                    m_MuM1.push_back(muM1);
                    m_MuLogA.push_back(muLoga);
                    m_MuQ.push_back(muQ);
                    m_CovM1.push_back(covM1);
                    m_CovLogA.push_back(covLogA);
                    m_CovQ.push_back(covQ);
                }
            }

            if (numGaussians <= 0 || (unsigned int)numGaussians != m_MuM1.size()) { // count doesn't match that specified in the file
                SHOW_WARN(ERROR::OUT_OF_BOUNDS, "AIS Gaussians file count error");  // show warning
            }
        }
	}

    SAY(" ------------------------------------------------------------- ");
    SAY(" Started Adaptive Importance Sampling (step 2)");
    SAY(" COMPAS will use an instrumental distribution based on " << m_CovQ.size() << " Gaussian distributions, i.e. 'hits'.  ");
    SAY(" ------------------------------------------------------------- ");
}


/*
 * Update the fraction that should be spent on the AIS exploratory phase
 *
 * Broekgaarden+18, eq 16
 *
 * Modifies class member variable m_FractionExploratory (uses previous value to calculate new value)
 *
 * void UpdateExploratoryPhaseFraction(const int p_PopulationSize)
 *
 * @param   [IN]    p_PopulationSize            The number of binaries in the current population
 */
void AIS::UpdateExploratoryPhaseFraction(const int p_PopulationSize) {

    if (OPTIONS->NBatchesUsed() > 0) {                                                                                  // check user defined number of batches
        double z1             = double(m_CounterDCOsAIS) / double(p_PopulationSize);                                    // the estimated weight of the target population region
        double z2             = 1.0 / (m_FractionExploratory * OPTIONS->NBatchesUsed() * double(OPTIONS->NBinaries())); // estimated rate of unidentified region
        double _1_z1          = 1.0 - z1;
        double sqrt_1_z1      = sqrt(_1_z1);

        double numerator      = z1 * (sqrt_1_z1 - sqrt(z2));                                                            // numerator of Broekgaarden+18, eq 16
        double denominator    = sqrt_1_z1 * (sqrt(z2 * _1_z1) + z1);                                                    // dominator of Broekgaarden+18, eq 16

        m_FractionExploratory = 1.0 - (float(numerator) / denominator);
    }
    else {                                                                                                              // nBatches must be positive
        SHOW_WARN(ERROR::OUT_OF_BOUNDS, "nBatches < 0!");                                                               // show warning
    }
}


/*
 * Determine whether the AIS exploratory phase should stop
 *
 * Updates the exporatory phase fraction (calls UpdateExploratoryPhaseFraction()) if 2 or more DCOs of
 * interest have been found
 *
 * bool ShouldStopExploratoryPhase(const int p_PopulationSize)
 *
 * @param   [IN]    p_PopulationSize            The number of binaries in the current population
 * @return                                      True if the exploratory phase should stop, else false
 */
bool AIS::ShouldStopExploratoryPhase(const int p_PopulationSize) {

    bool shouldStop = false;                                                                                // default is don't stop evolution

    if (m_CounterDCOsAIS >= 2) {                                                                            // if we have found at least 2 hits in total (to overcome possibility of 100% of samples are hits)
        UpdateExploratoryPhaseFraction(p_PopulationSize);                                                   // update fexplAIS to estimate how long we should be spending on sampling from exploratory phase
    }                                                                                                       // Floor: we could update this only every other 10 runs in the future..

    m_FractionSampled = double(p_PopulationSize) / OPTIONS->NBinaries();                                    // calculate fraction so far spent on exploratory phase

    if (utils::Compare(m_FractionExploratory, 1.0) != 0 && m_FractionSampled >= m_FractionExploratory) {    // if the fraction of total samples that we spend is larger than the exploratory fraction...
        shouldStop = true;                                                                                  // ... we should stop and switch to refinement phase
    }

    return shouldStop;
}
