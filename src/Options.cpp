#include "Options.h"
#include "changelog.h"

Options* Options::m_Instance = nullptr;

namespace po  = boost::program_options;
namespace cls = po::command_line_style;
namespace fs  = boost::filesystem;


// this is required to set default value for boost program options of type vector<string>
namespace std
{
  std::ostream& operator<<(std::ostream &os, const std::vector<string> &vec) {    
    for (auto item : vec) os << item << " ";
    return os; 
  }
} 


Options* Options::Instance() {
    if (!m_Instance) {
        m_Instance = new Options();
    }
    return m_Instance;
}


bool Options::ProgramOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription) {

    bool ok = true;                             // status - unless a problem occurs

    // create default strings for vector<string> types (too hard to do inline)

    std::ostringstream ss;

    // debug classes
    string defaultDebugClasses;
    ss << "";
    for (auto debugClass = p_Options->m_DebugClasses.begin(); debugClass != p_Options->m_DebugClasses.end(); ++debugClass) ss << *debugClass << ",";
    defaultDebugClasses = ss.str();
    if (defaultDebugClasses.size() > 0) defaultDebugClasses.erase(defaultDebugClasses.size() - 1);

    // log classes
    string defaultLogClasses;
    ss << "";
    for (auto logClass = p_Options->m_LogClasses.begin(); logClass != p_Options->m_LogClasses.end(); ++logClass) ss << *logClass << ",";
    defaultLogClasses = ss.str();
    if (defaultLogClasses.size() > 0) defaultLogClasses.erase(defaultLogClasses.size() - 1);


    // add options

    try {

        p_OptionsDescription->add_options()     // begin the list of options to be added - boost syntactic sugar

        // there is no good way of formatting these - the boost syntax doesn't help that much
        // there is just a boatload of options, so this function (and similar functions) are just going to be long...


        // switches

        (
            "help,h",                                                      
            po::bool_switch(), "Print this help message"
        )
        (
            "version,v",                                                   
            po::bool_switch(), "Print COMPAS version string"
        )


        // boolean options - alphabetically

        // Floor
        /*
        (
            "AIS-exploratory-phase",                                       
            po::value<bool>(&p_Options->m_AISexploratoryPhase)->default_value(p_Options->m_AISexploratoryPhase)->implicit_value(true),                                                            
            ("Run exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISexploratoryPhase ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-Hubble",                                                  
            po::value<bool>(&p_Options->m_AIShubble)->default_value(p_Options->m_AIShubble)->implicit_value(true),                                                                                
            ("Excluding not in Hubble time mergers selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AIShubble ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-Pessimistic",                                             
            po::value<bool>(&p_Options->m_AISpessimistic)->default_value(p_Options->m_AISpessimistic)->implicit_value(true),                                                                      
            ("Optimistic or Pessimistic selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISpessimistic ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-refinement-phase",                                        
            po::value<bool>(&p_Options->m_AISrefinementPhase)->default_value(p_Options->m_AISrefinementPhase)->implicit_value(true),                                                              
            ("Run main sampling phase (step2) of STROOPWAFEL (default = " + std::string(p_Options->m_AISrefinementPhase ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "AIS-RLOF",                                                    
            po::value<bool>(&p_Options->m_AISrlof)->default_value(p_Options->m_AISrlof)->implicit_value(true),                                                                                    
            ("RLOFSecondaryZAMS selection in exploratory phase of STROOPWAFEL (default = " + std::string(p_Options->m_AISrlof ? "TRUE" : "FALSE") + ")").c_str()
       )
        */

        (
            "BSEswitchLog",                                                
            po::value<bool>(&p_Options->m_BSEswitchLog)->default_value(p_Options->m_BSEswitchLog)->implicit_value(true),                                                                          
            ("Print BSE switch log to file (default = " + std::string(p_Options->m_BSEswitchLog ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "debug-to-file",                                               
            po::value<bool>(&p_Options->m_DebugToFile)->default_value(p_Options->m_DebugToFile)->implicit_value(true),                                                                            
            ("Write debug statements to file (default = " + std::string(p_Options->m_DebugToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "detailedOutput",                                              
            po::value<bool>(&p_Options->m_DetailedOutput)->default_value(p_Options->m_DetailedOutput)->implicit_value(true),                                                                      
            ("Print detailed output to file (default = " + std::string(p_Options->m_DetailedOutput ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "enable-warnings",                                             
            po::value<bool>(&p_Options->m_EnableWarnings)->default_value(p_Options->m_EnableWarnings)->implicit_value(true),                                                                      
            ("Display warning messages to stdout (default = " + std::string(p_Options->m_EnableWarnings ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "errors-to-file",                                              
            po::value<bool>(&p_Options->m_ErrorsToFile)->default_value(p_Options->m_ErrorsToFile)->implicit_value(true),                                                                          
            ("Write error messages to file (default = " + std::string(p_Options->m_ErrorsToFile ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "populationDataPrinting",                                      
            po::value<bool>(&p_Options->m_PopulationDataPrinting)->default_value(p_Options->m_PopulationDataPrinting)->implicit_value(true),                                                      
            ("Print details of population (default = " + std::string(p_Options->m_PopulationDataPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "quiet",                                                       
            po::value<bool>(&p_Options->m_Quiet)->default_value(p_Options->m_Quiet)->implicit_value(true),                                                                                        
            ("Suppress printing (default = " + std::string(p_Options->m_Quiet ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "RLOFPrinting",                                                
            po::value<bool>(&p_Options->m_RlofPrinting)->default_value(p_Options->m_RlofPrinting)->implicit_value(true),                                                                          
            ("Enable output parameters before/after RLOF (default = " + std::string(p_Options->m_RlofPrinting ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "single-star",                                                 
            po::value<bool>(&p_Options->m_SingleStar)->default_value(p_Options->m_SingleStar)->implicit_value(true),                                                                              
            ("Evolve single star(s) (default = " + std::string(p_Options->m_SingleStar ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "SSEswitchLog",                                                
            po::value<bool>(&p_Options->m_SSEswitchLog)->default_value(p_Options->m_SSEswitchLog)->implicit_value(true),                                                                          
            ("Print SSE switch log to file (default = " + std::string(p_Options->m_SSEswitchLog ? "TRUE" : "FALSE") + ")").c_str()
        )


        // numerical options - alphabetically grouped by type 

        // int

        (
            "debug-level",                                                 
            po::value<int>(&p_Options->m_DebugLevel)->default_value(p_Options->m_DebugLevel),                                                                                                     
            ("Determines which print statements are displayed for debugging (default = " + std::to_string(p_Options->m_DebugLevel) + ")").c_str()
        )
        (
            "log-level",                                                   
            po::value<int>(&p_Options->m_LogLevel)->default_value(p_Options->m_LogLevel),                                                                                                         
            ("Determines which print statements are included in the logfile (default = " + std::to_string(p_Options->m_LogLevel) + ")").c_str()
        )

        // Floor
        /*
        (
            "nbatches-used",                                               
            po::value<int>(&p_Options->m_nBatchesUsed)->default_value(p_Options->m_nBatchesUsed),                                                                                                 
            ("Number of batches used, for STROOPWAFEL (AIS), -1 = not required (default = " + std::to_string(p_Options->m_nBatchesUsed) + ")").c_str()
        )
        */

        (
            "number-of-binaries,n",                                        
            po::value<int>(&p_Options->m_nBinaries)->default_value(p_Options->m_nBinaries),                                                                                                       
            ("Specify the number of binaries to simulate (default = " + std::to_string(p_Options->m_nBinaries) + ")").c_str()
        )


        // double

        // Floor
        /*
        (
            "kappa-gaussians",                                             
            po::value<double>(&p_Options->m_KappaGaussians)->default_value(p_Options->m_KappaGaussians),                                                                                          
            ("Scaling factor for the width of the Gaussian distributions in STROOPWAFEL main sampling phase (default = " + std::to_string(p_Options->m_KappaGaussians) + ")").c_str()
        )
        */


        // string options - alphabetically

        // Floor
        /*
        (
            "AIS-DCOtype",                                                 
            po::value<string>(&p_Options->m_AISDCOtypeString)->default_value(p_Options->m_AISDCOtypeString),                                                                                      
            ("DCO type selection in exploratory phase of STROOPWAFEL, (options: ALL, BBH, BNS or BHNS), default = " + p_Options->m_AISDCOtypeString + ")").c_str()
        )
        */

        (
            "grid",                                                        
            po::value<string>(&p_Options->m_GridFilename)->default_value(p_Options->m_GridFilename)->implicit_value(""),                                                                      
            ("Grid filename (default = " + p_Options->m_GridFilename + ")").c_str()
        )

        // Serena
        /*
        (
            "logfile-BSE-be-binaries",                                     
            po::value<string>(&p_Options->m_LogfileBSEBeBinaries)->default_value(p_Options->m_LogfileBSEBeBinaries),                                                                              
            ("Filename for BSE Be Binaries logfile (default = " + p_Options->m_LogfileBSEBeBinaries + ")").c_str()
        )
        */

        (
            "logfile-BSE-rlof-parameters",                                 
            po::value<string>(&p_Options->m_LogfileBSERLOFParameters)->default_value(p_Options->m_LogfileBSERLOFParameters),                                                                      
            ("Filename for BSE RLOF Parameters logfile ( default = " + p_Options->m_LogfileBSERLOFParameters + ")").c_str()
        )
        (
            "logfile-BSE-common-envelopes",                                
            po::value<string>(&p_Options->m_LogfileBSECommonEnvelopes)->default_value(p_Options->m_LogfileBSECommonEnvelopes),                                                                    
            ("Filename for BSE Common Envelopes logfile (default = " + p_Options->m_LogfileBSECommonEnvelopes + ")").c_str()
        )
        (
            "logfile-BSE-detailed-output",                                 
            po::value<string>(&p_Options->m_LogfileBSEDetailedOutput)->default_value(p_Options->m_LogfileBSEDetailedOutput),                                                                      
            ("Filename for BSE Detailed Output logfile (default = " + p_Options->m_LogfileBSEDetailedOutput + ")").c_str()
        )
        (
            "logfile-BSE-double-compact-objects",                          
            po::value<string>(&p_Options->m_LogfileBSEDoubleCompactObjects)->default_value(p_Options->m_LogfileBSEDoubleCompactObjects),                                                          
            ("Filename for BSE Double Compact Objects logfile (default = " + p_Options->m_LogfileBSEDoubleCompactObjects + ")").c_str()
        )
        (
            "logfile-BSE-pulsar-evolution",                                
            po::value<string>(&p_Options->m_LogfileBSEPulsarEvolution)->default_value(p_Options->m_LogfileBSEPulsarEvolution),                                                                    
            ("Filename for BSE Pulsar Evolution logfile (default = " + p_Options->m_LogfileBSEPulsarEvolution + ")").c_str()
        )
        (
            "logfile-BSE-supernovae",                                      
            po::value<string>(&p_Options->m_LogfileBSESupernovae)->default_value(p_Options->m_LogfileBSESupernovae),                                                                              
            ("Filename for BSE Supernovae logfile (default = " + p_Options->m_LogfileBSESupernovae + ")").c_str()
        )
        (
            "logfile-BSE-switch-log",                                      
            po::value<string>(&p_Options->m_LogfileBSESwitchLog)->default_value(p_Options->m_LogfileBSESwitchLog),                                                                                
            ("Filename for BSE Switch Log logfile (default = " + p_Options->m_LogfileBSESwitchLog + ")").c_str()
        )
        (
            "logfile-BSE-system-parameters",                               
            po::value<string>(&p_Options->m_LogfileBSESystemParameters)->default_value(p_Options->m_LogfileBSESystemParameters),                                                                  
            ("Filename for BSE System Parameters logfile (default = " + p_Options->m_LogfileBSESystemParameters + ")").c_str()
        )
        (
            "logfile-definitions",                                         
            po::value<string>(&p_Options->m_LogfileDefinitionsFilename)->default_value(p_Options->m_LogfileDefinitionsFilename)->implicit_value(""),                                              
            ("Filename for logfile record definitions (default = " + p_Options->m_LogfileDefinitionsFilename + ")").c_str()
        )
        (
            "logfile-delimiter",                                           
            po::value<string>(&p_Options->m_LogfileDelimiterString)->default_value(p_Options->m_LogfileDelimiterString),                                                                          
            ("Field delimiter for logfile records (default = " + p_Options->m_LogfileDelimiterString + ")").c_str()
        )
        (
            "logfile-name-prefix",                                         
            po::value<string>(&p_Options->m_LogfileNamePrefix)->default_value(p_Options->m_LogfileNamePrefix)->implicit_value(""),                                                                
            ("Prefix for logfile names (default = " + p_Options->m_LogfileNamePrefix + ")").c_str()
        )
        (
            "logfile-SSE-parameters",                                      
            po::value<string>(&p_Options->m_LogfileSSEParameters)->default_value(p_Options->m_LogfileSSEParameters),                                                                              
            ("Filename for SSE Parameters logfile (default = " + p_Options->m_LogfileSSEParameters + ")").c_str()
        )
        (
            "logfile-SSE-supernova",                                       
            po::value<string>(&p_Options->m_LogfileSSESupernova)->default_value(p_Options->m_LogfileSSESupernova),                                                                                
            ("Filename for SSE Supernova logfile (default = " + p_Options->m_LogfileSSESupernova + ")").c_str()
        )
        (
            "logfile-SSE-switch-log",                                      
            po::value<string>(&p_Options->m_LogfileSSESwitchLog)->default_value(p_Options->m_LogfileSSESwitchLog),                                                                                
            ("Filename for SSE Switch Log logfile (default = " + p_Options->m_LogfileSSESwitchLog + ")").c_str()
        )
        (
            "output-container,c",                                          
            po::value<string>(&p_Options->m_OutputContainerName)->default_value(p_Options->m_OutputContainerName)->implicit_value(""),                                                            
            ("Container (directory) name for output files (default = " + p_Options->m_OutputContainerName + ")").c_str()
        )
        (
            "outputPath,o",                                                
            po::value<string>(&p_Options->m_OutputPathString)->default_value(p_Options->m_OutputPathString)->implicit_value(""),                                                                  
            ("Directory for output (default = " + p_Options->m_OutputPathString + ")").c_str()
        )


        // vector (list) options - alphabetically

        (
            "debug-classes",                                               
            po::value<vector<string>>(&p_Options->m_DebugClasses)->multitoken()->default_value(p_Options->m_DebugClasses),                                                                        
            ("Debug classes enabled (default = " + defaultDebugClasses + ")").c_str()
        )
        (
            "log-classes",                                                 
            po::value<vector<string>>(&p_Options->m_LogClasses)->multitoken()->default_value(p_Options->m_LogClasses),                                                                            
            ("Logging classes enabled (default = " + defaultLogClasses + ")").c_str()
        )
    
        ;   // end the list of options to be added

    } catch (...) {     // unhandled exception - something wrong with one of the declarations above
        ok = false;     // set status
    }

    return ok;
}


bool Options::ObjectOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription) {

    bool ok = true;                             // status - unless a problem occurs

    // add options

    try {

        p_OptionsDescription->add_options()     // begin the list of options to be added - boost syntactic sugar

        // there is no good way of formatting these - the boost syntz doesn't help that much
        // there is just a boatload of options, so this function (and similar functions) are just going to be long...
    

        // boolean options - alphabetically 

        (
            "allow-rlof-at-birth",                                         
            po::value<bool>(&p_Options->m_AllowRLOFAtBirth)->default_value(p_Options->m_AllowRLOFAtBirth)->implicit_value(true),                                                                  
            ("Allow binaries that have one or both stars in RLOF at birth to evolve (default = " + std::string(p_Options->m_AllowRLOFAtBirth ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "allow-touching-at-birth",                                     
            po::value<bool>(&p_Options->m_AllowTouchingAtBirth)->default_value(p_Options->m_AllowTouchingAtBirth)->implicit_value(true),                                                          
            ("Allow binaries that are touching at birth to evolve (default = " + std::string(p_Options->m_AllowTouchingAtBirth ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "angularMomentumConservationDuringCircularisation",            
            po::value<bool>(&p_Options->m_AngularMomentumConservationDuringCircularisation)->default_value(p_Options->m_AngularMomentumConservationDuringCircularisation)->implicit_value(true),  
            ("Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = " + std::string(p_Options->m_AngularMomentumConservationDuringCircularisation ? "TRUE" : "FALSE") + ")").c_str()
        )

        // Serena
        /* 
        (
            "BeBinaries",                                                  
            po::value<bool>(&p_Options->m_BeBinaries)->default_value(p_Options->m_BeBinaries)->implicit_value(true),                                                                              
            ("Enable Be Binaries study (default = " + std::string(p_Options->m_BeBinaries ? "TRUE" : "FALSE") + ")").c_str()
        )
        */

        (
            "circulariseBinaryDuringMassTransfer",                         
            po::value<bool>(&p_Options->m_CirculariseBinaryDuringMassTransfer)->default_value(p_Options->m_CirculariseBinaryDuringMassTransfer)->implicit_value(true),                            
            ("Circularise binary when it enters a Mass Transfer episode (default = " + std::string(p_Options->m_CirculariseBinaryDuringMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "common-envelope-allow-main-sequence-survive",                 
            po::value<bool>(&p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->default_value(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope)->implicit_value(true),          
            ("Allow main sequence stars to survive common envelope evolution (default = " + std::string(p_Options->m_AllowMainSequenceStarToSurviveCommonEnvelope ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "evolve-pulsars",                                              
            po::value<bool>(&p_Options->m_EvolvePulsars)->default_value(p_Options->m_EvolvePulsars)->implicit_value(true),                                                                        
            ("Evolve pulsars (default = " + std::string(p_Options->m_EvolvePulsars ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "evolve-unbound-systems",                                      
            po::value<bool>(&p_Options->m_EvolveUnboundSystems)->default_value(p_Options->m_EvolveUnboundSystems)->implicit_value(true),                                                          
            ("Continue evolving stars even if the binary is disrupted (default = " + std::string(p_Options->m_EvolveUnboundSystems ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "massTransfer",                                                
            po::value<bool>(&p_Options->m_UseMassTransfer)->default_value(p_Options->m_UseMassTransfer)->implicit_value(true),                                                                    
            ("Enable mass transfer (default = " + std::string(p_Options->m_UseMassTransfer ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pair-instability-supernovae",                                 
            po::value<bool>(&p_Options->m_UsePairInstabilitySupernovae)->default_value(p_Options->m_UsePairInstabilitySupernovae)->implicit_value(true),                                          
            ("Enable pair instability supernovae (PISN) (default = " + std::string(p_Options->m_UsePairInstabilitySupernovae ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "print-bool-as-string",                                        
            po::value<bool>(&p_Options->m_PrintBoolAsString)->default_value(p_Options->m_PrintBoolAsString)->implicit_value(true),                                                                
            ("Print boolean properties as 'TRUE' or 'FALSE' (default = " + std::string(p_Options->m_PrintBoolAsString ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "pulsational-pair-instability",                                
            po::value<bool>(&p_Options->m_UsePulsationalPairInstability)->default_value(p_Options->m_UsePulsationalPairInstability)->implicit_value(true),                                        
            ("Enable mass loss due to pulsational-pair-instability (PPI) (default = " + std::string(p_Options->m_UsePulsationalPairInstability ? "TRUE" : "FALSE") + ")").c_str()
        )
        (
            "revised-energy-formalism-Nandez-Ivanova",                     
            po::value<bool>(&p_Options->m_RevisedEnergyFormalismNandezIvanova)->default_value(p_Options->m_RevisedEnergyFormalismNandezIvanova)->implicit_value(true),                            
            ("Enable revised energy formalism (default = " + std::string(p_Options->m_RevisedEnergyFormalismNandezIvanova ? "TRUE" : "FALSE") + ")").c_str()
        )

        (
            "use-mass-loss",                                               
            po::value<bool>(&p_Options->m_UseMassLoss)->default_value(p_Options->m_UseMassLoss)->implicit_value(true),                                                                            
            ("Enable mass loss (default = " + std::string(p_Options->m_UseMassLoss ? "TRUE" : "FALSE") + ")").c_str()
        )


        // numerical options - alphabetically grouped by type

        // unsigned long

        (
            "random-seed",                                                 
            po::value<unsigned long>(&p_Options->m_RandomSeed)->default_value(p_Options->m_RandomSeed),                                                                                           
            ("Random seed to use (default = " + std::to_string(p_Options->m_RandomSeed) + ")").c_str()
        )


        // int

        (
            "maximum-number-timestep-iterations",                          
            po::value<int>(&p_Options->m_MaxNumberOfTimestepIterations)->default_value(p_Options->m_MaxNumberOfTimestepIterations),                                                               
            ("Maximum number of timesteps to evolve binary before giving up (default = " + std::to_string(p_Options->m_MaxNumberOfTimestepIterations) + ")").c_str()
        )
        (
            "single-star-mass-steps",                                      
            po::value<int>(&p_Options->m_SingleStarMassSteps)->default_value(p_Options->m_SingleStarMassSteps),                                                                                   
            ("Specify the number of mass steps for single star evolution (default = " + std::to_string(p_Options->m_SingleStarMassSteps) + ")").c_str()
        )


        // double

        (
            "common-envelope-alpha",                                       
            po::value<double>(&p_Options->m_CommonEnvelopeAlpha)->default_value(p_Options->m_CommonEnvelopeAlpha),                                                                                
            ("Common Envelope efficiency alpha (default = " + std::to_string(p_Options->m_CommonEnvelopeAlpha) + ")").c_str()
        )
        (
            "common-envelope-alpha-thermal",                               
            po::value<double>(&p_Options->m_CommonEnvelopeAlphaThermal)->default_value(p_Options->m_CommonEnvelopeAlphaThermal),                                                                  
            ("Defined such that lambda = alpha_th * lambda_b + (1.0 - alpha_th) * lambda_g (default = " + std::to_string(p_Options->m_CommonEnvelopeAlphaThermal) + ")").c_str()
        )
        (
            "common-envelope-lambda",                                      
            po::value<double>(&p_Options->m_CommonEnvelopeLambda)->default_value(p_Options->m_CommonEnvelopeLambda),                                                                              
            ("Common Envelope lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeLambda) + ")").c_str()
        )
        (
            "common-envelope-lambda-multiplier",                           
            po::value<double>(&p_Options->m_CommonEnvelopeLambdaMultiplier)->default_value(p_Options->m_CommonEnvelopeLambdaMultiplier),                                                          
            ("Multiply lambda by some constant (default = " + std::to_string(p_Options->m_CommonEnvelopeLambdaMultiplier) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-constant",                     
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionConstant)->default_value(p_Options->m_CommonEnvelopeMassAccretionConstant),                                                
            ("Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass (common-envelope-mass-accretion-prescription CONSTANT). Ignored otherwise (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionConstant) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-max",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMax)->default_value(p_Options->m_CommonEnvelopeMassAccretionMax),                                                          
            ("Maximum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMax) + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-min",                          
            po::value<double>(&p_Options->m_CommonEnvelopeMassAccretionMin)->default_value(p_Options->m_CommonEnvelopeMassAccretionMin),                                                          
            ("Minimum amount of mass accreted by NS/BHs during common envelope evolution in solar masses (default = " + std::to_string(p_Options->m_CommonEnvelopeMassAccretionMin) + ")").c_str()
        )
        (
            "common-envelope-recombination-energy-density",                
            po::value<double>(&p_Options->m_CommonEnvelopeRecombinationEnergyDensity)->default_value(p_Options->m_CommonEnvelopeRecombinationEnergyDensity),                                      
            ("Recombination energy density in erg/g (default = " + std::to_string(p_Options->m_CommonEnvelopeRecombinationEnergyDensity) + ")").c_str()
        )
        (
            "common-envelope-slope-Kruckow",                               
            po::value<double>(&p_Options->m_CommonEnvelopeSlopeKruckow)->default_value(p_Options->m_CommonEnvelopeSlopeKruckow),                                                                  
            ("Common Envelope slope for Kruckow lambda (default = " + std::to_string(p_Options->m_CommonEnvelopeSlopeKruckow) + ")").c_str()
        )

        // AVG
        /*
        (
            "critical-mass-ratio-giant-degenerate-accretor",               
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor),              
            ("Critical mass ratio for MT from a giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioGiantDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-giant-non-degenerate-accretor",           
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor),        
            ("Critical mass ratio for MT from a giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor) + ") Specify both giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-degenerate-accretor",        
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor),  
            ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-giant-non-degenerate-accretor",    
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a helium giant star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor) + ") Specify both helium giant flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-HG-degenerate-accretor",           
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor),        
            ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-HG-non-degenerate-accretor",       
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor),  
            ("Critical mass ratio for MT from a helium HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor) + ") Specify both helium HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-MS-degenerate-accretor",           
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor),        
            ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-helium-MS-non-degenerate-accretor",       
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor),  
            ("Critical mass ratio for MT from a helium MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor) + ") Specify both helium MS flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-HG-degenerate-accretor",                  
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor),                    
            ("Critical mass ratio for MT from a HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHGDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-HG-non-degenerate-accretor",              
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor),              
            ("Critical mass ratio for MT from a HG star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioHGNonDegenerateAccretor) + ") Specify both HG flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-high-mass-degenerate-accretor",        
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor),    
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor) + " Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-high-mass-non-degenerate-accretor",    
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor) + ") Specify both MS high mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-low-mass-degenerate-accretor",         
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor),      
            ("Critical mass ratio for MT from a MS star to a degenerate accretor (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor) + " Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-MS-low-mass-non-degenerate-accretor",     
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a MS star (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor) + ") Specify both MS low mass flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-degenerate-accretor",         
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor),    
            ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        (
            "critical-mass-ratio-white-dwarf-non-degenerate-accretor",     
            po::value<double>(&p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor)->default_value(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor), 
            ("Critical mass ratio for MT from a white dwarf (default = " + std::to_string(p_Options->m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor) + ") Specify both white dwarf flags to use. 0 is always stable, <0 is disabled").c_str()
        )
        */

        (
            "eccentricity-max",                                            
            po::value<double>(&p_Options->m_EccentricityDistributionMax)->default_value(p_Options->m_EccentricityDistributionMax),                                                                
            ("Maximum eccentricity to generate (default = " + std::to_string(p_Options->m_EccentricityDistributionMax) + ")").c_str()
        )
        (
            "eccentricity-min",                                            
            po::value<double>(&p_Options->m_EccentricityDistributionMin)->default_value(p_Options->m_EccentricityDistributionMin),                                                                
            ("Minimum eccentricity to generate (default = " + std::to_string(p_Options->m_EccentricityDistributionMin) + ")").c_str()
        )
        (
            "eddington-accretion-factor",                                  
            po::value<double>(&p_Options->m_EddingtonAccretionFactor)->default_value(p_Options->m_EddingtonAccretionFactor),                                                                      
            ("Multiplication factor for eddington accretion for NS & BH, i.e. >1 is super-eddington and 0. is no accretion (default = " + std::to_string(p_Options->m_EddingtonAccretionFactor) + ")").c_str()
        )

        (
            "fix-dimensionless-kick-magnitude",                            
            po::value<double>(&p_Options->m_FixedUK)->default_value(p_Options->m_FixedUK),                                                                                                        
            ("Fix dimensionless kick magnitude uk to this value (default = " + std::to_string(p_Options->m_FixedUK) + ", -ve values false, +ve values true)").c_str()
        )

        (
            "initial-mass-max",                                            
            po::value<double>(&p_Options->m_InitialMassFunctionMax)->default_value(p_Options->m_InitialMassFunctionMax),                                                                          
            ("Maximum mass (in Msol) to generate using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionMax) + ")").c_str()
        )
        (
            "initial-mass-min",                                            
            po::value<double>(&p_Options->m_InitialMassFunctionMin)->default_value(p_Options->m_InitialMassFunctionMin),                                                                          
            ("Minimum mass (in Msol) to generate using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionMin) + ")").c_str()
        )
        (
            "initial-mass-power",                                          
            po::value<double>(&p_Options->m_InitialMassFunctionPower)->default_value(p_Options->m_InitialMassFunctionPower),                                                                      
            ("Single power law power to generate primary mass using given IMF (default = " + std::to_string(p_Options->m_InitialMassFunctionPower) + ")").c_str()
        )

        (
            "kick-direction-power",                                        
            po::value<double>(&p_Options->m_KickDirectionPower)->default_value(p_Options->m_KickDirectionPower),                                                                                  
            ("Power for power law kick direction distribution (default = " + std::to_string(p_Options->m_KickDirectionPower) + " = isotropic, +ve = polar, -ve = in plane)").c_str()
        )
        (
            "kick-magnitude-max",                                          
            po::value<double>(&p_Options->m_KickMagnitudeDistributionMaximum)->default_value(p_Options->m_KickMagnitudeDistributionMaximum),                                                      
            ("Maximum drawn kick magnitude in km s^-1. Ignored if < 0. Must be > 0 if using kick-magnitude-distribution=FLAT (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionMaximum) + ")").c_str()
        )
        (
            "kick-magnitude",                                          
            po::value<double>(&p_Options->m_KickMagnitude)->default_value(p_Options->m_KickMagnitude),                                                      
            ("The magnitude of the kick velocity the star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-1",                                          
            po::value<double>(&p_Options->m_KickMagnitude1)->default_value(p_Options->m_KickMagnitude1),                                                      
            ("The magnitude of the kick velocity the primary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude1) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-2",                                          
            po::value<double>(&p_Options->m_KickMagnitude2)->default_value(p_Options->m_KickMagnitude2),                                                      
            ("The magnitude of the kick velocity the secondary star receives during the a supernova (default = " + std::to_string(p_Options->m_KickMagnitude2) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-random",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom)->default_value(p_Options->m_KickMagnitudeRandom),                                                      
            ("Number used to choose the kick velocity magnitude for the star during the a supernova (default = uniform random number [0.0, 1.0))").c_str()
        )
        (
            "kick-magnitude-random-1",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom1)->default_value(p_Options->m_KickMagnitudeRandom1),                                                      
            ("Number used to choose the kick velocity magnitude for the primary star during the a supernova (default = uniform random number [0.0, 1.0))").c_str()
        )
        (
            "kick-magnitude-random-2",                                          
            po::value<double>(&p_Options->m_KickMagnitudeRandom2)->default_value(p_Options->m_KickMagnitudeRandom2),                                                      
            ("Number used to choose the kick velocity magnitude for the secondary during the a supernova (default = uniform random number [0.0, 1.0))").c_str()
        )
        (
            "kick-magnitude-sigma-CCSN-BH",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH),                                            
            ("Sigma for chosen kick magnitude distribution for black holes (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_BH) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-CCSN-NS",                                
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS)->default_value(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS),                                            
            ("Sigma for chosen kick magnitude distribution for neutron stars (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaCCSN_NS) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-ECSN",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForECSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForECSN),                                            
            ("Sigma for chosen kick magnitude distribution for ECSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForECSN) + " km s^-1 )").c_str()
        )
        (
            "kick-magnitude-sigma-USSN",                                   
            po::value<double>(&p_Options->m_KickMagnitudeDistributionSigmaForUSSN)->default_value(p_Options->m_KickMagnitudeDistributionSigmaForUSSN),                                            
            ("Sigma for chosen kick magnitude distribution for USSN (default = " + std::to_string(p_Options->m_KickMagnitudeDistributionSigmaForUSSN) + " km s^-1 )").c_str()
        )
        (
            "kick-mean-anomaly-1",
            po::value<double>(&p_Options->m_KickMeanAnomaly1)->default_value(p_Options->m_SNmeanAnomaly1),                                                                                  
            ("Mean anomaly for the primary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))").c_str()
        )
        (
            "kick-mean-anomaly-2",
            po::value<double>(&p_Options->m_KickMeanAnomaly2)->default_value(p_Options->m_SNmeanAnomaly2),                                                                                  
            ("Mean anomaly for the secondary star at instantaneous time of the supernova (default = uniform random number [0.0, 2pi))").c_str()
        )
        (
            "kick-phi-1",
            po::value<double>(&p_Options->m_KickPhi1)->default_value(p_Options->m_SNphi1),                                                                                  
            ("Angle between 'x' and 'y', both in the orbital plane of the supernovae vector, for the primary star (default = drawn from kick direction distribution)").c_str()
        )
        (
            "kick-phi-2",
            po::value<double>(&p_Options->m_KickPhi2)->default_value(p_Options->m_SNphi2),                                                                                  
            ("Angle between 'x' and 'y', both in the orbital plane of the supernovae vector, for the secondary star (default = drawn from kick direction distribution)").c_str()
        )
        (
            "kick-scaling-factor",                                         
            po::value<double>(&p_Options->m_KickScalingFactor)->default_value(p_Options->m_KickScalingFactor),                                                                                    
            ("Arbitrary factor used to scale kicks (default = " + std::to_string(p_Options->m_KickScalingFactor) + ")").c_str()
        )
        (
            "kick-theta-1",                                        
            po::value<double>(&p_Options->m_KickTheta1)->default_value(p_Options->m_KickTheta1),                                                                                  
            ("Angle between the orbital plane and the 'z' axis of the supernovae vector, for the primary star (default = drawn from kick direction distribution)").c_str()
        )
        (
            "kick-theta-2",                                        
            po::value<double>(&p_Options->m_KickTheta2)->default_value(p_Options->m_KickTheta2),                                                                                  
            ("Angle between the orbital plane and the 'z' axis of the supernovae vector, for the secondary star (default = drawn from kick direction distribution)").c_str()
        )

        (
            "luminous-blue-variable-multiplier",                           
            po::value<double>(&p_Options->m_LuminousBlueVariableFactor)->default_value(p_Options->m_LuminousBlueVariableFactor),                                                                  
            ("Multiplicitive constant for LBV mass loss (default = " + std::to_string(p_Options->m_LuminousBlueVariableFactor) + ", use 10 for Mennekens & Vanbeveren 2014)").c_str()
        )

        (
            "mass-ratio-max",                                              
            po::value<double>(&p_Options->m_MassRatioDistributionMax)->default_value(p_Options->m_MassRatioDistributionMax),                                                                      
            ("Maximum mass ratio m2/m1 to generate (default = " + std::to_string(p_Options->m_MassRatioDistributionMax) + ")").c_str()
        )
        (
            "mass-ratio-min",                                              
            po::value<double>(&p_Options->m_MassRatioDistributionMin)->default_value(p_Options->m_MassRatioDistributionMin),                                                                      
            ("Minimum mass ratio m2/m1 to generate (default = " + std::to_string(p_Options->m_MassRatioDistributionMin) + ")").c_str()
        )
        (
            "mass-transfer-fa",                                            
            po::value<double>(&p_Options->m_MassTransferFractionAccreted)->default_value(p_Options->m_MassTransferFractionAccreted),                                                              
            ("Mass Transfer fraction accreted in FIXED prescription (default = " + std::to_string(p_Options->m_MassTransferFractionAccreted) + ", fully conservative)").c_str()
        )
        (
            "mass-transfer-jloss",                                         
            po::value<double>(&p_Options->m_MassTransferJloss)->default_value(p_Options->m_MassTransferJloss),                                                                                    
            ("Specific angular momentum with which the non-accreted system leaves the system (default = " + std::to_string(p_Options->m_MassTransferJloss) + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-C",                               
            po::value<double>(&p_Options->m_MassTransferCParameter)->default_value(p_Options->m_MassTransferCParameter),                                                                          
            ("Mass Transfer Thermal rate factor fo the accretor (default = " + std::to_string(p_Options->m_MassTransferCParameter) + ")").c_str()
        )
        (
            "maximum-evolution-time",                                      
            po::value<double>(&p_Options->m_MaxEvolutionTime)->default_value(p_Options->m_MaxEvolutionTime),                                                                                      
            ("Maximum time to evolve binaries in Myrs (default = " + std::to_string(p_Options->m_MaxEvolutionTime) + ")").c_str()
        )
        (
            "maximum-mass-donor-Nandez-Ivanova",                           
            po::value<double>(&p_Options->m_MaximumMassDonorNandezIvanova)->default_value(p_Options->m_MaximumMassDonorNandezIvanova),                                                            
            ("Maximum donor mass allowed for the revised common envelope formalism in Msol (default = " + std::to_string(p_Options->m_MaximumMassDonorNandezIvanova) + ")").c_str()
        )
        (
            "maximum-neutron-star-mass",                                   
            po::value<double>(&p_Options->m_MaximumNeutronStarMass)->default_value(p_Options->m_MaximumNeutronStarMass),                                                                          
            ("Maximum mass of a neutron star (default = " + std::to_string(p_Options->m_MaximumNeutronStarMass) + ")").c_str()
        )
        (
            "MCBUR1",                                                      
            po::value<double>(&p_Options->m_mCBUR1)->default_value(p_Options->m_mCBUR1),                                                                                                          
            ("MCBUR1: Min core mass at BAGB to avoid fully degenerate CO core  (default = " + std::to_string(p_Options->m_mCBUR1) + ")").c_str()
        )
        (
            "metallicity,z",                                               
            po::value<double>(&p_Options->m_Metallicity)->default_value(p_Options->m_Metallicity),                                                                                                
            ("Metallicity to use (default " + std::to_string(p_Options->m_Metallicity) + " Zsol)").c_str()
        )
        (
            "minimum-secondary-mass",                                      
            po::value<double>(&p_Options->m_MinimumMassSecondary)->default_value(p_Options->m_MinimumMassSecondary),                                                                              
            ("Minimum mass of secondary to generate in Msol (default = " + std::to_string(p_Options->m_MinimumMassSecondary) + ")").c_str()
        )

        (
            "neutrino-mass-loss-bh-formation-value",                       
            po::value<double>(&p_Options->m_NeutrinoMassLossValueBH)->default_value(p_Options->m_NeutrinoMassLossValueBH),                                                                        
            ("Value corresponding to neutrino mass loss assumption (default = " + std::to_string(p_Options->m_NeutrinoMassLossValueBH) + ")").c_str()
        )

        (
            "orbital-period-max",                                          
            po::value<double>(&p_Options->m_PeriodDistributionMax)->default_value(p_Options->m_PeriodDistributionMax),                                                                            
            ("Maximum period in days to generate (default = " + std::to_string(p_Options->m_PeriodDistributionMax) + ")").c_str()
        )
        (
            "orbital-period-min",                                          
            po::value<double>(&p_Options->m_PeriodDistributionMin)->default_value(p_Options->m_PeriodDistributionMin),                                                                            
            ("Minimum period in days to generate (default = " + std::to_string(p_Options->m_PeriodDistributionMin) + ")").c_str()
        )

        (
            "PISN-lower-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityLowerLimit)->default_value(p_Options->m_PairInstabilityLowerLimit),                                                                    
            ("Minimum core mass for PISN (default = " + std::to_string(p_Options->m_PairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "PISN-upper-limit",                                            
            po::value<double>(&p_Options->m_PairInstabilityUpperLimit)->default_value(p_Options->m_PairInstabilityUpperLimit),                                                                    
            ("Maximum core mass for PISN (default = " + std::to_string(p_Options->m_PairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "PPI-lower-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityLowerLimit)->default_value(p_Options->m_PulsationalPairInstabilityLowerLimit),                                              
            ("Minimum core mass for PPI (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityLowerLimit) + ")").c_str()
        )
        (
            "PPI-upper-limit",                                             
            po::value<double>(&p_Options->m_PulsationalPairInstabilityUpperLimit)->default_value(p_Options->m_PulsationalPairInstabilityUpperLimit),                                              
            ("Maximum core mass for PPI (default = " + std::to_string(p_Options->m_PulsationalPairInstabilityUpperLimit) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-max",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMax)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMax),                                        
            ("Maximum (log10) pulsar birth magnetic field (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-magnetic-field-distribution-min",                
            po::value<double>(&p_Options->m_PulsarBirthMagneticFieldDistributionMin)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionMin),                                        
            ("Minimum (log10) pulsar birth magnetic field) (default = " + std::to_string(p_Options->m_PulsarBirthMagneticFieldDistributionMin) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-max",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMax)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMax),                                              
            ("Maximum pulsar birth spin period in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMax) + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution-min",                   
            po::value<double>(&p_Options->m_PulsarBirthSpinPeriodDistributionMin)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionMin),                                              
            ("Minimum pulsar birth spin period in ms (default = " + std::to_string(p_Options->m_PulsarBirthSpinPeriodDistributionMin) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-massscale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayMassscale)->default_value(p_Options->m_PulsarMagneticFieldDecayMassscale),                                                    
            ("Mass scale on which magnetic field decays during accretion in solar masses (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayMassscale) + ")").c_str()
        )
        (
            "pulsar-magnetic-field-decay-timescale",                       
            po::value<double>(&p_Options->m_PulsarMagneticFieldDecayTimescale)->default_value(p_Options->m_PulsarMagneticFieldDecayTimescale),                                                    
            ("Timescale on which magnetic field decays in Myrs (default = " + std::to_string(p_Options->m_PulsarMagneticFieldDecayTimescale) + ")").c_str()
        )
        (
            "pulsar-minimum-magnetic-field",                               
            po::value<double>(&p_Options->m_PulsarLog10MinimumMagneticField)->default_value(p_Options->m_PulsarLog10MinimumMagneticField),                                                        
            ("log10 of the minimum pulsar magnetic field in Gauss (default = " + std::to_string(p_Options->m_PulsarLog10MinimumMagneticField) + ")").c_str()
        )

        (
            "semi-major-axis-max",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMax)->default_value(p_Options->m_SemiMajorAxisDistributionMax),                                                              
            ("Maximum semi major axis in AU to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMax) + ")").c_str()
        )
        (
            "semi-major-axis-min",                                         
            po::value<double>(&p_Options->m_SemiMajorAxisDistributionMin)->default_value(p_Options->m_SemiMajorAxisDistributionMin),                                                              
            ("Minimum semi major axis in AU to generate (default = " + std::to_string(p_Options->m_SemiMajorAxisDistributionMin) + ")").c_str()
        )
        (
            "single-star-mass-max",                                        
            po::value<double>(&p_Options->m_SingleStarMassMax)->default_value(p_Options->m_SingleStarMassMax),                                                                                    
            ("Maximum mass (in Msol) for single star evolution (default = " + std::to_string(p_Options->m_SingleStarMassMax) + ")").c_str()
        )
        (
            "single-star-mass-min",                                        
            po::value<double>(&p_Options->m_SingleStarMassMin)->default_value(p_Options->m_SingleStarMassMin),                                                                                    
            ("Minimum mass (in Msol) for single star evolution (default = " + std::to_string(p_Options->m_SingleStarMassMin) + ")").c_str()
        )

        (
            "wolf-rayet-multiplier",                                       
            po::value<double>(&p_Options->m_WolfRayetFactor)->default_value(p_Options->m_WolfRayetFactor),                                                                                        
            ("Multiplicitive constant for WR winds (default = " + std::to_string(p_Options->m_WolfRayetFactor) + ")").c_str()
        )

        (
            "zeta-adiabatic-arbitrary",                                    
            po::value<double>(&p_Options->m_ZetaAdiabaticArbitrary)->default_value(p_Options->m_ZetaAdiabaticArbitrary),                                                                          
            ("Value of mass-radius exponent zeta adiabatic (default = " + std::to_string(p_Options->m_ZetaAdiabaticArbitrary) + ")").c_str()
        )
        (
            "zeta-main-sequence",                                          
            po::value<double>(&p_Options->m_ZetaMainSequence)->default_value(p_Options->m_ZetaMainSequence),                                                                                      
            ("Value of mass-radius exponent zeta on the main sequence (default = " + std::to_string(p_Options->m_ZetaMainSequence) + ")").c_str()
        )
        (
            "zeta-radiative-envelope-giant",                               
            po::value<double>(&p_Options->m_ZetaRadiativeEnvelopeGiant)->default_value(p_Options->m_ZetaRadiativeEnvelopeGiant),                                                                  
            ("Value of mass-radius exponent zeta for radiative envelope giants (default = " + std::to_string(p_Options->m_ZetaRadiativeEnvelopeGiant) + ")").c_str()
        )


        // string options - alphabetically

        (
            "black-hole-kicks",                                            
            po::value<string>(&p_Options->m_BlackHoleKicksOptionString)->default_value(p_Options->m_BlackHoleKicksOptionString),                                                                              
            ("Black hole kicks relative to NS kicks (options: FULL, REDUCED, ZERO, FALLBACK), default = " + p_Options->m_BlackHoleKicksOptionString + ")").c_str()
        )

        (
            "case-bb-stability-prescription",                              
            po::value<string>(&p_Options->m_CaseBBStabilityPrescriptionString)->default_value(p_Options->m_CaseBBStabilityPrescriptionString),                                                    
            ("Case BB/BC mass transfer stability prescription (options: ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE), default = " + p_Options->m_CaseBBStabilityPrescriptionString + ")").c_str()
        )
        (
            "chemically-homogeneous-evolution",                            
            po::value<string>(&p_Options->m_CheString)->default_value(p_Options->m_CheString),                                                                                                    
            ("Chemically Homogeneous Evolution (options: NONE, OPTIMISTIC, PESSIMISTIC), default = " + p_Options->m_CheString + ")").c_str()
        )
        (
            "common-envelope-lambda-prescription",                         
            po::value<string>(&p_Options->m_CommonEnvelopeLambdaPrescriptionString)->default_value(p_Options->m_CommonEnvelopeLambdaPrescriptionString),                                          
            ("CE lambda prescription (options: LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI), default = " + p_Options->m_CommonEnvelopeLambdaPrescriptionString + ")").c_str()
        )
        (
            "common-envelope-mass-accretion-prescription",                 
            po::value<string>(&p_Options->m_CommonEnvelopeMassAccretionPrescriptionString)->default_value(p_Options->m_CommonEnvelopeMassAccretionPrescriptionString),                            
            ("Assumption about whether NS/BHs can accrete mass during common envelope evolution (options: ZERO, CONSTANT, UNIFORM, MACLEOD), default = " + p_Options->m_CommonEnvelopeMassAccretionPrescriptionString + ")").c_str()
        )
        
        (
            "eccentricity-distribution,e",                                 
            po::value<string>(&p_Options->m_EccentricityDistributionString)->default_value(p_Options->m_EccentricityDistributionString),                                                          
            ("Initial eccentricity distribution, e (options: ZERO, FIXED, FLAT, THERMALISED, GELLER+2013), default = " + p_Options->m_EccentricityDistributionString + ")").c_str()
        )
        (
            "envelope-state-prescription",                                 
            po::value<string>(&p_Options->m_EnvelopeStatePrescriptionString)->default_value(p_Options->m_EnvelopeStatePrescriptionString),                                                        
            ("Prescription for whether the envelope is radiative or convective (options: LEGACY, HURLEY, FIXED_TEMPERATURE), default = " + p_Options->m_EnvelopeStatePrescriptionString + ")").c_str()
        )

        (
            "fryer-supernova-engine",                                      
            po::value<string>(&p_Options->m_FryerSupernovaEngineString)->default_value(p_Options->m_FryerSupernovaEngineString),                                                                  
            ("If using Fryer et al 2012 fallback prescription, select between 'delayed' and 'rapid' engines (default = " + p_Options->m_FryerSupernovaEngineString + ")").c_str()
        )

        (
            "initial-mass-function,i",                                     
            po::value<string>(&p_Options->m_InitialMassFunctionString)->default_value(p_Options->m_InitialMassFunctionString),                                                                    
            ("Initial mass function (options: SALPETER, POWERLAW, UNIFORM, KROUPA), default = " + p_Options->m_InitialMassFunctionString + ")").c_str()
        )

        (
            "kick-direction",                                              
            po::value<string>(&p_Options->m_KickDirectionDistributionString)->default_value(p_Options->m_KickDirectionDistributionString),                                                        
            ("Natal kick direction distribution (options: ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES), default = " + p_Options->m_KickDirectionDistributionString + ")").c_str()
        )
        (
            "kick-magnitude-distribution",                                 
            po::value<string>(&p_Options->m_KickMagnitudeDistributionString)->default_value(p_Options->m_KickMagnitudeDistributionString),                                                        
            ("Natal kick magnitude distribution (options: ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL), default = " + p_Options->m_KickMagnitudeDistributionString + ")").c_str()
        )

        (
            "mass-loss-prescription",                                      
            po::value<string>(&p_Options->m_MassLossPrescriptionString)->default_value(p_Options->m_MassLossPrescriptionString),                                                                  
            ("Mass loss prescription (options: NONE, HURLEY, VINK), default = " + p_Options->m_MassLossPrescriptionString + ")").c_str()
        )
        (
            "mass-ratio-distribution,q",                                   
            po::value<string>(&p_Options->m_MassRatioDistributionString)->default_value(p_Options->m_MassRatioDistributionString),                                                                
            ("Initial mass ratio distribution for q=m2/m1 (options: FLAT, DuquennoyMayor1991, SANA2012), default = " + p_Options->m_MassRatioDistributionString + ")").c_str()
        )
        (
            "mass-transfer-accretion-efficiency-prescription",             
            po::value<string>(&p_Options->m_MassTransferAccretionEfficiencyPrescriptionString)->default_value(p_Options->m_MassTransferAccretionEfficiencyPrescriptionString),                    
            ("Mass Transfer Accretion Efficiency prescription (options: THERMAL, FIXED), default = " + p_Options->m_MassTransferAccretionEfficiencyPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-angular-momentum-loss-prescription",            
            po::value<string>(&p_Options->m_MassTransferAngularMomentumLossPrescriptionString)->default_value(p_Options->m_MassTransferAngularMomentumLossPrescriptionString),                    
            ("Mass Transfer Angular Momentum Loss prescription (options: JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY), default = " + p_Options->m_MassTransferAngularMomentumLossPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-rejuvenation-prescription",                     
            po::value<string>(&p_Options->m_MassTransferRejuvenationPrescriptionString)->default_value(p_Options->m_MassTransferRejuvenationPrescriptionString),                                  
            ("Mass Transfer Rejuvenation prescription (options: NONE, STARTRACK), default = " + p_Options->m_MassTransferRejuvenationPrescriptionString + ")").c_str()
        )
        (
            "mass-transfer-thermal-limit-accretor",                        
            po::value<string>(&p_Options->m_MassTransferThermallyLimitedVariationString)->default_value(p_Options->m_MassTransferThermallyLimitedVariationString),                                
            ("Mass Transfer Thermal Accretion limit (default = " + p_Options->m_MassTransferThermallyLimitedVariationString + ")").c_str()
        )

        (
            "neutrino-mass-loss-bh-formation",                             
            po::value<string>(&p_Options->m_NeutrinoMassLossAssumptionBHString)->default_value(p_Options->m_NeutrinoMassLossAssumptionBHString),                                                  
            ("Assumption about neutrino mass loss during BH formation (options: FIXED_FRACTION, FIXED_MASS), default = " + p_Options->m_NeutrinoMassLossAssumptionBHString + ")").c_str()
        )
        (
            "neutron-star-equation-of-state",                              
            po::value<string>(&p_Options->m_NeutronStarEquationOfStateString)->default_value(p_Options->m_NeutronStarEquationOfStateString),                                                      
            ("Neutron star equation of state to use (options: SSE, ARP3), default = " + p_Options->m_NeutronStarEquationOfStateString + ")").c_str()
        )

        (
            "pulsar-birth-magnetic-field-distribution",                    
            po::value<string>(&p_Options->m_PulsarBirthMagneticFieldDistributionString)->default_value(p_Options->m_PulsarBirthMagneticFieldDistributionString),                                  
            ("Pulsar Birth Magnetic Field distribution (options: ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL), default = " + p_Options->m_PulsarBirthMagneticFieldDistributionString + ")").c_str()
        )
        (
            "pulsar-birth-spin-period-distribution",                       
            po::value<string>(&p_Options->m_PulsarBirthSpinPeriodDistributionString)->default_value(p_Options->m_PulsarBirthSpinPeriodDistributionString),                                        
            ("Pulsar Birth Spin Period distribution (options: ZERO, FIXED, UNIFORM, NORMAL), default = " + p_Options->m_PulsarBirthSpinPeriodDistributionString + ")").c_str()
        )
        (
            "pulsational-pair-instability-prescription",                   
            po::value<string>(&p_Options->m_PulsationalPairInstabilityPrescriptionString)->default_value(p_Options->m_PulsationalPairInstabilityPrescriptionString),                              
            ("Pulsational Pair Instability prescription (options: COMPAS, STARTRACK, MARCHANT), default = " + p_Options->m_PulsationalPairInstabilityPrescriptionString + ")").c_str()
        )

        (
            "remnant-mass-prescription",                                   
            po::value<string>(&p_Options->m_RemnantMassPrescriptionString)->default_value(p_Options->m_RemnantMassPrescriptionString),                                                            
            ("Choose remnant mass prescription (options: HURLEY2000, BELCZYNSKI2002, FRYER2012, MULLER2016, MULLERMANDEL), default = " + p_Options->m_RemnantMassPrescriptionString + ")").c_str()
        )
        (
            "rotational-velocity-distribution",                            
            po::value<string>(&p_Options->m_RotationalVelocityDistributionString)->default_value(p_Options->m_RotationalVelocityDistributionString),                                              
            ("Initial rotational velocity distribution (options: ZERO, HURLEY, VLTFLAMES), default = " + p_Options->m_RotationalVelocityDistributionString + ")").c_str()
        )

        (
            "semi-major-axis-distribution,a",                              
            po::value<string>(&p_Options->m_SemiMajorAxisDistributionString)->default_value(p_Options->m_SemiMajorAxisDistributionString),                                                        
            ("Initial semi-major axis distribution, a (options: FLATINLOG, CUSTOM, DuquennoyMayor1991, SANA2012), default = " + p_Options->m_SemiMajorAxisDistributionString + ")").c_str()
        )        
        (
            "stellar-zeta-prescription",                                   
            po::value<string>(&p_Options->m_StellarZetaPrescriptionString)->default_value(p_Options->m_StellarZetaPrescriptionString),                                                            
            ("Prescription for stellar zeta (default = " + p_Options->m_StellarZetaPrescriptionString + ")").c_str()
        )
   
        ;   // end the list of options to be added

    } catch (...) {     // unhandled exception - something wrong with one of the declarations above
        ok = false;     // set status
    }

    return ok;
}


string Options::ProgramOptionDetails(const OptionValues *p_Options, const po::variables_map p_VM) {

    TYPENAME tipe = TYPENAME::NONE;      // for grid file processing

    std::ostringstream ss;

    ss << "COMMAND LINE OPTIONS\n-------------------\n\n";

    for (po::variables_map::const_iterator it = p_VM.begin(); it != p_VM.end(); it++) {

        // option name
        ss << it->first << " = ";

        if (((boost::any)it->second.value()).empty()) ss << "<EMPTY_OPTION>\n";
        else {

            // determine if option values was supplied, or whether the default was used

            string valueSource;
            if (p_VM[it->first].defaulted() || it->second.defaulted()) valueSource = "DEFAULT_USED";
            else                                                       valueSource = "USER_SUPPLIED";

            // find data type & print value
            // handles most data types - add others if they cause problems

            bool isCharPtr = false;
            bool isStr     = false;

            // (pre)check for data type = charPtr
            try {
                boost::any_cast<const char *>(it->second.value());
                isCharPtr = true;
            } catch (const boost::bad_any_cast &) {
                isCharPtr = false;
            }

            if (!isCharPtr) {
                // (pre)check for data type = string
                try {
                    boost::any_cast<std::string>(it->second.value());
                    isStr = true;
                } catch (const boost::bad_any_cast &) {
                    isStr = false;
                }
            }

            // find other data types
            // it's not pretty, but it works

            if (isCharPtr) { 
                tipe = TYPENAME::STRING; 
                ss << p_VM[it->first].as<const char *>() << ", " << valueSource << ", CONST_CHAR_*";
            }

            else if (isStr) {
                tipe = TYPENAME::STRING; 
                std::string tmp = p_VM[it->first].as<std::string>();
                if (tmp.size()) ss << "'" << tmp << "'";
                else            ss << "''";
                ss << ", " << valueSource << ", STRING";
            }

            else if (((boost::any)it->second.value()).type() == typeid(signed                )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed                >() << ", " << valueSource << ", SIGNED"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned              )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned              >() << ", " << valueSource << ", UNSIGNED"; }

            else if (((boost::any)it->second.value()).type() == typeid(short                 )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<short                 >() << ", " << valueSource << ", SHORT"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed short          )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed short          >() << ", " << valueSource << ", SIGNED_SHORT"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned short        )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned short        >() << ", " << valueSource << ", UNSIGNED_SHORT"; }

            else if (((boost::any)it->second.value()).type() == typeid(short int             )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<short int             >() << ", " << valueSource << ", SHORT_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed short int      )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed short int      >() << ", " << valueSource << ", SIGNED_SHORT_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned short int    )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned short int    >() << ", " << valueSource << ", UNSIGNED_SHORT_INT"; }

            else if (((boost::any)it->second.value()).type() == typeid(int                   )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<int                   >() << ", " << valueSource << ", INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed int            )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed int            >() << ", " << valueSource << ", SIGNED_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned int          )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned int          >() << ", " << valueSource << ", UNSIGNED_INT"; }

            else if (((boost::any)it->second.value()).type() == typeid(long                  )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<long                  >() << ", " << valueSource << ", LONG"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed long           )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed long           >() << ", " << valueSource << ", SIGNED_LONG"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long         )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned long         >() << ", " << valueSource << ", UNSIGNED_LONG"; }

            else if (((boost::any)it->second.value()).type() == typeid(long int              )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<long int              >() << ", " << valueSource << ", LONG_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed long int       )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed long int       >() << ", " << valueSource << ", SIGNED_LONG_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long int     )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned long int     >() << ", " << valueSource << ", UNSIGNED_LONG_INT"; }

            else if (((boost::any)it->second.value()).type() == typeid(long long             )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<long long             >() << ", " << valueSource << ", LONG_LONG"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed long long      )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed long long      >() << ", " << valueSource << ", SIGNED_LONG_LONG"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long long    )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned long long    >() << ", " << valueSource << ", UNSIGNED_LONG_LONG"; }

            else if (((boost::any)it->second.value()).type() == typeid(long long int         )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<long long int         >() << ", " << valueSource << ", LONG_LONG_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(signed long long int  )) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<signed long long int  >() << ", " << valueSource << ", SIGNED_LONG_LONG_INT"; }
            else if (((boost::any)it->second.value()).type() == typeid(unsigned long long int)) { tipe = TYPENAME::INT;   ss << p_VM[it->first].as<unsigned long long int>() << ", " << valueSource << ", UNSIGNED_LONG_LONG_INT"; }

            else if (((boost::any)it->second.value()).type() == typeid(float                 )) { tipe = TYPENAME::FLOAT; ss << p_VM[it->first].as<float                 >() << ", " << valueSource << ", FLOAT"; }
            else if (((boost::any)it->second.value()).type() == typeid(double                )) { tipe = TYPENAME::FLOAT; ss << p_VM[it->first].as<double                >() << ", " << valueSource << ", DOUBLE"; }
            else if (((boost::any)it->second.value()).type() == typeid(long double           )) { tipe = TYPENAME::FLOAT; ss << p_VM[it->first].as<long double           >() << ", " << valueSource << ", LONG_DOUBLE"; }

            else if (((boost::any)it->second.value()).type() == typeid(char                  )) ss << p_VM[it->first].as<char                  >() << ", " << valueSource << ", CHAR";
            else if (((boost::any)it->second.value()).type() == typeid(signed char           )) ss << p_VM[it->first].as<signed char           >() << ", " << valueSource << ", SIGNED_CHAR";
            else if (((boost::any)it->second.value()).type() == typeid(unsigned char         )) ss << p_VM[it->first].as<unsigned char         >() << ", " << valueSource << ", UNSIGNED_CHAR";

            else if (((boost::any)it->second.value()).type() == typeid(bool)) {
                tipe = TYPENAME::BOOL;
                bool v = p_VM[it->first].as<bool>();
                ss << (v ? "TRUE" : "FALSE") << ", " << valueSource << ", BOOL";
            } 

            else {  // Assume vector<string>
                try {
                    std::ostringstream elemsSS;
                    elemsSS << "{ ";
                    vector<string> tmp = p_VM[it->first].as<vector<string>>();
                    for (std::vector<string>::iterator elem=tmp.begin(); elem != tmp.end(); elem++) {
                        elemsSS << "'" << (*elem) << "', ";
                    }
                    string elems = elemsSS.str();
                    if (elems.size() > 2) elems.erase(elems.size() - 2);
                    else if (elems.size() == 2) elems.erase(elems.size() - 1);
                    elems += " }";
                    ss << elems << ", " << valueSource << ", VECTOR<STRING>";
                } catch (const boost::bad_any_cast &) {
                    ss << "<UNKNOWN_DATA_TYPE>, " << valueSource << ", <UNKNOWN_DATA_TYPE>";
                }
            }

//            // get short option name
//            po::option_description const& opt = p_OptionDescriptions.find(it->first, false, false, false);
//            std::string shortOptionName = opt.canonical_display_name(cls::allow_dash_for_short);
//            if ((shortOptionName != it->first) && (shortOptionName[0] == '-')) shortOptionName.erase(0, 1);
                 
//            // add option to option map (for grid file processing)
//            m_OptionMap[it->first] = std::make_tuple(tipe, shortOptionName);
//            if (shortOptionName != it->first) m_OptionMap[shortOptionName] = std::make_tuple(tipe, it->first);
        }

        ss << "\n";
    }
  
    ss << "\n\nOTHER PARAMETERS\n----------------\n\n";

    ss << "fixedMetallicity   = " << (p_Options->m_FixedMetallicity ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                           // fixedMetallicity
    ss << "useFixedUK         = " << (p_Options->m_UseFixedUK ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                                 // useFixedUK
    ss << "outputPath         = " << p_Options->m_OutputPath.string() << ", CALCULATED, STRING\n";                                           // outputPath (fully qualified)
    ss << "fixedRandomSeed    = " << (p_Options->m_FixedRandomSeed ? "TRUE" : "FALSE") << ", CALCULATED, BOOL\n";                            // fixedRandomSeed




//    std::cout << "Using Boost "     
//              << BOOST_VERSION / 100000     << "."  // major version
//              << BOOST_VERSION / 100 % 1000 << "."  // minor version
//              << BOOST_VERSION % 100                // patch level
//              << std::endl;


    // GSL library version?????????????????????????????


    return ss.str();
}


    
            // sanity check options and option values
            //
            // The boost library provides a mechanism to define a callback for each of the
            // options - boost calls them 'notifier functions'.  If any options have a notifier
            // function specified, calling po::notify(vm) here would cause each of the notifier
            // functions to be called - the idea is that rather than having a block of code
            // that acts on each option in turn, we can instead write small pieces of code for
            // each option and encapsulate the code for each option in its own function.  We  YES WE DO !!! HERE WE SET THE VALUE!!!!!!!!!!
            // don't act on options that way, but we could use the functionality to do the sanity
            // checks here, rather than have a big block of code that sanity checks each otion   WE WANT TO SET TWICE SO CAN'T USE - ONE FOR GLOBAL , ONE FOR OBJECT
            // in turn.  However, we have so many options all we'd be doing by using that
            // functionality is replacing a monolithic block of code that checks many options
            // with many functions - I'm not sure which is worse/better... Using notifier functions
            // is more in keeping with the OO pattern, and it would encapsulate knowledge of each
            // option more neatly, but haveing so many notifier functions is a bit overwhelming.
            //
            // So I'm leaving this as a big block of code that checks each option in turn, but 
            // wouldn't be averse to eventually implementing notifier functions.
            //
            // Note that if we don't implement notifier functions we don't need to call po::notify(vm)
            // If we do decide to implement notifier functions, we could put the code above that
            // acts on --help and --version into notifier functions - probably not necessaryfor now.

std::string Options::OptionValues::CheckAndSetOptions() {

    // MASS, or MASS1 + MASS2 ARE REQUIRED OPTIONS IN A GRID FILE!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    m_FixedRandomSeed  = !m_VM["random-seed"].defaulted();                                              // use random seed if it is provided by the user
    m_FixedMetallicity = !m_VM["metallicity"].defaulted();                                              // determine if user supplied a metallicity value
    m_UseFixedUK       = !m_VM["fix-dimensionless-kick-magnitude"].defaulted() && (m_FixedUK >= 0.0);   // determine if user supplied a valid kick magnitude


    // check & set prescriptions, distributions, assumptions etc. options - alphabetically

    bool found;

    // Floor
    /*
    if (!vmCmdLine["AIS-DCOtype"].defaulted()) {                                                        // Adaptive Importance Sampling DCO type
        std::tie(found, p_OptionValues->m_AISDCOtype) = utils::GetMapKey(p_OptionValues->m_AISDCOtypeString, AIS_DCO_LABEL, p_OptionValues->m_AISDCOtype);
        return "Unknown AIS DCO Type";
    }
    */

    if (!m_VM["black-hole-kicks"].defaulted()) {                                                        // black hole kicks option
        std::tie(found, m_BlackHoleKicksOption) = utils::GetMapKey(m_BlackHoleKicksOptionString, BLACK_HOLE_KICK_OPTION_LABEL, m_BlackHoleKicksOption);
        if (!found) return "Unknown Black Hole Kicks Option";
    }

    if (!m_VM["case-bb-stability-prescription"].defaulted()) {                                          //case BB/BC mass transfer stability prescription
        std::tie(found, m_CaseBBStabilityPrescription) = utils::GetMapKey(m_CaseBBStabilityPrescriptionString, CASE_BB_STABILITY_PRESCRIPTION_LABEL, m_CaseBBStabilityPrescription);
        if (!found) return "Unknown Case BB/BC Mass Transfer Stability Prescription";
    }
           
    if (!m_VM["chemically-homogeneous-evolution"].defaulted()) {                                        // Chemically Homogeneous Evolution
        std::tie(found, m_CheOption) = utils::GetMapKey(m_CheString, CHE_OPTION_LABEL, m_CheOption);
        if (!found) return "Unknown Chemically Homogeneous Evolution Option";
    }

    if (!m_VM["common-envelope-lambda-prescription"].defaulted()) {                                     // common envelope lambda prescription
        std::tie(found, m_CommonEnvelopeLambdaPrescription) = utils::GetMapKey(m_CommonEnvelopeLambdaPrescriptionString, CE_LAMBDA_PRESCRIPTION_LABEL, m_CommonEnvelopeLambdaPrescription);
        if (!found) return "Unknown CE Lambda Prescription";
    }

    if (!m_VM["common-envelope-mass-accretion-prescription"].defaulted()) {                             // common envelope mass accretion prescription
        std::tie(found, m_CommonEnvelopeMassAccretionPrescription) = utils::GetMapKey(m_CommonEnvelopeMassAccretionPrescriptionString, CE_ACCRETION_PRESCRIPTION_LABEL, m_CommonEnvelopeMassAccretionPrescription);
        if (!found) return "Unknown CE Mass Accretion Prescription";
    }
            
    if (!m_VM["envelope-state-prescription"].defaulted()) {                                             // envelope state prescription
        std::tie(found, m_EnvelopeStatePrescription) = utils::GetMapKey(m_EnvelopeStatePrescriptionString, ENVELOPE_STATE_PRESCRIPTION_LABEL, m_EnvelopeStatePrescription);
        if (!found) return "Unknown Envelope State Prescription";
    }

    if (!m_VM["stellar-zeta-prescription"].defaulted()) {                                               // common envelope zeta prescription
        std::tie(found, m_StellarZetaPrescription) = utils::GetMapKey(m_StellarZetaPrescriptionString, ZETA_PRESCRIPTION_LABEL, m_StellarZetaPrescription);
        if (!found) return "Unknown stellar Zeta Prescription";
    }

    if (!m_VM["eccentricity-distribution"].defaulted()) {                                               // eccentricity distribution
        std::tie(found, m_EccentricityDistribution) = utils::GetMapKey(m_EccentricityDistributionString, ECCENTRICITY_DISTRIBUTION_LABEL, m_EccentricityDistribution);
        if (!found) return "Unknown Eccentricity Distribution";
    }

    if (!m_VM["fryer-supernova-engine"].defaulted()) {                                                  // Fryer et al. 2012 supernova engine
        std::tie(found, m_FryerSupernovaEngine) = utils::GetMapKey(m_FryerSupernovaEngineString, SN_ENGINE_LABEL, m_FryerSupernovaEngine);
        if (!found) return "Unknown Fryer et al. Supernova Engine";
    }

    if (!m_VM["initial-mass-function"].defaulted()) {                                                   // initial mass function
        std::tie(found, m_InitialMassFunction) = utils::GetMapKey(m_InitialMassFunctionString, INITIAL_MASS_FUNCTION_LABEL, m_InitialMassFunction);
        if (!found) return "Unknown Initial Mass Function";
    }

    if (!m_VM["kick-direction"].defaulted()) {                                                          // kick direction
        std::tie(found, m_KickDirectionDistribution) = utils::GetMapKey(m_KickDirectionDistributionString, KICK_DIRECTION_DISTRIBUTION_LABEL, m_KickDirectionDistribution);
        if (!found) return "Unknown Kick Direction Distribution";
    }

    if (!m_VM["kick-magnitude-distribution"].defaulted()) {                                             // kick magnitude
        std::tie(found, m_KickMagnitudeDistribution) = utils::GetMapKey(m_KickMagnitudeDistributionString, KICK_MAGNITUDE_DISTRIBUTION_LABEL, m_KickMagnitudeDistribution);
        if (!found) return "Unknown Kick Magnitude Distribution";
    }

    // set values for m_KickPhi[1/2] and m_KickTheta[1/2] here
    // we now have the kick direction distribution and kick direction power (exponent) required by the user (either default or specified)

    bool phi1Defaulted   = m_VM["kick-phi-1"].defaulted();
    bool theta1Defaulted = m_VM["kick-theta-1"].defaulted();

    if (phi1Defaulted || theta1Defaulted) {
        double phi1, theta1;
        std::tie(phi1, theta1) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);
        if (phi1Defaulted  ) m_KickPhi1   = phi1;
        if (theta1Defaulted) m_KickTheta1 = theta1;
    }

    bool phi2Defaulted   = m_VM["kick-phi-2"].defaulted();
    bool theta2Defaulted = m_VM["kick-theta-2"].defaulted();

    if (phi2Defaulted || theta2Defaulted) {
        double phi2, theta2;
        std::tie(phi2, theta2) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);
        if (phi2Defaulted  ) m_KickPhi2   = phi2;
        if (theta2Defaulted) m_KickTheta2 = theta2;
    }

    if (!m_VM["logfile-delimiter"].defaulted()) {                                                       // logfile field delimiter
        std::tie(found, m_LogfileDelimiter) = utils::GetMapKey(m_LogfileDelimiterString, DELIMITERLabel, m_LogfileDelimiter);
        if (!found) return "Unknown Logfile Delimiter";
    }

    if (!m_VM["mass-loss-prescription"].defaulted()) {                                                  // mass loss prescription
        std::tie(found, m_MassLossPrescription) = utils::GetMapKey(m_MassLossPrescriptionString, MASS_LOSS_PRESCRIPTION_LABEL, m_MassLossPrescription);
        if (!found) return "Unknown Mass Loss Prescription";
    }

    if (!m_VM["mass-ratio-distribution"].defaulted()) {                                                 // mass ratio distribution
        std::tie(found, m_MassRatioDistribution) = utils::GetMapKey(m_MassRatioDistributionString, MASS_RATIO_DISTRIBUTION_LABEL, m_MassRatioDistribution);
        if (!found) return "Unknown Mass Ratio Distribution";
    }


    if (m_UseMassTransfer && !m_VM["mass-transfer-accretion-efficiency-prescription"].defaulted()) {      // mass transfer accretion efficiency prescription
        std::tie(found, m_MassTransferAccretionEfficiencyPrescription) = utils::GetMapKey(m_MassTransferAccretionEfficiencyPrescriptionString, MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL, m_MassTransferAccretionEfficiencyPrescription);
        if (!found) return "Unknown Mass Transfer Angular Momentum Loss Prescription";
    }

    if (m_UseMassTransfer && !m_VM["mass-transfer-angular-momentum-loss-prescription"].defaulted()) {     // mass transfer angular momentum loss prescription
        std::tie(found, m_MassTransferAngularMomentumLossPrescription) = utils::GetMapKey(m_MassTransferAngularMomentumLossPrescriptionString, MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL, m_MassTransferAngularMomentumLossPrescription);
        if (!found) return "Unknown Mass Transfer Angular Momentum Loss Prescription";
    }

    if (m_UseMassTransfer && !m_VM["mass-transfer-rejuvenation-prescription"].defaulted()) {              // mass transfer rejuvenation prescription
        std::tie(found, m_MassTransferRejuvenationPrescription) = utils::GetMapKey(m_MassTransferRejuvenationPrescriptionString, MT_REJUVENATION_PRESCRIPTION_LABEL, m_MassTransferRejuvenationPrescription);
        if (!found) return "Unknown Mass Transfer Rejuvenation Prescription";
    }

    if (m_UseMassTransfer && !m_VM["mass-transfer-thermal-limit-accretor"].defaulted()) {                 // mass transfer accretor thermal limit
        std::tie(found, m_MassTransferThermallyLimitedVariation) = utils::GetMapKey(m_MassTransferThermallyLimitedVariationString, MT_THERMALLY_LIMITED_VARIATION_LABEL, m_MassTransferThermallyLimitedVariation);
        if (!found) return "Unknown Mass Transfer Accretor Thermal Limit";

        if (m_MassTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::C_FACTOR) {
            m_MassTransferCParameter = m_VM["mass-transfer-thermal-limit-C"].defaulted() ? 10.0 : m_MassTransferCParameter;
        }

        if (m_MassTransferThermallyLimitedVariation == MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE) {
            m_MassTransferCParameter = m_VM["mass-transfer-thermal-limit-C"].defaulted() ? 1.0 : m_MassTransferCParameter;
        }
    }

    if (!m_VM["neutrino-mass-loss-bh-formation"].defaulted()) {                                         // neutrino mass loss assumption
        std::tie(found, m_NeutrinoMassLossAssumptionBH) = utils::GetMapKey(m_NeutrinoMassLossAssumptionBHString, NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL, m_NeutrinoMassLossAssumptionBH);
        if (!found) return "Unknown Neutrino Mass Loss Assumption";
    }

    if (!m_VM["neutron-star-equation-of-state"].defaulted()) {                                          // neutron star equation of state
        std::tie(found, m_NeutronStarEquationOfState) = utils::GetMapKey(m_NeutronStarEquationOfStateString, NS_EOSLabel, m_NeutronStarEquationOfState);
        if (!found) return "Unknown Neutron Star Equation of State";
    }

    if (!m_VM["pulsar-birth-magnetic-field-distribution"].defaulted()) {                                // pulsar birth magnetic field distribution
        std::tie(found, m_PulsarBirthMagneticFieldDistribution) = utils::GetMapKey(m_PulsarBirthMagneticFieldDistributionString, PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL, m_PulsarBirthMagneticFieldDistribution);
        if (!found) return "Unknown Pulsar Birth Magnetic Field Distribution";
    }

    if (!m_VM["pulsar-birth-spin-period-distribution"].defaulted()) {                                   // pulsar birth spin period distribution
        std::tie(found, m_PulsarBirthSpinPeriodDistribution) = utils::GetMapKey(m_PulsarBirthSpinPeriodDistributionString, PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL, m_PulsarBirthSpinPeriodDistribution);
        if (!found) return "Unknown Pulsar Birth Spin Period Distribution";
    }

    if (!m_VM["pulsational-pair-instability-prescription"].defaulted()) {                               // pulsational pair instability prescription
        std::tie(found, m_PulsationalPairInstabilityPrescription) = utils::GetMapKey(m_PulsationalPairInstabilityPrescriptionString, PPI_PRESCRIPTION_LABEL, m_PulsationalPairInstabilityPrescription);
        if (!found) return "Unknown Pulsational Pair Instability Prescription";
    }

    if (!m_VM["remnant-mass-prescription"].defaulted()) {                                               // remnant mass prescription
        std::tie(found, m_RemnantMassPrescription) = utils::GetMapKey(m_RemnantMassPrescriptionString, REMNANT_MASS_PRESCRIPTION_LABEL, m_RemnantMassPrescription);
        if (!found) return "Unknown Remnant Mass Prescription";
    }

    if (!m_VM["rotational-velocity-distribution"].defaulted()) {                                        // rotational velocity distribution
        std::tie(found, m_RotationalVelocityDistribution) = utils::GetMapKey(m_RotationalVelocityDistributionString, ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL, m_RotationalVelocityDistribution);
        if (!found) return "Unknown Rotational Velocity Distribution";
    }

    if (!m_VM["semi-major-axis-distribution"].defaulted()) {                                            // semi-major axis distribution
        std::tie(found, m_SemiMajorAxisDistribution) = utils::GetMapKey(m_SemiMajorAxisDistributionString, SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL, m_SemiMajorAxisDistribution);
        if (!found) return "Unknown Semi-Major Axis Distribution";
    }

    // constraint/value/range checks - alphabetically (where possible)

    if (!m_VM["common-envelope-alpha"].defaulted() && m_CommonEnvelopeAlpha < 0.0) return "CE alpha (--common-envelope-alpha) < 0";
    if (!m_VM["common-envelope-alpha-thermal"].defaulted() && (m_CommonEnvelopeAlphaThermal < 0.0 || m_CommonEnvelopeAlphaThermal > 1.0)) return "CE alpha thermal (--common-envelope-alpha-thermal) must be between 0 and 1";
    if (!m_VM["common-envelope-lambda-multiplier"].defaulted() && m_CommonEnvelopeLambdaMultiplier < 0.0) return "CE lambda multiplie (--common-envelope-lambda-multiplier < 0";
    if (!m_VM["common-envelope-mass-accretion-constant"].defaulted() && m_CommonEnvelopeMassAccretionConstant < 0.0) return "CE mass accretion constant (--common-envelope-mass-accretion-constant) < 0";
    if (!m_VM["common-envelope-mass-accretion-max"].defaulted() && m_CommonEnvelopeMassAccretionMax < 0.0) return "Maximum accreted mass (--common-envelope-mass-accretion-max) < 0";
    if (!m_VM["common-envelope-mass-accretion-min"].defaulted() && m_CommonEnvelopeMassAccretionMin < 0.0) return "Minimum accreted mass (--common-envelope-mass-accretion-min) < 0";

    if (m_DebugLevel < 0) return "Debug level (--debug-level) < 0";

    if (m_EccentricityDistributionMin < 0.0 || m_EccentricityDistributionMin > 1.0) return "Minimum eccentricity (--eccentricity-min) must be between 0 and 1";
    if (m_EccentricityDistributionMax < 0.0 || m_EccentricityDistributionMax > 1.0) return "Maximum eccentricity (--eccentricity-max) must be between 0 and 1";
    if (m_EccentricityDistributionMax <= m_EccentricityDistributionMin) return "Maximum eccentricity (--eccentricity-max) must be > Minimum eccentricity (--eccentricity-min)";

    if (m_InitialMassFunctionMin < 0.0) return "Minimum initial mass (--initial-mass-min) < 0";
    if (m_InitialMassFunctionMax < 0.0) return "Maximum initial mass (--initial-mass-max) < 0";
    if (m_InitialMassFunctionMax <= m_InitialMassFunctionMin) return "Maximum initial mass (--initial-mass-max) must be > Minimum initial mass (--initial-mass-min)";

    if (m_KickMagnitudeDistribution == KICK_MAGNITUDE_DISTRIBUTION::FLAT) {
        if (m_KickMagnitudeDistributionMaximum <= 0.0) return "User specified --kick-magnitude-distribution = FLAT with Maximum kick magnitude (--kick-magnitude-max) <= 0.0";
    }

    if (m_LogLevel < 0) return "Logging level (--log-level) < 0";
 
    if (m_LuminousBlueVariableFactor < 0.0) return "LBV multiplier (--luminous-blue-variable-multiplier) < 0";

    if (m_MassRatioDistributionMin < 0.0 || m_MassRatioDistributionMin > 1.0) return "Minimum mass ratio (--mass-ratio-min) must be between 0 and 1";
    if (m_MassRatioDistributionMax < 0.0 || m_MassRatioDistributionMax > 1.0) return "Maximum mass ratio (--mass-ratio-max) must be between 0 and 1";
    if (m_MassRatioDistributionMax <= m_MassRatioDistributionMin) return "Maximum mass ratio (--mass-ratio-max) must be > Minimum mass ratio (--mass-ratio-min)";

    if (m_MaxEvolutionTime <= 0.0) return "Maximum evolution time in Myr (--maxEvolutionTime) must be > 0";

    if (m_Metallicity < 0.0 || m_Metallicity > 1.0) return "Metallicity (--metallicity) should be absolute metallicity and must be between 0 and 1";

    if (m_MinimumMassSecondary < 0.0) return "Seconday minimum mass (--minimum-secondary-mass) must be >= 0";
    if (m_MinimumMassSecondary > m_InitialMassFunctionMax) return "Seconday minimum mass (--minimum-secondary-mass) must be <= Maximum initial mass (--initial-mass-max)";

    if (m_NeutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS) {
        if (m_NeutrinoMassLossValueBH < 0.0) return "Neutrino mass loss value < 0";
    }

    if (m_nBinaries <= 0) return "Number of binaries requested <= 0";

    if (m_NeutrinoMassLossAssumptionBH == NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION) {
        if (m_NeutrinoMassLossValueBH < 0.0 || m_NeutrinoMassLossValueBH > 1.0) return "Neutrino mass loss must be between 0 and 1";
    }

    if (!m_VM["outputPath"].defaulted()) {                                                              // user specified output path?
                                                                                                        // yes
        fs::path userPath = m_OutputPathString;                                                         // user-specifed path
        if (fs::is_directory(userPath)) {                                                               // valid directory?
            m_OutputPath = userPath;                                                                    // yes - set outputPath to user-specified path
        }
        else {                                                                                          // not a valid directory
            m_OutputPath = m_DefaultOutputPath;                                                         // use default path = CWD
        }
    }

    if (m_PeriodDistributionMin < 0.0) return "Minimum orbital period (--orbital-period-min) < 0";
    if (m_PeriodDistributionMax < 0.0) return "Maximum orbital period (--orbital-period-max) < 0";

    if (!m_VM["pulsar-magnetic-field-decay-timescale"].defaulted() && m_PulsarMagneticFieldDecayTimescale <= 0.0) return "Pulsar magnetic field decay timescale (--pulsar-magnetic-field-decay-timescale) <= 0";
    if (!m_VM["pulsar-magnetic-field-decay-massscale"].defaulted() && m_PulsarMagneticFieldDecayMassscale <= 0.0) return "Pulsar Magnetic field decay massscale (--pulsar-magnetic-field-decay-massscale) <= 0";

    if (m_SemiMajorAxisDistributionMin < 0.0) return "Minimum semi-major Axis (--semi-major-axis-min) < 0";
    if (m_SemiMajorAxisDistributionMax < 0.0) return "Maximum semi-major Axis (--semi-major-axis-max) < 0";

    if (m_SingleStarMassMax <= 0.0) return "Single star mass maximum (--single-star-mass-max) <= 0";
    if (m_SingleStarMassSteps > 1 && (m_SingleStarMassMax <= m_SingleStarMassMin)) return "Single star mass maximum (--single-star-mass-max) <= minimum (--single-star-mass-min)";
    if (m_SingleStarMassMin <= 0.0) return "Single star mass minimum (--single-star-mass-min) <= 0";
    if (m_SingleStarMassSteps <= 0) return "Single star mass steps (--single-star-mass-steps) <= 0";

    // check illegal combinations

    if (m_SingleStar && m_BSEswitchLog) return "--BSEswitchLog does not apply to Single Star evolution";
    if (!m_SingleStar && m_SSEswitchLog) return "--SSEswitchLog does not apply to Binary Star evolution";

 
    if (m_WolfRayetFactor < 0.0) return "WR multiplier (--wolf-rayet-multiplier) < 0";

    return "";  // all good
}


std::string Options::OptionValues::Initialise(const po::variables_map p_VM) {

    // First set all options to their default values

    // flags

    m_AllowRLOFAtBirth                                              = false;
    m_AllowTouchingAtBirth                                          = false;

    m_DebugToFile                                                   = false;
    m_ErrorsToFile                                                  = false;

    m_EnableWarnings                                                = false;

    m_SingleStar                                                    = false;

	m_BeBinaries                                                    = false;
    m_EvolvePulsars                                                 = false;
	m_EvolveUnboundSystems                                          = false;

    m_DetailedOutput                                                = false;
    m_PopulationDataPrinting                                        = false;
    m_PrintBoolAsString                                             = false;
    m_Quiet                                                         = false;
    m_RlofPrinting                                                  = false;
    m_BSEswitchLog                                                  = false;
    m_SSEswitchLog                                                  = false;

    m_nBatchesUsed                                                  = -1;


    // Public population synthesis variables
    m_nBinaries                                                     = 10;
    m_FixedRandomSeed                                               = false;                                                // TRUE if --random-seed is passed on command line
    m_RandomSeed                                                    = 0;

    // Specify how long to evolve binaries for
    m_MaxEvolutionTime                                              = 13700.0;
    m_MaxNumberOfTimestepIterations                                 = 99999;


    // Initial mass options
    m_InitialMassFunction                                           = INITIAL_MASS_FUNCTION::KROUPA;
    m_InitialMassFunctionString                                     = INITIAL_MASS_FUNCTION_LABEL.at(m_InitialMassFunction);
    m_InitialMassFunctionMin                                        = 8.0;
    m_InitialMassFunctionMax                                        = 100.0; 
    m_InitialMassFunctionPower                                      = -2.3;


    // Initial mass ratios
    m_MassRatioDistribution                                         = MASS_RATIO_DISTRIBUTION::FLAT;                        // Most likely want FLAT or SANA2012
    m_MassRatioDistributionString                                   = MASS_RATIO_DISTRIBUTION_LABEL.at(m_MassRatioDistribution);
    m_MassRatioDistributionMin                                      = 0.0;
    m_MassRatioDistributionMax                                      = 1.0;

    m_MinimumMassSecondary                                          = 0.0;


    // Initial orbit options
    m_SemiMajorAxisDistribution                                     = SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG;              // Most likely want FLATINLOG or SANA2012
    m_SemiMajorAxisDistributionString                               = SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL.at(SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG);
    m_SemiMajorAxisDistributionMin                                  = 0.1;
    m_SemiMajorAxisDistributionMax                                  = 1000.0;
    m_SemiMajorAxisDistributionPower                                = -1.0; 


    // Initial orbital period
    m_PeriodDistributionMin                                         = 1.1;
    m_PeriodDistributionMax                                         = 1000.0;

    // Eccentricity
    m_EccentricityDistribution                                      = ECCENTRICITY_DISTRIBUTION::ZERO; 
    m_EccentricityDistributionString                                = ECCENTRICITY_DISTRIBUTION_LABEL.at(m_EccentricityDistribution);
    m_EccentricityDistributionMin                                   = 0.0;
    m_EccentricityDistributionMax                                   = 1.0;

    // Kick options
    m_KickMagnitudeDistribution                                     = KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN;
    m_KickMagnitudeDistributionString                               = KICK_MAGNITUDE_DISTRIBUTION_LABEL.at(m_KickMagnitudeDistribution);
    m_KickMagnitudeDistributionSigmaCCSN_NS                         = 250;
    m_KickMagnitudeDistributionSigmaCCSN_BH                         = 250;
    m_KickMagnitudeDistributionMaximum                              = -1.0; 
    m_KickMagnitudeDistributionSigmaForECSN                         = 30.0;
    m_KickMagnitudeDistributionSigmaForUSSN   	                    = 30.0;
	m_KickScalingFactor						                        = 1.0;

    // Kick direction option
    m_KickDirectionDistribution                                     = KICK_DIRECTION_DISTRIBUTION::ISOTROPIC;
    m_KickDirectionDistributionString                               = KICK_DIRECTION_DISTRIBUTION_LABEL.at(m_KickDirectionDistribution);
    m_KickDirectionPower                                            = 0.0;

    // Kick magnitude
    m_KickMagnitude                                                 = 0.0;
    m_KickMagnitude1                                                = 0.0;
    m_KickMagnitude2                                                = 0.0;                               

    // Kick magnitude random number (used to draw kick magnitude if necessary)
    m_KickMagnitudeRandom                                           = RAND->Random();
    m_KickMagnitudeRandom1                                          = RAND->Random();
    m_KickMagnitudeRandom2                                          = RAND->Random();

    // Mean anomaly
    m_KickMeanAnomaly1                                              = RAND->Random(0.0, _2_PI);
    m_KickMeanAnomaly2                                              = RAND->Random(0.0, _2_PI);

    // Phi
    m_KickPhi1                                                      = 0.0;                                              // actual value set later
    m_KickPhi2                                                      = 0.0;                                              // actual value set later

    // Theta
    m_KickTheta1                                                    = 0.0;                                              // actual value set later 
    m_KickTheta2                                                    = 0.0;                                              // actual value set later

    // Black hole kicks
    m_BlackHoleKicksOption                                          = BLACK_HOLE_KICK_OPTION::FALLBACK;
    m_BlackHoleKicksOptionString                                    = BLACK_HOLE_KICK_OPTION_LABEL.at(m_BlackHoleKicksOption);


    // Chemically Homogeneous Evolution
    m_CheOption                                                     = CHE_OPTION::NONE;
    m_CheString                                                     = CHE_OPTION_LABEL.at(m_CheOption);


    // Supernova remnant mass prescription options
    m_RemnantMassPrescription                                       = REMNANT_MASS_PRESCRIPTION::FRYER2012;
    m_RemnantMassPrescriptionString                                 = REMNANT_MASS_PRESCRIPTION_LABEL.at(m_RemnantMassPrescription);

    m_FryerSupernovaEngine                                          = SN_ENGINE::DELAYED;
    m_FryerSupernovaEngineString                                    = SN_ENGINE_LABEL.at(m_FryerSupernovaEngine);

    m_NeutrinoMassLossAssumptionBH                                  = NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION;
    m_NeutrinoMassLossAssumptionBHString                            = NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL.at(m_NeutrinoMassLossAssumptionBH);
    m_NeutrinoMassLossValueBH                                       = 0.1;


    // Fixed uk options
    m_UseFixedUK                                                    = false;
    m_FixedUK                                                       = -1.0;


    // Pair instability and pulsational pair instability mass loss
    m_UsePairInstabilitySupernovae                                  = false;
    m_PairInstabilityLowerLimit                                     = 60.0;                                                 // Belczynski+ 2016 is 65 Msol
    m_PairInstabilityUpperLimit                                     = 135.0;                                                // Belczynski+ 2016 is 135 Msol

    m_UsePulsationalPairInstability                                 = false;
    m_PulsationalPairInstabilityLowerLimit                          = 35.0;                                                 // Belczynski+ 2016 is 45 Msol
    m_PulsationalPairInstabilityUpperLimit                          = 60.0;                                                 // Belczynski+ 2016 is 65 Msol

    m_PulsationalPairInstabilityPrescription                        = PPI_PRESCRIPTION::COMPAS;
    m_PulsationalPairInstabilityPrescriptionString                  = PPI_PRESCRIPTION_LABEL.at(m_PulsationalPairInstabilityPrescription);

	m_MaximumNeutronStarMass                                        = 3.0;                                                  // StarTrack is 3.0
    
    m_mCBUR1                                                        = MCBUR1HURLEY;                                         // MHurley value, Fryer+ and Belczynski+ use 1.83


    // Output path
    m_OutputPathString                                              = ".";
    m_DefaultOutputPath                                             = boost::filesystem::current_path();
    m_OutputPath                                                    = m_DefaultOutputPath;
    m_OutputContainerName                                           = DEFAULT_OUTPUT_CONTAINER_NAME;
    

    // Mass loss options
    m_UseMassLoss                                                   = false;

    m_MassLossPrescription                                          = MASS_LOSS_PRESCRIPTION::VINK;
    m_MassLossPrescriptionString                                    = MASS_LOSS_PRESCRIPTION_LABEL.at(m_MassLossPrescription);


    // Wind mass loss multiplicitive constants
    m_LuminousBlueVariableFactor                                    = 1.5;
    m_WolfRayetFactor                                               = 1.0;


    // Mass transfer options
    m_UseMassTransfer                                               = true;
	m_CirculariseBinaryDuringMassTransfer         	                = false;
	m_AngularMomentumConservationDuringCircularisation              = false;

    // Case BB/BC mass transfer stability prescription
    m_CaseBBStabilityPrescription                                   = CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE;
    m_CaseBBStabilityPrescriptionString                             = CASE_BB_STABILITY_PRESCRIPTION_LABEL.at(m_CaseBBStabilityPrescription);

    // Options adaptive Roche Lobe Overflow prescription
    m_MassTransferAdaptiveAlphaParameter                            = 0.5;
    m_MaxPercentageAdaptiveMassTransfer                             = 0.01;

    // Options for mass transfer accretion efficiency
    m_MassTransferAccretionEfficiencyPrescription                   = MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED;
    m_MassTransferAccretionEfficiencyPrescriptionString             = MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL.at(m_MassTransferAccretionEfficiencyPrescription);

    m_MassTransferFractionAccreted                                  = 1.0;
    m_MassTransferCParameter                                        = 10.0;
    m_EddingtonAccretionFactor                                      = 1;                                                    // >1 is super-eddington, 0 is no accretion

    // Mass transfer thermally limited options
	m_MassTransferThermallyLimitedVariation                         = MT_THERMALLY_LIMITED_VARIATION::C_FACTOR;
	m_MassTransferThermallyLimitedVariationString                   = MT_THERMALLY_LIMITED_VARIATION_LABEL.at(m_MassTransferThermallyLimitedVariation);

    // Mass transfer angular momentum loss prescription options
    m_MassTransferJloss                                             = 1.0;
    m_MassTransferAngularMomentumLossPrescription                   = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION;
    m_MassTransferAngularMomentumLossPrescriptionString             = MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL.at(m_MassTransferAngularMomentumLossPrescription);

    // Mass transfer rejuvenation prescriptions
    m_MassTransferRejuvenationPrescription                          = MT_REJUVENATION_PRESCRIPTION::NONE;
    m_MassTransferRejuvenationPrescriptionString                    = MT_REJUVENATION_PRESCRIPTION_LABEL.at(m_MassTransferRejuvenationPrescription);

    // Mass transfer critical mass ratios
    m_MassTransferCriticalMassRatioMSLowMass                        = false;
    m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor   = 1.44;                                                 // Claeys+ 2014 = 1.44
    m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor      = 1.0;                                                  // Claeys+ 2014 = 1.0

    m_MassTransferCriticalMassRatioMSHighMass                       = false;
    m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor  = 0.625;                                                // Claeys+ 2014 = 0.625
    m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor     = 0.0;

    m_MassTransferCriticalMassRatioHG                               = false;
    m_MassTransferCriticalMassRatioHGNonDegenerateAccretor          = 0.40;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHGDegenerateAccretor             = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioGiant                            = false;
    m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor       = 0.0;
    m_MassTransferCriticalMassRatioGiantDegenerateAccretor          = 0.87;                                                 // Claeys+ 2014 = 0.81

    m_MassTransferCriticalMassRatioHeliumMS                         = false;
    m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor    = 0.625;
    m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor       = 0.0;

    m_MassTransferCriticalMassRatioHeliumHG                         = false;
    m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor    = 0.25;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor       = 0.21;                                                 // Claeys+ 2014 = 0.21

    m_MassTransferCriticalMassRatioHeliumGiant                      = false;
    m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor = 1.28;                                                 // Claeys+ 2014 = 0.25
    m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor    = 0.87;

    m_MassTransferCriticalMassRatioWhiteDwarf                       = false;
	m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor  = 0.0;
    m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor     = 1.6;                                                  // Claeys+ 2014 = 1.6


    // Common Envelope options
    m_CommonEnvelopeAlpha                                           = 1.0;
    m_CommonEnvelopeLambda                                          = 0.1;
	m_CommonEnvelopeSlopeKruckow                                    = -4.0 / 5.0;
	m_CommonEnvelopeAlphaThermal                                    = 1.0;
    m_CommonEnvelopeLambdaMultiplier                                = 1.0;
    m_AllowMainSequenceStarToSurviveCommonEnvelope                  = false;

    // Prescription for envelope state (radiative or convective)
    m_EnvelopeStatePrescription                                     = ENVELOPE_STATE_PRESCRIPTION::LEGACY;
    m_EnvelopeStatePrescriptionString                               = ENVELOPE_STATE_PRESCRIPTION_LABEL.at(m_EnvelopeStatePrescription);

    // Accretion during common envelope
    m_CommonEnvelopeMassAccretionPrescription                       = CE_ACCRETION_PRESCRIPTION::ZERO;
    m_CommonEnvelopeMassAccretionPrescriptionString                 = CE_ACCRETION_PRESCRIPTION_LABEL.at(m_CommonEnvelopeMassAccretionPrescription);
    
    m_CommonEnvelopeMassAccretionMin                                = 0.04;
    m_CommonEnvelopeMassAccretionMax                                = 0.1;
    m_CommonEnvelopeMassAccretionConstant                           = 0.0;

	// Common envelope lambda prescription
	m_CommonEnvelopeLambdaPrescription                              = CE_LAMBDA_PRESCRIPTION::NANJING;
	m_CommonEnvelopeLambdaPrescriptionString                        = CE_LAMBDA_PRESCRIPTION_LABEL.at(m_CommonEnvelopeLambdaPrescription);

	// Common envelope Nandez and Ivanova energy formalism
	m_RevisedEnergyFormalismNandezIvanova	                        = false;
	m_MaximumMassDonorNandezIvanova                                 = 2.0;
	m_CommonEnvelopeRecombinationEnergyDensity                      = 1.5E13;


    // Adaptive Importance Sampling options
    m_AISexploratoryPhase                                           = false;
    m_AISDCOtype                                                    = AIS_DCO::ALL;
    m_AISDCOtypeString                                              = AIS_DCO_LABEL.at(AIS_DCO::ALL);
    m_AIShubble                                                     = false;
    m_AISpessimistic                                                = false;
    m_AISrefinementPhase                                            = false;
    m_AISrlof                                                       = false;
    m_KappaGaussians                                                = 2;


	// Zetas
	m_StellarZetaPrescription                                       = ZETA_PRESCRIPTION::SOBERMAN;
	m_StellarZetaPrescriptionString                                 = ZETA_PRESCRIPTION_LABEL.at(m_StellarZetaPrescription);
	m_ZetaAdiabaticArbitrary                                        = 10000.0;                                              // large value favours stable MT
    m_ZetaMainSequence 	                                            = 2.0;
	m_ZetaRadiativeEnvelopeGiant	                                = 6.5;


    // Metallicity options
    m_Metallicity                                                   = ZSOL;
    m_FixedMetallicity                                              = true;


    // Neutron star equation of state
    m_NeutronStarEquationOfState                                    = NS_EOS::SSE;
    m_NeutronStarEquationOfStateString                              = NS_EOSLabel.at(NS_EOS::SSE);


    // Pulsar birth magnetic field distribution
    m_PulsarBirthMagneticFieldDistribution                          = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO;
    m_PulsarBirthMagneticFieldDistributionString                    = PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL.at(m_PulsarBirthMagneticFieldDistribution);
    m_PulsarBirthMagneticFieldDistributionMin                       = 11.0;
    m_PulsarBirthMagneticFieldDistributionMax                       = 13.0;


    // Pulsar birth spin period distribution string
    m_PulsarBirthSpinPeriodDistribution                             = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO;
    m_PulsarBirthSpinPeriodDistributionString                       = PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL.at(m_PulsarBirthSpinPeriodDistribution);
    m_PulsarBirthSpinPeriodDistributionMin                          = 0.0;
    m_PulsarBirthSpinPeriodDistributionMax                          = 100.0;

    m_PulsarMagneticFieldDecayTimescale                             = 1000.0;
    m_PulsarMagneticFieldDecayMassscale                             = 0.025;
    m_PulsarLog10MinimumMagneticField                               = 8.0;


    // Rotational velocity distribution options
    m_RotationalVelocityDistribution                                = ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO;
    m_RotationalVelocityDistributionString                          = ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL.at(m_RotationalVelocityDistribution);


	// grids

	m_GridFilename                                                  = "";


    // debug and logging options

    m_DebugLevel                                                    = 0;
    m_DebugClasses.clear();

    m_LogLevel                                                      = 0;
    m_LogClasses.clear();
    
    m_LogfileNamePrefix                                             = "";
    m_LogfileDelimiter                                              = DELIMITER::TAB;
    m_LogfileDelimiterString                                        = DELIMITERLabel.at(m_LogfileDelimiter);

    m_LogfileDefinitionsFilename                                    = "";


    // SSE options
    m_LogfileSSEParameters                                          = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_PARAMETERS));
    m_LogfileSSESupernova                                           = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVA));
    m_LogfileSSESwitchLog                                           = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG));

    m_SingleStarMassSteps                                           = 100;
    m_SingleStarMassMin                                             = 5.0;
    m_SingleStarMassMax                                             = 100.0;


    // BSE options
    m_LogfileBSESystemParameters                                    = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS));
    m_LogfileBSEDetailedOutput                                      = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT));
    m_LogfileBSEDoubleCompactObjects                                = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS));
    m_LogfileBSESupernovae                                          = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE));
    m_LogfileBSECommonEnvelopes                                     = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_COMMON_ENVELOPES));
    m_LogfileBSERLOFParameters                                      = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS));
    m_LogfileBSEBeBinaries                                          = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES));
    m_LogfileBSEPulsarEvolution                                     = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION));
    m_LogfileBSESwitchLog                                           = get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG));


    // set variable map
    m_VM = p_VM;

    // now check user-supplied options and set values as appropriate
    return CheckAndSetOptions();
}


bool Options::Initialise(int p_OptionCount, char *p_OptionStrings[]) {

    bool ok = true;                                                                                                         // status - unless something changes

    try {

        // initialise program options - these are the options that are specified only on the 
        // commandline and cannot be specified in a grid file on a per object (star/binary) basis
        po::options_description programInstanceOptions("ProgramOptions");
        ok = ProgramOptions(&m_CmdLine, &programInstanceOptions);                                                  // boost options descriptions object for per program options
        // CHECK OK HERE!!!!!!!!!!!!!!!!!!!!

        // initialise program options - these are the options that can be specified on the 
        // commandline, but can also be specified in a grid file on a per object (star/binary)
        // basis.  The values of options specified in a grid file take precedence over the
        // values of the same options specified on the commandline - only for the object (star/
        // binary) corresponding to the grid file record.
        po::options_description objectInstanceOptions("ObjectOptions");
        ok = ObjectOptions(&m_CmdLine, &objectInstanceOptions);                                                   // boost options descriptions object for per object (star/binary) options
        // CHECK OK HERE!!!!!!!!!!!!!!!!!!!!

        // populate the commandline options
        // this stays static throughout the life of the program

        po::options_description cmdlineOptions;                                                                             // boost options descriptions object for options available on commandline
        cmdlineOptions.add(programInstanceOptions).add(objectInstanceOptions);                                              // both program and object options are available on the commandline

        po::variables_map vm;                                                                                               // boost variables map - populated by parse

        po::parsed_options const parsedOptions = po::parse_command_line(p_OptionCount, p_OptionStrings, cmdlineOptions);    // parse user-supplied options
        po::store(parsedOptions, vm);                                                                                       // store parsed options into variable map

        std::string errStr = m_CmdLine.Initialise(vm);                                                                                      // initialise commandline option variables and variable map
        // CHECK ERRSTR HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        m_OptionsDetails = ProgramOptionDetails(&m_CmdLine, vm);                                                            // construct options details string for output to run details file

    } catch (...) {                                                                                                         // unhandled exception - something wrong with boost
        ok = false;                                                                                                         // set status
    }
    
    return ok;
}


bool Options::InitialiseObject(const std::string p_OptionsString) {

    // parse the option string

    std::vector<std::string> parsedStrings;                                     // parsed option strings

    size_t start      = 0;                                                      // start position of parsed option string
    size_t end        = 0;                                                      // end position of parsed option strinf
    std::string delim = " -";                                                   // delimiter
    while (end != string::npos) {                                               // iterate over input string
        end = p_OptionsString.find(delim, start);                               // find delimiter
        parsedStrings.push_back(p_OptionsString.substr(start, end - start));    // grab option string
        start = end + delim.length();                                           // new start position
    }
    
    std::vector<char const*> args {"placeHolder"};                              // place-holder - boost expects command name as argv[0]
    for (auto& arg : parsedStrings)                                             // iterate over the parsed strings
        args.push_back(arg.c_str());                                            // and grab the c_str



    m_EvolvingObject = m_CmdLine;                                               // copy commandline options

    bool ok;

    // initialise program options - these are the options that can be specified on the 
    // commandline, but can also be specified in a grid file on a per object (star/binary)
    // basis.  The values of options specified in a grid file take precedence over the
    // values of the same options specified on the commandline - only for the object (star/
    // binary) corresponding to the grid file record.
    OptionValues thisObject;
    po::options_description objectInstanceOptions("thisObjectOptions");
    ok = ObjectOptions(&thisObject, &objectInstanceOptions);                                                   // boost options descriptions object for per object (star/binary) options
    // CHECK OK HERE!!!!!!!!!!!!!!!!!!!!

    // populate the commandline options
    // this stays static throughout the life of the program

    po::options_description objectOptions;                                                                             // boost options descriptions object for options available on commandline
    objectOptions.add(objectInstanceOptions);                                              // both program and object options are available on the commandline

    po::variables_map vm;                                                                                               // boost variables map - populated by parse

    po::parsed_options const parsedOptions = po::parse_command_line(args.size(), args.data(), objectOptions);    // parse user-supplied options
    po::store(parsedOptions, vm);                                                                                       // store parsed options into variable map

//    std::string errStr = m_CmdLine.Initialise(vm);                                                                                      // initialise commandline option variables and variable map
    // CHECK ERRSTR HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        // initialise any options that rely on the value of other options
// ONLY IF NECESSARY - CHECK WHAT OPTIONS ARE SPECIFIED HERE!!!!
//        std::tie(m_KickTheta, m_KickPhi) = utils::DrawKickDirection(m_KickDirectionDistribution, m_KickDirectionPower);


    return true;


}


/*
 * Determine the value of the requested program option
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function handles properties of type PROGRAM_OPTION
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  The type of the
 * functional return is a tuple: std::tuple<bool, COMPAS_VARIABLE_TYPE>.  This type
 * is COMPAS_VARIABLE by typedef.
 *
 * The bool returned indicates whether the property value was retrieved ok: true = yes, fales = no
 * The COMPAS_VARIABLE_TYPE variable returned is a boost variant variable, the value of which is the
 * value of the underlying primitive variable.
 *
 *
 * COMPAS_VARIABLE OptionValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE Options::OptionValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = true;                                                                                                     // status - unless a problem occurs

    COMPAS_VARIABLE_TYPE value;                                                                                         // default property value

    PROGRAM_OPTION property = boost::get<PROGRAM_OPTION>(p_Property);                                                   // get property
                                                                                                                        // get property value
    switch (property) {

        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH:  value = KickMagnitudeDistributionSigmaCCSN_BH(); break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS:  value = KickMagnitudeDistributionSigmaCCSN_NS(); break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN: value = KickMagnitudeDistributionSigmaForECSN(); break;
        case PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN: value = KickMagnitudeDistributionSigmaForUSSN(); break;
        case PROGRAM_OPTION::RANDOM_SEED:                                value = RandomSeed();                            break;

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            std::cerr << ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION) << std::endl;                                           // show warning (don't have logging or errors here...)
    }

    return std::make_tuple(ok, value);
}


/*
 * Open the SSE grid file, read and parse the header record   FIX THIS DOCUMENTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
 *
 * std::tuple<int, std::vector<std::string>, std::string> OpenSSEGridFile(std::ifstream &p_Grid, const std::string p_Filename)
 *
 * @param   [IN/OUT]    p_Grid                  The std::ifstream to which the file should be opened
 * @param   [IN]        p_Filename              The filename of the Grid file
 * @return                                      tuple containing the line number of the next line to be read and a vector of header strings
 */
void Options::OpenGridFile(const std::string p_GridFilename) {

    std::ifstream GridFile; // JRFIX !!!!!!!!!!!!!!!

    int tokenCount         = 0;                                                                                                 // token count - to check if header should be present
//    int currentPos         = GridFile.tellg();                                                                                    // current file position - to rewind if no header

    int lineNo = 1;                                                                                                             // line number - for error messages

    if (!p_GridFilename.empty()) {                                                                                     // have grid filename?

        GridFile.open(p_GridFilename);                                                                                   // yes - open the file
        if (GridFile.fail()) {                                                                                                    // open ok?
            ANNOUNCE(ERR_MSG(ERROR::FILE_OPEN_ERROR) << " " << p_GridFilename);                                             // no - show error
        }
        else {                                                                                                                  // file open ok

            bool readFailed  = false;                                                                                           // record read ok?
            bool emptyRecord = true;                                                                                            // record empty?
            std::string record;
            while (emptyRecord && !readFailed) {                                                                                // get the header - first non-empty record

//                currentPos = GridFile.tellg();                                                                                    // current file position
                if (std::getline(GridFile, record)) {                                                                             // read next record
                    readFailed = false;                                                                                         // ok

                    size_t hashPos = record.find("#");                                                                          // find first occurrence of "#"
                    if (hashPos != std::string::npos) record.erase(hashPos, record.size() - hashPos);                           // if "#" found, prune it and everything after it (ignore comments)

                    while (record.size() > 0 && (record[0] == ' ' || record[0] == '\t')) record.erase(0, 1);                    // strip any leading ' ' or TAB ('\t') characters from the record

                    while (record.size() > 0               &&
                          (record[record.size()-1] == '\r' ||
                           record[record.size()-1] == '\n' ||
                           record[record.size()-1] == '\t' ||
                           record[record.size()-1] == ' '  )) record.erase(record.size()-1, 1);                                 // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the record

                    // check if record is effectively empty
                    // significant characters are characters other than {' ', ',', '\t', '\r', '\n'}
                    int significantChars = 0;                                                                                   // initially none
                    for (size_t i = 0; i < record.length(); i++) {                                                              // for each character in the record
                        if (record[i] != ' '  && record[i] != ',' &&  record[i] != '\t' && record[i] != '\r' && record[i] != '\n') {
                            significantChars++;
                            break;                                                                                              // only need one to be not (effectively) empty
                        }
                    }
                    if (significantChars == 0) record = "";                                                                     // make it empty if nothing significant

                    if (record.empty()) {
                        lineNo++;                                                                                               // increment line number
                        continue;                                                                                               // skip empty record
                    }

                    emptyRecord = false;                                                                                        // have non-empty record

                    record = utils::ToUpper(record);                                                                            // upshift - we ignore case here

                    std::stringstream recordSS(record);                                                                         // open stream on record
                    std::string token;
                    while (std::getline(recordSS, token, ',')) {                                                                // get token from record read

                        tokenCount++;                                                                                           // increment token count - even if it's empty, it's a token

                        while (token.size() > 0 && (token[0] == ' ' || token[0] == '\t')) token.erase(0, 1);                    // strip any leading ' ' or TAB ('\t') characters from the token

                        while (token.size() > 0              &&
                              (token[token.size()-1] == '\r' ||
                               token[token.size()-1] == '\n' ||
                               token[token.size()-1] == '\t' ||
                               token[token.size()-1] == ' ')) token.erase(token.size()-1, 1);                                   // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the token

                        if (token.empty()) {                                                                                    // empty token?
//                            ANNOUNCE(ERR_MSG(ERROR::GRID_FILE_EMPTY_HEADER));                                                        // show error
                            continue;                                                                                           // ignore it
                        }

                        // have a token - process it                      
                        // tokens here are expected to be of the form:
                        //     'option' - for specifying BOOLEAN options with value = TRUE
                        //     'option = value' (spaces around '=' not significant) - for fully specifying the opion and its value

                        // first, split into up to two tokens: optionName and, otionally, optionValue
                        // if '=' is present in the stringm then 'optionValue' must also be present
                        // if '=' is not present in the string, then:
                        //     'optionValue' is assumed not to be present in the string
                        //     'optionName' must specify a BOOLEAN option, and its value is assumed to be TRUE
                        // BOOLEAN options can be present, and their valuesspecified (so '=' will be present)
                        // BOOLEAN option values can be specified as:
                        //      'TRUE' (case not significant), or '1', or
                        //      'FASLE' (case not significant), or '0'

                        std::string optionName  = "";
                        std::string optionValue = "";

                        size_t loc = token.find("=");                                                                           // find '='
                        if (loc == string::npos) {                                                                              // '=' not found
                            optionName  = token;                                                                                // option name - boolean option
                            optionValue = "TRUE";                                                                               // default value for boolean options
                        }
                        else {                                                                                                  // '=' found
                            optionName  = token.substr (0, loc);                                                                // option name
                            optionValue = token.substr (0, loc);                                                                // option value
                        }
                        optionName  = utils::trim(optionName);                                                                  // trim leading and trailing whitespace
                        optionValue = utils::trim(optionName);                                                                  // trim leading and trailing whitespace

                        // neither optionName nor optionValue may have spaces
                        loc = optionName.find(" ");                                                                             // look for ' '
                        if (loc != string::npos) {                                                                              // found - that's not ok
//                            ERROR error = ERROR::GRID_FILE_INVALID_OPTION_NAME;                                                 // set error
//                            SAY(ERR_MSG(error));                                                                                // announce error
                            break;                                                                                              // and bail out
                        }

                        loc = optionValue.find(" ");                                                                            // look for ' '
                        if (loc != string::npos) {                                                                              // found - that's not ok
//                            ERROR error = ERROR::GRID_FILE_INVALID_OPTION_VALUE;                                                // set error
//                            SAY(ERR_MSG(error));                                                                                // announce error
                            break;                                                                                              // and bail out
                        }

                        // have option name and option value (if required) - process them

                        // we expect the supplied option name to be one of the COMPAS defined option names
                        bool validOptionName = false;                                                                           // default is not valid
                        /*
                        for (po::variables_map::const_iterator it = vm.begin(); it != vm.end(); it++) {                     // for all options
                            std::string thisOptionName = utils::ToUpper(it->first);                                             // upshifted option name
//                            po::option_description const& opt = desc.find(it->first.c_str(), false, false, false);
//                            std::string thisOptionShortName = "!"; //opt.canonical_display_name(cls::allow_dash_for_short);
//                            std::cout << "JRPINT, thisOptionName = " << thisOptionName << ", short: " << thisOptionShortName << "\n";
                            if (optionName == thisOptionName) {                                                                 // supplied option name = this option name?
                                validOptionName = true;                                                                         // yes - option name is valid
                                optionName = it->first;                                                                         // set name to actual COMPAS option name
                                break;                                                                                          // and we're done
                            }
                        }
                        */
                        if (!validOptionName) {                                                                                 // supplied option name a valid COMPAS option name?
//                            ERROR error = ERROR::GRID_FILE_INVALID_OPTION_NAME;                                                 // no - set error (we could just warn and continue, but I think error and stop is better)
//                            SAY(ERR_MSG(error));                                                                                // announce error
                            break;                                                                                              // and bail out
                        }






                    }

                    lineNo++;                                                                                                   // increment line number
                }
                else readFailed = true;                                                                                         // read failed - may be EOF
            }
        }
    }

}
