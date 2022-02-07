BSE system parameters
=====================

Default record definition for the BSE System Parameters log file::

    const ANY_PROPERTY_VECTOR BSE_SYSTEM_PARAMETERS_REC = {
        BINARY_PROPERTY::RANDOM_SEED,                                    
        STAR_1_PROPERTY::MZAMS,                                          
        STAR_2_PROPERTY::MZAMS,                                          
        STAR_1_PROPERTY::INITIAL_RADIUS,                                  
        STAR_2_PROPERTY::INITIAL_RADIUS,                                  
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,                         
        BINARY_PROPERTY::ECCENTRICITY_INITIAL,                            
        STAR_1_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,          
        STAR_2_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,          
        STAR_1_PROPERTY::OMEGA_ZAMS,                                      
        STAR_2_PROPERTY::OMEGA_ZAMS,                                      
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,        
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,        
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,       
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,       
        PROGRAM_OPTION::LBV_FACTOR,                                       
        PROGRAM_OPTION::WR_FACTOR,                                        
        PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA,                            
        STAR_1_PROPERTY::METALLICITY,                                     
        STAR_2_PROPERTY::METALLICITY,                                     
        BINARY_PROPERTY::UNBOUND,                                         
        BINARY_PROPERTY::STELLAR_MERGER,                                  
        BINARY_PROPERTY::STELLAR_MERGER_AT_BIRTH,                         
        BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH,                    
        STAR_1_PROPERTY::INITIAL_STELLAR_TYPE,                            
        STAR_1_PROPERTY::STELLAR_TYPE,                                    
        STAR_2_PROPERTY::INITIAL_STELLAR_TYPE,                            
        STAR_2_PROPERTY::STELLAR_TYPE,                                    
        STAR_1_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,            
        STAR_2_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,            
        BINARY_PROPERTY::ERROR,                                       
        PROGRAM_OPTION::NOTES
    };
