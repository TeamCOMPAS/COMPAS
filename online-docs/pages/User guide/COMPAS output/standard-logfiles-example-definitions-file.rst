Example log file definitions file
=================================

Following is an example log file definitions file. COMPAS can be configured to use this file via the ``--logfile-definitions`` program option.

This file (``COMPAS_Output_Definitions.txt``) is also delivered as part of the COMPAS github repository.

::

    # sample standard log file specifications file
    
    # the ’#’ character and anything following it on a single line is considered a comment
    # (so, lines starting with ’#’ are comment lines)
    
    # case is not significant
    # specifications can span several lines
    # specifications for the same log file are cumulative
    # if a log file is not specified in this file, the default specification is used
    
    
    # SSE Parameters
    
    # start with the default SSE Parameters specification and add ENV_MASS
    
    sse_sysparms_rec += { STAR_PROPERTY::ENV_MASS }
    
    # take the updated SSE Parameters specification and add ANGULAR_MOMENTUM
    
    sse_sysparms_rec += { STAR_PROPERTY::ANGULAR_MOMENTUM }
    
    # take the updated SSE Parameters specification and subtract MASS_0 and MDOT
    
    sse_sysparms_rec -= { STAR_PROPERTY::MASS_0, STAR_PROPERTY::MDOT }
    
    
    # BSE System Parameters
    
    bse_sysparms_rec = {                  # set the BSE System Parameters specification to:
        BINARY_PROPERTY::ID,              # ID of the binary
        BINARY_PROPERTY::RANDOM_SEED,     # RANDOM_SEED for the binary
        STAR_1_PROPERTY::MZAMS,           # MZAMS for Star1
        STAR_2_PROPERTY::MZAMS            # MZAMS for Star2
    }
    
    # ADD to the BSE System Parameters specification:
    # SEMI_MAJOR_AXIS_INITIAL for the binary
    # ECCENTRICITY_INITIAL for the binary
    # SUPERNOVA_THETA for Star1 and SUPERNOVA_PHI for Star1
    
    bse_sysparms_rec += {
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,
        BINARY_PROPERTY::ECCENTRICITY_INITIAL,
        STAR_1_PROPERTY::SUPERNOVA_THETA, STAR_1_PROPERTY::SUPERNOVA_PHI
    }
    
    bse_sysparms_rec += {                          # ADD to the BSE System Parameters specification:
        SUPERNOVA_PROPERTY::IS_ECSN,               # IS_ECSN for the supernova star
        SUPERNOVA_PROPERTY::IS_SN,                 # IS_SN for the supernova star
        SUPERNOVA_PROPERTY::IS_USSN,               # IS_USSN for the supernova star
        SUPERNOVA_PROPERTY::EXPERIENCED_PISN,      # EXPERIENCED_PISN for the supernova star
        SUPERNOVA_PROPERTY::EXPERIENCED_PPISN,     # EXPERIENCED_PPISN for the supernova star
        BINARY_PROPERTY::UNBOUND,                  # UNBOUND for the binary
        SUPERNOVA_PROPERTY::MZAMS,                 # MZAMS for the supernova star
        COMPANION_PROPERTY::MZAMS                  # MZAMS for the companion star
    }
    
    # SUBTRACT from the BSE System Parameters specification:
    # RANDOM_SEED for the binary
    # ID for the binary
    # all ANNOTATIONS
    
    bse_sysparms_rec -= {                 # SUBTRACT from the BSE System Parameters specification:
        BINARY_PROPERTY::RANDOM_SEED,     # RANDOM_SEED for the binary
        BINARY_PROPERTY::ID,              # ID for the binary
        PROGRAM_OPTION::NOTES             # all ANNOTATIONS
    }
    
    bse_sysparms_rec += { PROGRAM_OPTION::NOTES[1], PROGRAM_OPTION::NOTES[3] } # ADD ANNOTATIONS 1 & 3 to the BSE System Parameters specification


    # BSE Double Compact Objects
    
    # set the BSE Double Compact Objects specification to MZAMS for Star1, and MZAMS for Star2
    
    BSE_DCO_Rec = { STAR_1_PROPERTY::MZAMS, STAR_2_PROPERTY::MZAMS }
    
    # set the BSE Double Compact Objects specification to empty - nothing will be printed
    # (file will not be created)
    
    BSE_DCO_Rec = {}
    
    
    # BSE Supernovae
    
    BSE_SNE_Rec = {}     # set spec empty - nothing will be printed (file will not be created)
    
    
    # BSE Common Envelopes
    
    BSE_CEE_Rec = {}     # set spec empty - nothing will be printed (file will not be created)
    
    
    # BSE Pulsars
    
    # line ignored (comment). BSE Pulsars specification will be default
    
    # BSE_Pulsars_Rec = { STAR_1_PROPERTY::MASS, STAR_2_PROPERTY::MASS }
    
    
    # BSE Detailed Output
    
    BSE_Detailed_Rec = {} # set spec empty - nothing will be printed (file will not be created)
