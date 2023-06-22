#ifndef __yaml_h__
#define __yaml_h__

#include <algorithm>
#include <chrono>

namespace yaml {

// YAML template rules (in no particular order):
//
//  1. The following records will be automatically written to the start of YAML file:
//         ##~!!~## COMPAS option values
//         ##~!!~## Created at ddd MMM DD HH:MM:SS YYYY by COMPAS vxx.yy.zz
//         ##~!!~## 
//         ##~!!~## The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has
//         ##~!!~## all COMPAS option entries commented so that the COMPAS default value for the
//         ##~!!~## option is used by default. To use a value other than the COMPAS default value,
//         ##~!!~## users must uncomment the entry and change the option value to the desired value.
//
//  2. Lines in the template beginning with "##~!!~##"" will not be preserved (these are assumed to be COMPAS generated headers, and will be rewritten by COMPAS).
//  3. Leading '#' characters on option definition lines in the template will not be preserved (but they may be rewritten by COMPAS).
//  4. Option comments in the template must be preceded by "# " or they will not be preserved.
//  5. Strings in the template beginning with "# Default: " and up to (but not including) the next '#' (or end of line if no #) will not be preserved.
//  6. Strings in the template beginning with "# Options: " and up to (but not including) the next '#' (or end of line if no #) will not be preserved.
//  7. Blank lines in the template will be preserved.
//  8. Option values in the template will not be preserved (but they may be rewritten by COMPAS).
//  9. Option values written by COMPAS will be the option default values unless COMPAS was run with command-line options set - if the user executed COMPAS and
//     specified options on the command line, the user-specified values will be written to the YAML file, and those option records in the YAML file will not be
//     commented.  This gives users the option of creating project-specific YAML files via this method.
// 10. Options present in the template that are not valid COMPAS options will be ignored and not written to the YAML file.
// 11. Any COMPAS options that are not present in the template will be written in alphabetical order at the end of the YAML file.
//
// In the following example template:
//
// 0001     ##~!!~## COMPAS option values
// 0002     ##~!!~## File Created Tue Feb 14 13:09:06 2023 by COMPAS v02.34.06
// 0003     ##~!!~## 
// 0004     ##~!!~## The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has
// 0005     ##~!!~## all COMPAS option entries commented so that the COMPAS default value for the
// 0006     ##~!!~## option is used by default. To use a value other than the COMPAS default value,
// 0007     ##~!!~## users must uncomment the entry and change the option value to the desired value.
// 0008
// 0009     # first comment
// 0010
// 0011     booleanChoices:
// 0012         ### BINARY PROPERTIES
// 0013     #    --allow-touching-at-birth          # Default: False                                        # second comment
// 0014
// 0015         ### STELLAR PROPERTIES
// 0016         --mass-loss-prescription: 'HURLEY'  # Default: 'VINK'  # Options: ['VINK','HURLEY','NONE']    third comment
//
// Lines 0001 - 0007 will not be preserved (but will be replaced by new COMPAS headers).
// The blank line at line 0008 will be preserved.
// The comment "first comment" (on line 0009) will be preserved.
// The blank line at line 0010 will be preserved.
// The header "booleanChoices:" on line 0011 will be preserved.
// The header "### BINARY PROPERTIES" on line 0012 will be preserved.
// The leading '#' on line 0013 will not be preserved (but may be rewritten by COMPAS if the option is set to default).
// The string beginning with "# Default: " and extending to the next '#' on line 0013 will not be preserved (but will be replaced by COMPAS).
// The comment "second comment" on line 0013 will be preserved.
// The blank line at line 0014 will be preserved.
// The header "### STELLAR PROPERTIES" on line 0015 will be preserved.
// The string beginning with "# Default: " and extending to the next '#' on line 0016 will not be preserved (but will be replaced by COMPAS).
// The string beginning with "# Options: " and extending to the next '#' (or, in this case because there is no subsequent #, the end of the 
// line) on line 0016 will not be preserved (but will be replaced by COMPAS).
// The comment "third comment" on line 0016 will not be preserved - there is no "# " prefix, so it will be subsumed by the "# Options: " string
// (which extends from "# Options: " to the end of the line).


// The default COMPAS YAML template follows    
    namespace {
        std::vector<std::string> yamlTemplate {

            "",
            "booleanChoices:",
            "",
            "    ### LOGISTICS",
            "    --debug-to-file",
            "    --detailed-output                                               # WARNING! this creates a data heavy file",
            "    --enable-warnings                                               # option to enable/disable warning messages",
            "    --errors-to-file",
            "    --evolve-unbound-systems",
            "    --population-data-printing",
            "    --print-bool-as-string",
            "    --quiet",
            "    --rlof-printing",
            "    --store-input-files",
            "    --switch-log",
            "",
            "    ### STELLAR PROPERTIES",
            "    --check-photon-tiring-limit",
            "    --use-mass-loss",
            "    --expel-convective-envelope-above-luminosity-threshold",
            "",
            "    ### BINARY PROPERTIES",
            "    --allow-touching-at-birth                                       # record binaries that have stars touching at birth in output files",
            "",
            "    ### MASS TRANSFER",
            "    --angular-momentum-conservation-during-circularisation",
            "    --allow-rlof-at-birth                                           # allow binaries that have one or both stars in RLOF at birth to evolve, particularly useful in the context of CHE binaries",
            "    --circularise-binary-during-mass-transfer",
            "    --hmxr-binaries",
            "    --mass-transfer",
            "    --retain-core-mass-during-caseA-mass-transfer",
            "",
            "    ### COMMON ENVELOPE",
            "    --common-envelope-allow-immediate-RLOF-post-CE-survive",
            "    --common-envelope-allow-main-sequence-survive                   # Allow main sequence stars to survive CE",
            "    --common-envelope-allow-radiative-envelope-survive",
            "    --common-envelope-lambda-nanjing-enhanced",
            "    --common-envelope-lambda-nanjing-interpolate-in-mass",
            "    --common-envelope-lambda-nanjing-interpolate-in-metallicity",
            "    --common-envelope-lambda-nanjing-use-rejuvenated-mass",
            "    --revised-energy-formalism-nandez-ivanova",
            "",
            "    ### SUPERNOVAE, KICKS AND REMNANTS",
            "    --allow-non-stripped-ECSN",
            "    --pair-instability-supernovae",
            "    --pulsational-pair-instability",
            "",
            "    ### PULSAR PARAMETERS",
            "    --evolve-pulsars",
            "",
            "",
            "numericalChoices:",
            "",
            "    ### LOGISTICS",
            "    --debug-level",
            "    --logfile-common-envelopes-record-types",
            "    --logfile-detailed-output-record-types",
            "    --logfile-double-compact-objects-record-types",
            "    --logfile-pulsar-evolution-record-types",
            "    --logfile-rlof-parameters-record-types",
            "    --logfile-supernovae-record-types",
            "    --logfile-system-parameters-record-types",
            "    --grid-lines-to-process",
            "    --grid-start-line",
            "    --hdf5-chunk-size",
            "    --hdf5-buffer-size",
            "    --log-level",
            "    --maximum-evolution-time                                        # maximum physical time a system can be evolved [Myr]",
            "    --maximum-number-timestep-iterations",
            "    --number-of-systems                                             # number of systems per batch",
            "    --timestep-multiplier                                           # optional multiplier relative to default time step duration",
            "",
            "    ### STELLAR PROPERTIES",
            "    --cool-wind-mass-loss-multiplier",
            "    --initial-mass                                                  # initial mass for SSE",
            "    --initial-mass-min                                              # use 5.0 for DCOs [Msol]",
            "    --initial-mass-max                                              # stellar tracks extrapolated above 50 Msol (Hurley+2000) [Msol]",
            "    --initial-mass-power",
            "    --luminosity-to-mass-threshold",
            "    --metallicity                                                   # metallicity for both SSE and BSE - Solar metallicity Asplund+2010",
            "    --metallicity-min",
            "    --metallicity-max",
            "    --luminous-blue-variable-multiplier",
            "    --overall-wind-mass-loss-multiplier",
            "    --random-seed",
            "    --rotational-frequency",
            "    --rotational-frequency-1",
            "    --rotational-frequency-2",
            "    --wolf-rayet-multiplier",
            "",
            "    ### BINARY PROPERTIES",
            "    --eccentricity                                                  # eccentricity for BSE",
            "    --eccentricity-min",
            "    --eccentricity-max",
            "    --initial-mass-1                                                # primary initial mass for BSE",
            "    --initial-mass-2                                                # secondary initial mass for BSE",
            "    --mass-ratio",
            "    --mass-ratio-min",
            "    --mass-ratio-max",
            "    --minimum-secondary-mass                                        # Brown dwarf limit [Msol]",
            "    --orbital-period                                                # orbital period for BSE",
            "    --orbital-period-min                                            # [days]",
            "    --orbital-period-max                                            # [days]",
            "    --semi-major-axis                                               # semi-major axis for BSE",
            "    --semi-major-axis-min                                           # [AU]",
            "    --semi-major-axis-max                                           # [AU]",
            "",
            "    ### MASS TRANSFER",
            "    --convective-envelope-temperature-threshold                     # Only if using envelope-state-prescription = 'FIXED_TEMPERATURE'",
            "    --critical-mass-ratio-HG-degenerate-accretor",
            "    --critical-mass-ratio-HG-non-degenerate-accretor",
            "    --critical-mass-ratio-MS-high-mass-degenerate-accretor",
            "    --critical-mass-ratio-MS-high-mass-non-degenerate-accretor",
            "    --critical-mass-ratio-MS-low-mass-degenerate-accretor",
            "    --critical-mass-ratio-MS-low-mass-non-degenerate-accretor",
            "    --critical-mass-ratio-giant-degenerate-accretor",
            "    --critical-mass-ratio-giant-non-degenerate-accretor",
            "    --critical-mass-ratio-helium-HG-degenerate-accretor",
            "    --critical-mass-ratio-helium-HG-non-degenerate-accretor",
            "    --critical-mass-ratio-helium-MS-degenerate-accretor",
            "    --critical-mass-ratio-helium-giant-degenerate-accretor",
            "    --critical-mass-ratio-helium-MS-non-degenerate-accretor",
            "    --critical-mass-ratio-helium-giant-non-degenerate-accretor",
            "    --critical-mass-ratio-white-dwarf-degenerate-accretor",
            "    --critical-mass-ratio-white-dwarf-non-degenerate-accretor",
            "    --mass-transfer-fa                                              # Only if using mass-transfer-accretion-efficiency-prescription = 'FIXED'",
            "    --mass-transfer-jloss                                           # Only if using mass-transfer-angular-momentum-loss-prescription = 'FIXED'",
            "    --mass-transfer-jloss-macleod-linear-fraction",
            "    --mass-transfer-thermal-limit-C",
            "    --zeta-adiabatic-arbitrary",
            "    --zeta-main-sequence",
            "    --zeta-radiative-envelope-giant",
            "",
            "    ### COMMON ENVELOPE",
            "    --common-envelope-alpha",
            "    --common-envelope-alpha-thermal                                 # lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g",
            "    --common-envelope-lambda                                        # Only if using 'LAMBDA_FIXED'",
            "    --common-envelope-lambda-multiplier                             # Multiply common envelope lambda by some constant",
            "    --common-envelope-mass-accretion-constant",
            "    --common-envelope-mass-accretion-max                            # For 'MACLEOD+2014' [Msol]",
            "    --common-envelope-mass-accretion-min                            # For 'MACLEOD+2014' [Msol]",
            "    --common-envelope-recombination-energy-density",
            "    --common-envelope-slope-kruckow",
            "    --maximum-mass-donor-nandez-ivanova",
            "",
            "    ### SUPERNOVAE, KICKS AND REMNANTS",
            "    --eddington-accretion-factor                                    # multiplication Factor for eddington accretion onto NS&BH",
            "    --fix-dimensionless-kick-magnitude",
            "    --fryer-22-fmix                                                 # parameter describing mixing growth time when using the 'FRYER2022' remnant mass prescription",
            "    --fryer-22-mcrit                                                # critical mass for BH formation when using the 'FRYER2022' remnant mass prescription",
            "    --kick-direction-power",
            "    --kick-magnitude-sigma-CCSN-NS                                  # [km/s]",
            "    --kick-magnitude-sigma-CCSN-BH                                  # [km/s]",
            "    --kick-magnitude-max",
            "    --kick-magnitude-random                                         # (SSE) used to draw the kick magnitude for the star should it undergo a supernova event",
            "    --kick-magnitude                                                # (SSE) (drawn) kick magnitude for the star should it undergo a supernova event [km/s]",
            "    --kick-magnitude-random-1                                       # (BSE) used to draw the kick magnitude for the primary star should it undergo a supernova event",
            "    --kick-magnitude-1                                              # (BSE) (drawn) kick magnitude for the primary star should it undergo a supernova event [km/s]",
            "    --kick-theta-1                                                  # (BSE) angle between the orbital plane and the 'z' axis of the supernova vector for the primary star should it undergo a supernova event [radians]",
            "    --kick-phi-1                                                    # (BSE) angle between 'x' and 'y', both in the orbital plane of the supernova vector, for the primary star should it undergo a supernova event [radians]",
            "    --kick-mean-anomaly-1                                           # (BSE) mean anomaly at the instant of the supernova for the primary star should it undergo a supernova event - should be uniform in [0, 2pi) [radians]",
            "    --kick-magnitude-random-2                                       # (BSE) used to draw the kick velocity for the secondary star should it undergo a supernova event",
            "    --kick-magnitude-2                                              # (BSE) (drawn) kick magnitude for the secondary star should it undergo a supernova event [km/s]",
            "    --kick-theta-2                                                  # (BSE) angle between the orbital plane and the 'z' axis of the supernova vector for the secondary star should it undergo a supernova event [radians]",
            "    --kick-phi-2                                                    # (BSE) angle between 'x' and 'y', both in the orbital plane of the supernova vector, for the secondary star should it undergo a supernova event [radians]",
            "    --kick-mean-anomaly-2                                           # (BSE) mean anomaly at the instant of the supernova for the secondary star should it undergo a supernova event - should be uniform in [0, 2pi) [radians]",
            "    --kick-magnitude-sigma-ECSN                                     # [km/s]",
            "    --kick-magnitude-sigma-USSN                                     # [km/s]",
            "    --kick-scaling-factor",
            "    --maximum-neutron-star-mass",
            "    --mcbur1",
            "    --muller-mandel-kick-multiplier-BH                              # scaling prefactor for BH kicks when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --muller-mandel-kick-multiplier-NS                              # scaling prefactor for NS kicks when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --muller-mandel-sigma-kick                                      # kick scatter when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --neutrino-mass-loss-BH-formation-value",
            "    --pisn-lower-limit                                              # Minimum core mass for PISN [Msol]",
            "    --pisn-upper-limit                                              # Maximum core mass for PISN [Msol]",
            "    --ppi-lower-limit                                               # Minimum core mass for PPI [Msol]",
            "    --ppi-upper-limit                                               # Maximum core mass for PPI [Msol]",
            "",
            "    ### PULSAR PARAMETERS",
            "    --pulsar-birth-magnetic-field-distribution-min                  # [log10(B/G)]",
            "    --pulsar-birth-magnetic-field-distribution-max                  # [log10(B/G)]",
            "    --pulsar-birth-spin-period-distribution-min                     # [ms]",
            "    --pulsar-birth-spin-period-distribution-max                     # [ms]",
            "    --pulsar-magnetic-field-decay-timescale                         # [Myr]",
            "    --pulsar-magnetic-field-decay-massscale                         # [Msol]",
            "    --pulsar-minimum-magnetic-field                                 # [log10(B/G)]",
            "",
            "",
            "stringChoices:",
            "",
            "    ### LOGISTICS",
            "    --add-options-to-sysparms",
            "    --grid                                                          # grid file name (e.g. 'mygrid.txt')",
            "    --mode                                                          # evolving single (SSE) or binary stars (BSE)",
            "    --notes",
            "    --notes-hdrs",
            "    --output-container",
            "",
            "    ### STELLAR PROPERTIES",
            "    --chemically-homogeneous-evolution                              # chemically homogeneous evolution",
            "    --envelope-state-prescription",
            "    --initial-mass-function",
            "    --luminous-blue-variable-prescription",
            "    --mass-loss-prescription",
            "    --metallicity-distribution",
            "    --pulsational-pair-instability-prescription",
            "",
            "    ### BINARY PROPERTIES",
            "    --eccentricity-distribution",
            "    --mass-ratio-distribution",
            "    --orbital-period-distribution",
            "    --rotational-velocity-distribution",
            "    --semi-major-axis-distribution",
            "",
            "    ### MASS TRANSFER",
            "    --case-BB-stability-prescription",
            "    --critical-mass-ratio-prescription",
            "    --stellar-zeta-prescription",
            "    --mass-transfer-angular-momentum-loss-prescription",
            "    --mass-transfer-accretion-efficiency-prescription",
            "    --mass-transfer-rejuvenation-prescription",
            "    --mass-transfer-thermal-limit-accretor",
            "",
            "    ### COMMON ENVELOPE",
            "    --common-envelope-formalism",
            "    --common-envelope-lambda-prescription                           # Xu & Li 2010",
            "    --common-envelope-mass-accretion-prescription",
            "",
            "    ### SUPERNOVAE, KICKS AND REMNANTS",
            "    --black-hole-kicks",
            "    --fryer-supernova-engine",
            "    --kick-magnitude-distribution",
            "    --kick-direction",
            "    --neutron-star-equation-of-state",
            "    --neutrino-mass-loss-BH-formation",
            "    --pulsar-birth-magnetic-field-distribution",
            "    --pulsar-birth-spin-period-distribution",
            "    --remnant-mass-prescription",
            "",
            "    ### LOGFILES AND OUTPUTS ",
            "    --logfile-type",
            "    --logfile-name-prefix",
            "    --logfile-definitions",
            "    --logfile-common-envelopes",
            "    --logfile-detailed-output",
            "    --logfile-double-compact-objects",
            "    --logfile-pulsar-evolution",
            "    --logfile-rlof-parameters",
            "    --logfile-supernovae",
            "    --logfile-switch-log",
            "    --logfile-system-parameters",
            "    --output-path",
            "",
            "",
            "listChoices: ",
            "",
            "    --log-classes",
            "    --debug-classes",
        };
    }


    // namespace functions

    void MakeYAMLfile(const std::string p_YAMLfilename, const std::string p_YAMLtemplate);
    int  ReadYAMLtemplate(const std::string p_YAMLtemplateName);
    int  WriteYAMLfile(const std::string p_YAMLname, const std::vector<std::string> p_YAMLcontent);

}

#endif // __yaml_h__
