#ifndef __Options_H__
#define __Options_H__

#define OPTIONS Options::Instance()

#include <iostream>
#include <string>
#include <sstream>

#include <boost/algorithm/string.hpp>   // Boost string manipulation
#include <boost/program_options.hpp>    // Boost command line options tools
#include <boost/filesystem.hpp>         // Boost filesystem tools for handling paths etc.

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

using std::string;
using std::vector;
using std::get;


// for convenience
#define ANNOUNCE(announceStr)           {                                                                                                           \
                                            std::stringstream _ss; _ss << announceStr; std::cerr << _ss.str() << std::endl;                         \
                                        }

#define WARN(warnStr)                   {                                                                                                           \
                                            std::stringstream _ss; _ss << warnStr; std::cerr << "OPTIONS Warning: " << _ss.str() << std::endl;      \
                                        }

#define COMPLAIN(complainStr)           {                                                                                                           \
                                            std::stringstream _ss; _ss << complainStr; std::cerr << "OPTIONS Error: " << _ss.str() << std::endl;    \
                                            throw("OptionError");                                                                                   \
                                        }

#define COMPLAIN_IF(cond, complainStr)  { if (cond) COMPLAIN(complainStr) }


// JR: todo: one day rename all member variables "m_..."


/*
 * Options Singleton
 *
 * Holds program options
 *
 * Singletons are sometimes frowned-upon, but doing it this way means
 * the program options don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Options {

private:

    Options() {};
    Options(Options const&) = delete;
    Options& operator = (Options const&) = delete;

    static Options* m_Instance;

    string m_OptionsDetails;

    void InitialiseMemberVariables(void);
    COMMANDLINE_STATUS CommandLineSorter(int argc, char * argv[]);

    string ProgramOptionDetails(const boost::program_options::variables_map p_VM);

public:

    static Options* Instance();

    COMMANDLINE_STATUS Initialise(int argc, char *argv[]);
    void SetToFiducialValues(void);

    AIS_DCO                                     AIS_DCOType() const                                                     { return AISDCOtype; }
    string                                      AIS_DCOTypeString() const                                               { return AISDCOtypeString; }
    bool                                        AIS_ExploratoryPhase() const                                            { return AISexploratoryPhase; }
    bool                                        AIS_Hubble() const                                                      { return AIShubble; }
    bool                                        AIS_Pessimistic() const                                                 { return AISpessimistic; }
    bool                                        AIS_RefinementPhase() const                                             { return AISrefinementPhase; }
    bool                                        AIS_RLOF() const                                                        { return AISrlof; }

    bool                                        AllowMainSequenceStarToSurviveCommonEnvelope() const                    { return allowMainSequenceStarToSurviveCommonEnvelope; }
    bool                                        AllowRLOFAtBirth() const                                                { return allowRLOFAtBirth; }
    bool                                        AllowTouchingAtBirth() const                                            { return allowTouchingAtBirth; }
    bool                                        AlwaysStableCaseBBBCFlag() const                                        { return alwaysStableCaseBBBCFlag; }
    bool                                        AngularMomentumConservationDuringCircularisation() const                { return angularMomentumConservationDuringCircularisation; }

    bool                                        BeBinaries() const                                                      { return beBinaries; }

    BLACK_HOLE_KICK_OPTION                      BlackHoleKicksOption() const                                            { return blackHoleKicksOption; }

    CHE_OPTION                                  CHE_Option() const                                                      { return cheOption; }

    bool                                        CirculariseBinaryDuringMassTransfer() const                             { return circulariseBinaryDuringMassTransfer; }

    double                                      CommonEnvelopeAlpha() const                                             { return commonEnvelopeAlpha; }
    double                                      CommonEnvelopeAlphaThermal() const                                      { return commonEnvelopeAlphaThermal; }
    COMMON_ENVELOPE_PRESCRIPTION                CommonEnvelopeHertzsprungGapDonor() const                               { return commonEnvelopeHertzsprungGapDonor; }
    double                                      CommonEnvelopeLambda() const                                            { return commonEnvelopeLambda; }
    double                                      CommonEnvelopeLambdaMultiplier() const                                  { return commonEnvelopeLambdaMultiplier; }
    CE_LAMBDA_PRESCRIPTION                      CommonEnvelopeLambdaPrescription() const                                { return commonEnvelopeLambdaPrescription; }
    double                                      CommonEnvelopeMassAccretionConstant() const                             { return commonEnvelopeMassAccretionConstant; }
    double                                      CommonEnvelopeMassAccretionMax() const                                  { return commonEnvelopeMassAccretionMax; }
    double                                      CommonEnvelopeMassAccretionMin() const                                  { return commonEnvelopeMassAccretionMin; }
    CE_ACCRETION_PRESCRIPTION                   CommonEnvelopeMassAccretionPrescription() const                         { return commonEnvelopeMassAccretionPrescription; }
    double                                      CommonEnvelopeRecombinationEnergyDensity() const                        { return commonEnvelopeRecombinationEnergyDensity; }
    double                                      CommonEnvelopeSlopeKruckow() const                                      { return commonEnvelopeSlopeKruckow; }

    CE_ZETA_PRESCRIPTION                        CommonEnvelopeZetaPrescription() const                                  { return commonEnvelopeZetaPrescription; }

    vector<string>                              DebugClasses() const                                                    { return debugClasses; }
    int                                         DebugLevel() const                                                      { return debugLevel; }
    bool                                        DebugToFile() const                                                     { return debugToFile; }
    bool                                        DetailedOutput() const                                                  { return detailedOutput; }
    bool                                        ErrorsToFile() const                                                    { return errorsToFile; }
    ECCENTRICITY_DISTRIBUTION                   EccentricityDistribution() const                                        { return eccentricityDistribution; }
    double                                      EccentricityDistributionMax() const                                     { return eccentricityDistributionMax; }
    double                                      EccentricityDistributionMin() const                                     { return eccentricityDistributionMin; }
    double                                      EddingtonAccretionFactor() const                                        { return eddingtonAccretionFactor; }
    bool                                        EvolvePulsars() const                                                   { return evolvePulsars; }
    bool                                        EvolveUnboundSystems() const                                            { return evolveUnboundSystems; }

    bool                                        FixedMetallicity() const                                                { return fixedMetallicity; }
    bool                                        FixedRandomSeed() const                                                 { return fixedRandomSeed; }
    double                                      FixedUK() const                                                         { return fixedUK; }                 // JR: todo: this isn't consistent naming see fixedMetallicity, fixedRandomSeed)
    bool                                        ForceCaseBBBCStabilityFlag() const                                      { return forceCaseBBBCStabilityFlag; }
    SN_ENGINE                                   FryerSupernovaEngine() const                                            { return fryerSupernovaEngine; }

    string                                      GridFilename() const                                                    { return gridFilename; }

    INITIAL_MASS_FUNCTION                       InitialMassFunction() const                                             { return initialMassFunction; }
    double                                      InitialMassFunctionMax() const                                          { return initialMassFunctionMax; }
    double                                      InitialMassFunctionMin() const                                          { return initialMassFunctionMin; }
    double                                      InitialMassFunctionPower() const                                        { return initialMassFunctionPower; }

    KICK_DIRECTION_DISTRIBUTION                 KickDirectionDistribution() const                                       { return kickDirectionDistribution; }
    double                                      KickDirectionPower() const                                              { return kickDirectionPower; }
    double                                      KickScalingFactor() const                                               { return kickScalingFactor; }
    KICK_VELOCITY_DISTRIBUTION                  KickVelocityDistribution() const                                        { return kickVelocityDistribution; }

    double                                      KickVelocityDistributionMaximum() const                                 { return kickVelocityDistributionMaximum; }

    double                                      KickVelocityDistributionSigmaCCSN_BH() const                            { return kickVelocityDistributionSigmaCCSN_BH; }
    double                                      KickVelocityDistributionSigmaCCSN_NS() const                            { return kickVelocityDistributionSigmaCCSN_NS; }
    double                                      KickVelocityDistributionSigmaForECSN() const                            { return kickVelocityDistributionSigmaForECSN; }
    double                                      KickVelocityDistributionSigmaForUSSN() const                            { return kickVelocityDistributionSigmaForUSSN; }

    bool                                        LambdaCalculationEveryTimeStep() const                                  { return lambdaCalculationEveryTimeStep; }

    vector<string>                              LogClasses() const                                                      { return logClasses; }
    string                                      LogfileBSEBeBinaries() const                                            { return logfileBSEBeBinaries; }
    string                                      LogfileBSECommonEnvelopes() const                                       { return logfileBSECommonEnvelopes; }
    string                                      LogfileBSEDetailedOutput() const                                        { return logfileBSEDetailedOutput; }
    string                                      LogfileBSEDoubleCompactObjects() const                                  { return logfileBSEDoubleCompactObjects; }
    string                                      LogfileBSEPulsarEvolution() const                                       { return logfileBSEPulsarEvolution; }
    string                                      LogfileBSESupernovae() const                                            { return logfileBSESupernovae; }
    string                                      LogfileBSESystemParameters() const                                      { return logfileBSESystemParameters; }
    string                                      LogfileDefinitionsFilename() const                                      { return logfileDefinitionsFilename; }
    DELIMITER                                   LogfileDelimiter() const                                                { return logfileDelimiter; }
    string                                      LogfileDelimiterString() const                                          { return logfileDelimiterString; }
    string                                      LogfileNamePrefix() const                                               { return logfileNamePrefix; }
    string                                      LogfileSSEParameters() const                                            { return logfileSSEParameters; }
    int                                         LogLevel() const                                                        { return logLevel; }

    double                                      LuminousBlueVariableFactor() const                                      { return luminousBlueVariableFactor; }

    MASS_LOSS_PRESCRIPTION                      MassLossPrescription() const                                            { return massLossPrescription; }

    MASS_RATIO_DISTRIBUTION                     MassRatioDistribution() const                                           { return massRatioDistribution; }
    double                                      MassRatioDistributionMax() const                                        { return massRatioDistributionMax; }
    double                                      MassRatioDistributionMin() const                                        { return massRatioDistributionMin; }

    MT_ACCRETION_EFFICIENCY_PRESCRIPTION        MassTransferAccretionEfficiencyPrescription() const                     { return massTransferAccretionEfficiencyPrescription; }
    MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION       MassTransferAngularMomentumLossPrescription() const                     { return massTransferAngularMomentumLossPrescription; }
    double                                      MassTransferCParameter() const                                          { return massTransferCParameter; }

    bool                                        MassTransferCriticalMassRatioGiant() const                              { return massTransferCriticalMassRatioGiant; }
    double                                      MassTransferCriticalMassRatioGiantDegenerateAccretor() const            { return massTransferCriticalMassRatioGiantDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioGiantNonDegenerateAccretor() const         { return massTransferCriticalMassRatioGiantNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHeliumGiant() const                        { return massTransferCriticalMassRatioHeliumGiant; }
    double                                      MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor() const      { return massTransferCriticalMassRatioHeliumGiantDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor() const   { return massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHeliumHG() const                           { return massTransferCriticalMassRatioHeliumHG; }
    double                                      MassTransferCriticalMassRatioHeliumHGDegenerateAccretor() const         { return massTransferCriticalMassRatioHeliumHGDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor() const      { return massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHeliumMS() const                           { return massTransferCriticalMassRatioHeliumMS; }
    double                                      MassTransferCriticalMassRatioHeliumMSDegenerateAccretor() const         { return massTransferCriticalMassRatioHeliumMSDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor() const      { return massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioHG() const                                 { return massTransferCriticalMassRatioHG; }
    double                                      MassTransferCriticalMassRatioHGDegenerateAccretor() const               { return massTransferCriticalMassRatioHGDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioHGNonDegenerateAccretor() const            { return massTransferCriticalMassRatioHGNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioMSHighMass() const                         { return massTransferCriticalMassRatioMSHighMass; }
    double                                      MassTransferCriticalMassRatioMSHighMassDegenerateAccretor() const       { return massTransferCriticalMassRatioMSHighMassDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor() const    { return massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioMSLowMass() const                          { return massTransferCriticalMassRatioMSLowMass; }
    double                                      MassTransferCriticalMassRatioMSLowMassDegenerateAccretor() const        { return massTransferCriticalMassRatioMSLowMassDegenerateAccretor; }
    double                                      MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor() const     { return massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor; }
    bool                                        MassTransferCriticalMassRatioWhiteDwarf() const                         { return massTransferCriticalMassRatioWhiteDwarf; }

    double                                      MassTransferFractionAccreted() const                                    { return massTransferFractionAccreted; }
    double                                      MassTransferJloss() const                                               { return massTransferJloss; }
    MT_PRESCRIPTION                             MassTransferPrescription() const                                        { return massTransferPrescription; }
    MT_REJUVENATION_PRESCRIPTION                MassTransferRejuvenationPrescription() const                            { return massTransferRejuvenationPrescription; }
    MT_THERMALLY_LIMITED_VARIATION              MassTransferThermallyLimitedVariation() const                           { return massTransferThermallyLimitedVariation; }
    double                                      MaxEvolutionTime() const                                                { return maxEvolutionTime; }
    int                                         MaxNumberOfTimestepIterations() const                                   { return maxNumberOfTimestepIterations; }
    double                                      MaxPercentageAdaptiveMassTransfer() const                               { return maxPercentageAdaptiveMassTransfer; }

    double                                      Metallicity() const                                                     { return metallicity; }

    double                                      MinimumMassSecondary() const                                            { return minimumMassSecondary; }
    double                                      MaximumNeutronStarMass() const                                          { return maximumNeutronStarMass; }

    int                                         NBatchesUsed() const                                                    { return nBatchesUsed; }
    int                                         NBinaries() const                                                       { return nBinaries; }

    NEUTRINO_MASS_LOSS_PRESCRIPTION             NeutrinoMassLossAssumptionBH() const                                    { return neutrinoMassLossAssumptionBH; }
    double                                      NeutrinoMassLossValueBH() const                                         { return neutrinoMassLossValueBH; }

    NS_EOS                                      NeutronStarEquationOfState() const                                      { return neutronStarEquationOfState; }

    bool                                        OptimisticCHE() const                                                   { return cheOption == CHE_OPTION::OPTIMISTIC; }

    string                                      OptionsDetails() const                                                  { return m_OptionsDetails; }

    string                                      OutputContainerName() const                                             { return outputContainerName; }
    string                                      OutputPathString() const                                                { return outputPath.string(); }

    double                                      PairInstabilityLowerLimit() const                                       { return pairInstabilityLowerLimit; }
    double                                      PairInstabilityUpperLimit() const                                       { return pairInstabilityUpperLimit; }

    double                                      PeriodDistributionMax() const                                           { return periodDistributionMax; }
    double                                      PeriodDistributionMin() const                                           { return periodDistributionMin; }

    bool                                        PopulationDataPrinting() const                                          { return populationDataPrinting; }
    bool                                        PrintBoolAsString() const                                               { return printBoolAsString; }

    PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION    PulsarBirthMagneticFieldDistribution() const                            { return pulsarBirthMagneticFieldDistribution; }
    double                                      PulsarBirthMagneticFieldDistributionMax() const                         { return pulsarBirthMagneticFieldDistributionMax; }
    double                                      PulsarBirthMagneticFieldDistributionMin() const                         { return pulsarBirthMagneticFieldDistributionMin; }

    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION       PulsarBirthSpinPeriodDistribution() const                               { return pulsarBirthSpinPeriodDistribution; }
    double                                      PulsarBirthSpinPeriodDistributionMax() const                            { return pulsarBirthSpinPeriodDistributionMax; }
    double                                      PulsarBirthSpinPeriodDistributionMin() const                            { return pulsarBirthSpinPeriodDistributionMin; }

    double                                      PulsarLog10MinimumMagneticField() const                                 { return pulsarLog10MinimumMagneticField; }

    double                                      PulsarMagneticFieldDecayMassscale() const                               { return pulsarMagneticFieldDecayMassscale; }
    double                                      PulsarMagneticFieldDecayTimescale() const                               { return pulsarMagneticFieldDecayTimescale; }

    PPI_PRESCRIPTION                            PulsationalPairInstabilityPrescription() const                          { return pulsationalPairInstabilityPrescription; }
    double                                      PulsationalPairInstabilityLowerLimit() const                            { return pulsationalPairInstabilityLowerLimit; }
    double                                      PulsationalPairInstabilityUpperLimit() const                            { return pulsationalPairInstabilityUpperLimit; }

    bool                                        Quiet() const                                                           { return quiet; }

    unsigned long int                           RandomSeed() const                                                      { return randomSeed; }

    REMNANT_MASS_PRESCRIPTION                   RemnantMassPrescription() const                                         { return remnantMassPrescription; }


    ROTATIONAL_VELOCITY_DISTRIBUTION            RotationalVelocityDistribution() const                                  { return rotationalVelocityDistribution; }

    bool                                        SampleCommonEnvelopeAlpha() const                                       { return sampleCommonEnvelopeAlpha; }
    double                                      SampleCommonEnvelopeAlphaMax() const                                    { return sampleCommonEnvelopeAlphaMax; }
    double                                      SampleCommonEnvelopeAlphaMin() const                                    { return sampleCommonEnvelopeAlphaMin; }

    bool                                        SampleLuminousBlueVariableMultiplier() const                            { return sampleLuminousBlueVariableMultiplier; }
    double                                      SampleLuminousBlueVariableMultiplierMax() const                         { return sampleLuminousBlueVariableMultiplierMax; }
    double                                      SampleLuminousBlueVariableMultiplierMin() const                         { return sampleLuminousBlueVariableMultiplierMin; }

    bool                                        SampleWolfRayetMultiplier() const                                       { return sampleWolfRayetMultiplier; }
    double                                      SampleWolfRayetMultiplierMax() const                                    { return sampleWolfRayetMultiplierMax; }
    double                                      SampleWolfRayetMultiplierMin() const                                    { return sampleWolfRayetMultiplierMin; }

    SEMI_MAJOR_AXIS_DISTRIBUTION                SemiMajorAxisDistribution() const                                       { return semiMajorAxisDistribution; }
    double                                      SemiMajorAxisDistributionMax() const                                    { return semiMajorAxisDistributionMax; }
    double                                      SemiMajorAxisDistributionMin() const                                    { return semiMajorAxisDistributionMin; }
    double                                      SemiMajorAxisDistributionPower() const                                  { return semiMajorAxisDistributionPower; }

    bool                                        SingleStar() const                                                      { return singleStar; }
    int                                         SingleStarMassSteps() const                                             { return singleStarMassSteps; }
    double                                      SingleStarMassMin() const                                               { return singleStarMassMin; }
    double                                      SingleStarMassMax() const                                               { return singleStarMassMax; }

    bool                                        UseFixedUK() const                                                      { return useFixedUK; }
    bool                                        UseMassLoss() const                                                     { return useMassLoss; }
    bool                                        UseMassTransfer() const                                                 { return useMassTransfer; }
    bool                                        UsePairInstabilitySupernovae() const                                    { return usePairInstabilitySupernovae; }
    bool                                        UsePulsationalPairInstability() const                                   { return usePulsationalPairInstability; }

    double                                      WolfRayetFactor() const                                                 { return wolfRayetFactor; }
    bool                                        ZetaCalculationEveryTimeStep() const                                    { return zetaCalculationEveryTimeStep; }
    double                                      ZetaHertzsprungGap() const                                              { return zetaHertzsprungGap; }
    double                                      ZetaMainSequence() const                                                { return zetaMainSequence; }
    double                                      ZetaAdiabaticArbitrary() const                                          { return zetaAdiabaticArbitrary; }
    double                                      ZetaThermalArbitrary() const                                            { return zetaThermalArbitrary; }


    COMPAS_VARIABLE OptionValue(const T_ANY_PROPERTY p_Property) const;

private:

    // member variables - alphabetically in groups (sort of...)

    bool                                        allowRLOFAtBirth;                                               // indicates whether binaries that have one or both stars in RLOF at birth are allowed to evolve
    bool                                        allowTouchingAtBirth;                                           // indicates whether binaries that are touching at birth are allowed to evolve

    bool                                        debugToFile;                                                    // flag used to determine whether debug statements should also be written to a log file
    bool                                        errorsToFile;                                                   // flag used to determine whether error statements should also be written to a log file

    bool                                        singleStar;                                                     // Whether to evolve a single star or a binary

	bool                                        lambdaCalculationEveryTimeStep;							        // Flag indicates if lambda is calculated at each timestep or no
	bool                                        zetaCalculationEveryTimeStep;

	bool                                        beBinaries;													    // Flag if we want to print BeBinaries (main.cpp)
    bool                                        evolvePulsars;                                                  // Whether to evolve pulsars or not
	bool                                        evolveUnboundSystems;							                // Option to chose if unbound systems are evolved until death or the evolution stops after the system is unbound during a SN.

    bool                                        detailedOutput;                                                 // Print detailed output details to file (default = false)
    bool                                        populationDataPrinting;                                         // Print certain data for small populations, but not for larger one
    bool                                        printBoolAsString;                                              // flag used to indicate that boolean properties should be printed as "TRUE" or "FALSE" (default is 1 or 0)
    bool                                        quiet;                                                          // suppress some output

    int                                         nBatchesUsed;                                                   // nr of batches used, only needed for STROOPWAFEL (AIS) (default = -1, not needed)


    // Code variables
    int                                         nBinaries;                                                      // Number of binaries to simulate (default = 10 for quick test)
    bool                                        fixedRandomSeed;                                                // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line)
    unsigned long                               randomSeed;                                                     // Random seed to use
    double                                      maxEvolutionTime;                                               // Maximum time to evolve a binary by
    int                                         maxNumberOfTimestepIterations;                                  // Maximum number of timesteps to evolve binary for before giving up

    // Initial distribution variables

    INITIAL_MASS_FUNCTION                       initialMassFunction;                                            // Which initial mass function to use (default="Kroupa")
    string                                      initialMassFunctionString;                                      // Which initial mass function to use (default="Kroupa")
    double                                      initialMassFunctionMin;                                         // Minimum mass to generate in Msol
    double                                      initialMassFunctionMax;                                         // Maximum mass to generate in Msol
    double                                      initialMassFunctionPower;                                       // single IMF power law set manually

    // Mass ratio
    MASS_RATIO_DISTRIBUTION massRatioDistribution;                                                              // Which mass ratio distribution to use (default = "Flat")
    string                                      massRatioDistributionString;                                    // Which mass ratio distribution to use (default = "Flat")
    double                                      massRatioDistributionMin;                                       // Minimum initial mass ratio when using a distribution
    double                                      massRatioDistributionMax;                                       // Maximum initial mass ratio when using a distribution

    double                                      minimumMassSecondary;                                           // Minimum mass of secondary to draw (in Msol)

    // Semi major axis
    SEMI_MAJOR_AXIS_DISTRIBUTION                semiMajorAxisDistribution;                                      // Which semi-major axis distribution to use (default = "FlatInLog")
    string                                      semiMajorAxisDistributionString;                                // Which semi-major axis distribution to use (default = "FlatInLog")
    double                                      semiMajorAxisDistributionMin;                                   // Minimum a in AU
    double                                      semiMajorAxisDistributionMax;                                   // Maximum a in AU
    double                                      semiMajorAxisDistributionPower;                                 // Set semi-major axis distribution power law slope by hand

    // Period
    double                                      periodDistributionMin;                                          // Minimum initial period in days
    double                                      periodDistributionMax;                                          // Maximum initial period in days

    // Eccentricity
    ECCENTRICITY_DISTRIBUTION                   eccentricityDistribution;                                       // Which eccentricity distribution to use (default = "Zero")
    string                                      eccentricityDistributionString;                                 // Which eccentricity distribution to use (default = "Zero")
    double                                      eccentricityDistributionMin;                                    // Minimum initial eccentricity when using a distribution
    double                                      eccentricityDistributionMax;                                    // Maximum initial eccentricity when using a distribution

    // Supernova variables
    KICK_VELOCITY_DISTRIBUTION                  kickVelocityDistribution;                                       // Which kick velocity distribution to use (default = "Maxwellian". Can also choose "flat")
    string                                      kickVelocityDistributionString;                                 // Which kick velocity distribution to use (default = "Maxwellian". Can also choose "flat")
    double                                      kickVelocityDistributionSigmaCCSN_NS;                           // Kick velocity sigma in km s^-1 for neutron stars (default = "250" )
    double                                      kickVelocityDistributionSigmaCCSN_BH;                           // Kick velocity sigma in km s^-1 for black holes (default = "250" )
    double                                      kickVelocityDistributionMaximum;                                // Maximum kick velocity to draw. If negative, no maximum
	double                                      kickVelocityDistributionSigmaForECSN;			                // Kick velocity sigma for ECSN in km s^-1 (default = "0" )
	double                                      kickVelocityDistributionSigmaForUSSN;			                // Kick velocity sigma for USSN in km s^-1 (default = "20" )
	double                                      kickScalingFactor;								                // Arbitrary factor for scaling kicks

    // Black hole kicks
    BLACK_HOLE_KICK_OPTION                      blackHoleKicksOption;
    string                                      blackHoleKicksString;                                           // Whether to use full black hole kicks (default = "Full")

    // CHE - Chemically Homogeneous Evolution
    CHE_OPTION                                  cheOption;                                                      // whether and how to apply Chemically Homogeneous Evolution
    string                                      cheString;

    // Supernova remnant mass
    REMNANT_MASS_PRESCRIPTION                   remnantMassPrescription;
    string                                      remnantMassPrescriptionString;                                  // Choose which prescription to use to calculate remnant mass
    string                                      fryerSupernovaEngineString;                                     // If using Fryer et al prescription, select which supernova engine to use
    SN_ENGINE                                   fryerSupernovaEngine;                                           // If using Fryer et al prescription, sets whether to use delayed or rapid engine; otherwise ignored (default = SN_DELAYED)

    string                                      neutrinoMassLossAssumptionBHString;                             // String for neutrino mass loss assumption
    NEUTRINO_MASS_LOSS_PRESCRIPTION             neutrinoMassLossAssumptionBH;                                   // Assumption to make about neutrino mass loss for BH formation
    double                                      neutrinoMassLossValueBH;                                        // Value (corresponding to assumption) for neutrino mass loss for BH formation

    // Fixed uk options
    bool                                        useFixedUK;                                                     // Whether to fix uk to a certain value (default is to NOT fix uk)
    double                                      fixedUK;                                                        // Dimensionless value to fix the kick velocity to

    // Kick direction options
    KICK_DIRECTION_DISTRIBUTION                 kickDirectionDistribution;                                      // Which distribution to use for the kick directions
    string                                      kickDirectionDistributionString;                                // Which kick velocity distribution to use (default = "Maxwellian". Can also choose "flat")
    double                                      kickDirectionPower;

    // Pair instability and pulsational pair instability mass loss
    bool                                        usePairInstabilitySupernovae;                                   // Whether to use pair instability supernovae (PISN)
    double                                      pairInstabilityUpperLimit;                                      // Maximum core mass leading to PISN
    double                                      pairInstabilityLowerLimit;                                      // Minimum core mass leading to PISN

    bool                                        usePulsationalPairInstability;                                  // Whether to use pulsational pair instability (PPI)
    double                                      pulsationalPairInstabilityUpperLimit;                           // Minimum core mass leading to PPI
    double                                      pulsationalPairInstabilityLowerLimit;                           // Maximum core mass leading to PPI

    string                                      pulsationalPairInstabilityPrescriptionString;                   // String for which PPI prescription to use
    PPI_PRESCRIPTION                            pulsationalPairInstabilityPrescription;                         // Prescription for PPI to use

	double                                      maximumNeutronStarMass;							                // Maximum mass of a neutron star allowed, set to default in StarTrack

    // Setup default output directory and desired output directory
    string                                      outputPathString;                                               // String to hold the output directory
    boost::filesystem::path                     defaultOutputPath;                                              // Default output location
    boost::filesystem::path                     outputPath;                                                     // Desired output location
    string                                      outputContainerName;                                            // Name of output container (directory)

    // Mass loss options
    bool                                        useMassLoss;                                                    // Whether to activate mass loss (default = True)
    // Can also have options for modifying strength of winds etc here

    MASS_LOSS_PRESCRIPTION                      massLossPrescription;                                           // Which mass loss prescription will be used by the code (default = MASS_LOSS_PRESCRIPTION_NONE)
    string                                      massLossPrescriptionString;                                     // String containing which mass loss prescription to use (default = "None")

    double                                      luminousBlueVariableFactor;                                     // Multiplicitive factor for luminous blue variable (LBV) mass loss rates
    double                                      wolfRayetFactor;                                                // Multiplicitive factor for Wolf-Rayet (WR) wind mass loss rates

    // Mass transfer options
    bool                                        useMassTransfer;                                                // Whether to use mass transfer (default = false)
	bool                                        circulariseBinaryDuringMassTransfer;						    // Whether to circularise binary when it starts (default = false)
	bool                                        forceCaseBBBCStabilityFlag;									    // Whether if all case BB/BC systems are forced to be stable or unstable (default = true)
	bool                                        alwaysStableCaseBBBCFlag;									    // Whether if case BB/BC is always stable (default = true)
	bool                                        angularMomentumConservationDuringCircularisation;			    // Whether to conserve angular momentum while circularising or circularise to periastron (default = false)

    MT_PRESCRIPTION                             massTransferPrescription;                                       // Which mass transfer prescription will be used by the code (default = MASS_TRANSFER_PRESCRIPTION_NONE)
    string                                      massTransferPrescriptionString;                                 // String containing which mass transfer prescription to use (defaul = "None")

    MT_ACCRETION_EFFICIENCY_PRESCRIPTION        massTransferAccretionEfficiencyPrescription;                    // Which mass transfer prescription will be used by the code (default = MASS_TRANSFER_PRESCRIPTION_NONE)
    string                                      massTransferAccretionEfficiencyPrescriptionString;              // String containing which mass transfer prescription to use (defaul = "None")

    double                                      massTransferFractionAccreted;                                   // In mass transfer, ammount of mass transferred that is accreted. 1 for conservative, 0 for fully-non conservative.
    double                                      massTransferCParameter;                                         // Detailed model parameter used in mass transfer
    double                                      eddingtonAccretionFactor;                                       // Multiplication factor for eddington accretion for NS & BH
                                                                                                                // i.e. >1 is super-eddington
                                                                                                                //       0. is no accretion

    double                                      massTransferAdaptiveAlphaParameter;                             // Parameter used in adaptive RLOF to avoid overshoot of the solution
	double                                      maxPercentageAdaptiveMassTransfer;                              // As used in binary_c for adaptive RLOF prescription in mass transfer, according to notes.

	MT_THERMALLY_LIMITED_VARIATION              massTransferThermallyLimitedVariation;                          // Choose how to deal with mass transfer if it is set as thermally limited.
	string                                      massTransferThermallyLimitedVariationString;			        // String


    double                                      massTransferJloss;                                              // Specific angular momentum of the material leaving the system (not accreted)
    MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION       massTransferAngularMomentumLossPrescription;                    // Which mass transfer angular momentum loss prescription will be used by the code (default = MASS_TRANSFER_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_NONE)
    string                                      massTransferAngularMomentumLossPrescriptionString;              // String containing which mass transfer angular momentum loss prescription to use (defaul = "None")

    // Mass transfer rejuvenation prescription
    MT_REJUVENATION_PRESCRIPTION                massTransferRejuvenationPrescription;                           // Which mass transfer rejuvenation prescription (default = MASS_TRANSFER_REJUVENATION_NONE)
    string                                      massTransferRejuvenationPrescriptionString;                     // String containing which mass transfer prescription to use (default = "NONE")

    // Mass transfer critical mass ratios
    bool                                        massTransferCriticalMassRatioMSLowMass;                         // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioMSLowMassNonDegenerateAccretor;    // Critical mass ratio for MT from a MS low mass star
    double                                      massTransferCriticalMassRatioMSLowMassDegenerateAccretor;       // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioMSHighMass;                        // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioMSHighMassNonDegenerateAccretor;   // Critical mass ratio for MT from a MS high mass star
    double                                      massTransferCriticalMassRatioMSHighMassDegenerateAccretor;      // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioHG;                                // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioHGNonDegenerateAccretor;           // Critical mass ratio for MT from a HG star
    double                                      massTransferCriticalMassRatioHGDegenerateAccretor;              // Critical mass ratio for MT from a HG star on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioGiant;                             // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioGiantNonDegenerateAccretor;        // Critical mass ratio for MT from a giant
    double                                      massTransferCriticalMassRatioGiantDegenerateAccretor;           // Critical mass ratio for MT from a giant on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioHeliumMS;                          // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioHeliumMSNonDegenerateAccretor;     // Critical mass ratio for MT from a Helium MS star
    double                                      massTransferCriticalMassRatioHeliumMSDegenerateAccretor;        // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioHeliumHG;                          // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioHeliumHGNonDegenerateAccretor;     // Critical mass ratio for MT from a Helium HG star
    double                                      massTransferCriticalMassRatioHeliumHGDegenerateAccretor;        // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioHeliumGiant;                       // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor;  // Critical mass ratio for MT from a helium giant
    double                                      massTransferCriticalMassRatioHeliumGiantDegenerateAccretor;     // Critical mass ratio for MT from a helium giant on to a degenerate accretor

    bool                                        massTransferCriticalMassRatioWhiteDwarf;                        // Whether to use critical mass ratios
    double                                      massTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor;   // Critical mass ratio for MT from a white dwarf
    double                                      massTransferCriticalMassRatioWhiteDwarfDegenerateAccretor;      // Critical mass ratio for MT from a white dwarf on to a degenerate accretor

    // Common Envelope options
    COMMON_ENVELOPE_PRESCRIPTION                commonEnvelopePrescriptionFlag;                                 // Which common envelope prescription to use
    double                                      commonEnvelopeAlpha;                                            // Common envelope efficiency alpha parameter (default = X)
    double                                      commonEnvelopeLambda;                                           // Common envelope Lambda parameter (default = X)
	double                                      commonEnvelopeSlopeKruckow;									    // Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1
    double                                      commonEnvelopeAlphaThermal;                                     // lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g
    double                                      commonEnvelopeLambdaMultiplier;                                 // Multiply common envelope lambda by some constant
    bool                                        allowMainSequenceStarToSurviveCommonEnvelope;                   // Whether or not to allow a main sequence star to survive a common envelope event

    // Accretion during common envelope
    CE_ACCRETION_PRESCRIPTION                   commonEnvelopeMassAccretionPrescription;
    string                                      commonEnvelopeMassAccretionPrescriptionString;
    double                                      commonEnvelopeMassAccretionMin;
    double                                      commonEnvelopeMassAccretionMax;
    double                                      commonEnvelopeMassAccretionConstant;

    // Common envelope dealing with HG stars
    COMMON_ENVELOPE_PRESCRIPTION                commonEnvelopeHertzsprungGapDonor;                              // Which prescription to use for Hertzsprung gap donors in a CE (default = OPTIMISTIC_HG_CE)
    string                                      commonEnvelopeHertzsprungGapDonorString;                        // String containing which prescription to use for Hertzsprung gap donors in a CE (default = "OPTIMISTIC_HG_CE")

	// Common envelope lambda prescription
	CE_LAMBDA_PRESCRIPTION                      commonEnvelopeLambdaPrescription;								// Which prescription to use for CE lambda (default = LAMBDA_FIXED)
	string                                      commonEnvelopeLambdaPrescriptionString;					        // String containing which prescription to use for CE lambda (default = "LAMBDA_FIXED")


	// Common envelope Nandez and Ivanova energy formalism
	bool                                        revisedEnergyFormalismNandezIvanova	= false;				    // Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	double                                      maximumMassDonorNandezIvanova;								    // Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	double                                      commonEnvelopeRecombinationEnergyDensity;					    // Factor using to calculate the binding energy depending on the mass of the envelope. (default = 1.5x10^13 ergs/g)


    //  Adaptive Importance Sampling options
    bool                                        AISexploratoryPhase;                                            // Flag if we want to run Exploratory phase of Adaptive Importance Sampling // Floor
    AIS_DCO                                     AISDCOtype;                                                     // Which prescription to use for DCO type (default = ALL)
    string                                      AISDCOtypeString;                                               // String containing which type of DCOs to focus on (default = "ALL")
    bool                                        AIShubble;                                                      // whether to exclude DCOs that not merge within Hubble
    bool                                        AISpessimistic;                                                 // whether to exclude Optimistic binaries
    bool                                        AISrefinementPhase;                                             // Flag if we want to run refinement phase of Adaptive Importance Sampling
    bool                                        AISrlof;                                                        // whether to exclude binaries that have RLOFSecondaryZAMS
    double                                      kappaGaussians;                                                 // Scaling factor for the width of the Gaussian distributions in AIS main sampling phase [should be in [0,1]]


    // Which prescription to use for calculating zetas
    CE_ZETA_PRESCRIPTION                        commonEnvelopeZetaPrescription;                                 // Which prescription to use for calculating CE zetas (default = ZETA_ADIABATIC)
    string                                      commonEnvelopeZetaPrescriptionString;                           // String containing which prescription to use for calculating CE zetas (default = STARTRACK)

	double                                      zetaAdiabaticArbitrary;
	double                                      zetaHertzsprungGap;
	double                                      zetaThermalArbitrary;
	double                                      zetaMainSequence;

    // Metallicity options
    bool                                        fixedMetallicity;                                               // Whether user has specified a metallicity to use
    double                                      metallicity;                                                    // Metallicity default = solar


    // Neutron star equation of state
    NS_EOS                                      neutronStarEquationOfState;                                     // which NS EOS to use
    string                                      neutronStarEquationOfStateString;                               // String for which NS EOS to use

    // Pulsar birth magnetic field distribution string
    PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION    pulsarBirthMagneticFieldDistribution;
    string                                      pulsarBirthMagneticFieldDistributionString;                     // Which birth magnetic field distribution to use for pulsars
    double                                      pulsarBirthMagneticFieldDistributionMin;                        // Minimum birth magnetic field (log10 B/G)
    double                                      pulsarBirthMagneticFieldDistributionMax;                        // Maximum birth magnetic field (log10 B/G)

    // Pulsar birth spin period distribution string
    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION       pulsarBirthSpinPeriodDistribution;
    string                                      pulsarBirthSpinPeriodDistributionString;                        // Which birth spin period distribution to use for pulsars
    double                                      pulsarBirthSpinPeriodDistributionMin;                           // Minimum birth spin period (ms)
    double                                      pulsarBirthSpinPeriodDistributionMax;                           // Maximum birth spin period (ms)

    double                                      pulsarMagneticFieldDecayTimescale;                              // Timescale on which magnetic field decays (Myrs)
    double                                      pulsarMagneticFieldDecayMassscale;                              // Mass scale on which magnetic field decays during accretion (solar masses)
    double                                      pulsarLog10MinimumMagneticField;                                // log10 of the minimum pulsar magnetic field in Gauss


    // Rotational Velocity distribution options
    ROTATIONAL_VELOCITY_DISTRIBUTION            rotationalVelocityDistribution;
    string                                      rotationalVelocityDistributionString;                           // Which rotational velocity distribution to use (default = "ZERO")


	bool                                        sampleCommonEnvelopeAlpha;
	double                                      sampleCommonEnvelopeAlphaMax;
	double                                      sampleCommonEnvelopeAlphaMin;

	bool                                        sampleKickDirectionPower;
	double                                      sampleKickDirectionPowerMax;
	double                                      sampleKickDirectionPowerMin;

	bool                                        sampleKickVelocitySigma;
	double                                      sampleKickVelocitySigmaMax;
	double                                      sampleKickVelocitySigmaMin;

	bool                                        sampleLuminousBlueVariableMultiplier;
	double                                      sampleLuminousBlueVariableMultiplierMax;
	double                                      sampleLuminousBlueVariableMultiplierMin;

	bool                                        sampleWolfRayetMultiplier;
	double                                      sampleWolfRayetMultiplierMax;
	double                                      sampleWolfRayetMultiplierMin;


	// grids

    string                                      gridFilename;                                                   // grid filename - SSE or BSE


    // debug and logging options

    int                                         debugLevel;                                                     // debug level - used to determine which debug statements are actually written
    vector<string>                              debugClasses;                                                   // debug classes - used to determine which debug statements are actually written

    int                                         logLevel;                                                       // logging level - used to determine which logging statements are actually written
    vector<string>                              logClasses;                                                     // logging classes - used to determine which logging statements are actually written

    string                                      logfileNamePrefix;                                              // prefix for log file names

    string                                      logfileDefinitionsFilename;                                     // filename for the logfile record definitions

    string                                      logfileDelimiterString;                                         // field delimiter for log file records (program option string)
    DELIMITER                                   logfileDelimiter;                                               // field delimiter for log file records

    // SSE options
    string                                      logfileSSEParameters;                                           // SSE output file name

    int                                         singleStarMassSteps;                                            // Number of stars of different masses to evolve
    double                                      singleStarMassMin;                                              // The minimum mass to use for SSE (i.e. the mass of the first star to be evolved)
    double                                      singleStarMassMax;                                              // The maximum mass to use for SSE (i.e. the mass of the last star to be evolved)


    // BSE options
    string                                      logfileBSESystemParameters;                                     // BSE output file name: system parameters
    string                                      logfileBSEDetailedOutput;                                       // BSE output file name: detailed output
    string                                      logfileBSEDoubleCompactObjects;                                 // BSE output file name: double compact objects
    string                                      logfileBSESupernovae;                                           // BSE output file name: supernovae
    string                                      logfileBSECommonEnvelopes;                                      // BSE output file name: common envelopes
    string                                      logfileBSEBeBinaries;                                           // BSE output file name: Be Binaries
    string                                      logfileBSEPulsarEvolution;                                      // BSE output file name: pulsar evolution

};


#endif // __Options_H__
