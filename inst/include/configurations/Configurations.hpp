
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
* @file Configurations.hpp
* @version 1.1.0
* @brief Contains the declaration of the Configurations class and its member functions.
* @author Mahmoud ElKarargy
* @author Sameh Abdulah
* @date 2024-02-04
**/

#ifndef EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_CONFIGURATIONS_HPP

#include <vector>
#include <unordered_map>
#include <any>

#include <common/Definitions.hpp>

/**
 * @brief Macro that generates a setter function for a member variable.
 * @details This macro generates a function named Set##name that takes an argument of
 * the specified type and sets the member variable with the specified name
 * to the value of the argument. The name of the member variable is used as
 * the key to set the corresponding value in the specified dictionary.
 * @param[in] name The name of the member variable to be set.
 * @param[in] type The data type of the member variable.
 * @param[in] argument_name The name of the argument to the generated function.
 * @param[in] dictionary_name The name of the dictionary to set the value in.
 *
 */
#define CREATE_SETTER_FUNCTION(name, type, argument_name, dictionary_name)  \
void Set##name(type argument_name)                                          \
{                                                                           \
    mDictionary[dictionary_name] = argument_name;                           \
}

/**
 * @brief Macro that generates a getter function for a member variable.
 * @details This macro generates a function named Get##name that returns the value of
 * the member variable with the specified name from the specified dictionary.
 * @param[in] name The name of the member variable to be retrieved.
 * @param[in] type The data type of the member variable.
 * @param[in] dictionary_name The name of the dictionary to retrieve the value from.
 *
 */
#define CREATE_GETTER_FUNCTION(name, type, dictionary_name)                                                 \
type Get##name()                                                                                            \
{                                                                                                           \
    if (mDictionary.find(dictionary_name) == mDictionary.end()) {                                           \
        throw std::range_error(std::string("Argument ").append(dictionary_name).append(" is not set!"));    \
    }                                                                                                       \
    return std::any_cast<type>(mDictionary[dictionary_name]);                                               \
}

namespace exageostat::configurations {
    /**
     * @class Configurations
     * @brief Contains methods to set and get.
     *
     */
    class Configurations {
    public:

        /**
         * @brief Constructor initializing a Configuration object with default values.
         *
         */
        Configurations();

        /**
         * @brief destructor to allow calls to the correct concrete destructor.
         *
         */
        ~Configurations();

        /**
         * @brief Initialize the module arguments.
         * @param[in] aArgC The number of arguments being passed into the program from the command line.
         * @param[in] apArgV The array of arguments.
         * @param[in] aEnableR check if R is enabled
         * @details This method initializes the command line arguments and set default values for unused args.
         * @return void
         *
         */
        void InitializeArguments(const int &aArgC, char **apArgV, const bool &aEnableR = false);

        /**
         * @brief Initialize the all theta arguments.
         * @return void
         *
         */
        void InitializeAllTheta();
        /**
         * @brief Print the usage and accepted Arguments.
         * @return void
         *
         */
        static void PrintUsage();

        /**
         * @brief Validate the config through a set of if/else.
         * @throw exception in case some if/else conditions are not met.
         * @return void
         *
         */
        void ValidateConfiguration();

        /** START OF THE COMMON ARGUMENTS BETWEEN ALL MODULES. **/

        CREATE_SETTER_FUNCTION(ProblemSize, int, aProblemSize, "n")
        CREATE_GETTER_FUNCTION(ProblemSize, int, "n")

        CREATE_SETTER_FUNCTION(KernelName, const std::string&, aKernel, "kernel")
        CREATE_GETTER_FUNCTION(KernelName, const std::string&, "kernel")

        CREATE_SETTER_FUNCTION(PGrid, int, aPGrid, "p")
        CREATE_GETTER_FUNCTION(PGrid, int, "p")

        CREATE_SETTER_FUNCTION(QGrid, int, aQGrid, "q")
        CREATE_GETTER_FUNCTION(QGrid, int, "q")

        CREATE_SETTER_FUNCTION(TimeSlot, int, aTimeSlot, "timeslot")
        CREATE_GETTER_FUNCTION(TimeSlot, int, "timeslot")

        CREATE_SETTER_FUNCTION(Computation, common::Computation, aComputation, "computation")
        CREATE_GETTER_FUNCTION(Computation, common::Computation, "computation")

        CREATE_SETTER_FUNCTION(Precision, common::Precision, aPrecision, "precision")
        CREATE_GETTER_FUNCTION(Precision, common::Precision, "precision")

        CREATE_SETTER_FUNCTION(CoresNumber, int, aCoresNumbers, "cores")
        CREATE_GETTER_FUNCTION(CoresNumber, int, "cores")

        CREATE_SETTER_FUNCTION(GPUsNumbers, int, aGPUsNumber, "gpus")
        CREATE_GETTER_FUNCTION(GPUsNumbers, int, "gpus")

        CREATE_SETTER_FUNCTION(DenseTileSize, int, aTileSize, "dts")
        CREATE_GETTER_FUNCTION(DenseTileSize, int, "dts")

        CREATE_SETTER_FUNCTION(LowTileSize, int, aTileSize, "lts")
        CREATE_GETTER_FUNCTION(LowTileSize, int, "lts")

        CREATE_SETTER_FUNCTION(Band, int, aBand, "band")
        CREATE_GETTER_FUNCTION(Band, int, "band")

        CREATE_SETTER_FUNCTION(MaxRank, int, aMaxRank, "maxrank")
        CREATE_GETTER_FUNCTION(MaxRank, int, "maxrank")

        CREATE_SETTER_FUNCTION(ActualObservationsFilePath, const std::string&, aActualObservationsFilePath, "observationsfile")
        CREATE_GETTER_FUNCTION(ActualObservationsFilePath, std::string, "observationsfile")

        CREATE_SETTER_FUNCTION(Seed, int, aSeed, "seed")
        CREATE_GETTER_FUNCTION(Seed, int, "seed")

        CREATE_SETTER_FUNCTION(LoggerPath, const std::string&, aLoggerPath, "logpath")
        CREATE_GETTER_FUNCTION(LoggerPath, std::string, "logpath")

        CREATE_SETTER_FUNCTION(InitialTheta, const std::vector<double>&, apTheta, "initialtheta")
        CREATE_GETTER_FUNCTION(InitialTheta, std::vector<double>&, "initialtheta")

        CREATE_SETTER_FUNCTION(IsOOC, bool, aIsOOC, "ooc")
        CREATE_GETTER_FUNCTION(IsOOC, bool, "ooc")

        CREATE_SETTER_FUNCTION(ApproximationMode, int, aApproximationMode, "approximationmode")
        CREATE_GETTER_FUNCTION(ApproximationMode, int, "approximationmode")

        CREATE_SETTER_FUNCTION(Logger, bool, aLogger, "log")
        CREATE_GETTER_FUNCTION(Logger, bool, "log")

        CREATE_SETTER_FUNCTION(LowerBounds, const std::vector<double>&, apTheta, "lb")
        CREATE_GETTER_FUNCTION(LowerBounds, std::vector<double>&, "lb")

        CREATE_SETTER_FUNCTION(UpperBounds, const std::vector<double>&, apTheta, "ub")
        CREATE_GETTER_FUNCTION(UpperBounds, std::vector<double>&, "ub")

        CREATE_SETTER_FUNCTION(EstimatedTheta, const std::vector<double>&, apTheta, "estimatedtheta")
        CREATE_GETTER_FUNCTION(EstimatedTheta, std::vector<double>&, "estimatedtheta")

        CREATE_SETTER_FUNCTION(StartingTheta, const std::vector<double>&, apTheta, "startingtheta")
        CREATE_GETTER_FUNCTION(StartingTheta, std::vector<double>&, "startingtheta")

        CREATE_SETTER_FUNCTION(IsNonGaussian, bool, aIsNonGaussian, "isnongaussian")
        CREATE_GETTER_FUNCTION(IsNonGaussian, bool, "isnongaussian")

        /**
         * @brief Getter for the verbosity.
         * @return The verbosity mode.
         *
         */
        static exageostat::common::Verbose GetVerbosity();

        static void SetVerbosity(const common::Verbose &aVerbose);

        /** END OF THE COMMON ARGUMENTS BETWEEN ALL MODULES. **/

        /** START OF THE HICMA-PARSEC SPECIFIC ARGUEMNTS. **/

        CREATE_SETTER_FUNCTION(DenseBandDP, int, aDenseBandDP, "banddense")
        CREATE_GETTER_FUNCTION(DenseBandDP, int, "banddense")

        CREATE_SETTER_FUNCTION(ObjectsNumber, int, aObjectsNumber, "objectsnumber")
        CREATE_GETTER_FUNCTION(ObjectsNumber, int, "objectsnumber")

        CREATE_SETTER_FUNCTION(AdaptiveDecision, int, aAdaptiveDecision, "adaptivedecision")
        CREATE_GETTER_FUNCTION(AdaptiveDecision, int, "adaptivedecision")

        CREATE_SETTER_FUNCTION(DiagonalAddition, int, aDiagonalAddition, "adddiagonal")
        CREATE_GETTER_FUNCTION(DiagonalAddition, int, "adddiagonal")

        CREATE_SETTER_FUNCTION(TimeSlotPerFile, int, aTimeSlotPerFile, "filetimeslot")
        CREATE_GETTER_FUNCTION(TimeSlotPerFile, int, "filetimeslot")

        CREATE_SETTER_FUNCTION(FileNumber, int, aFileNumber, "filenumber")
        CREATE_GETTER_FUNCTION(FileNumber, int, "filenumber")

        CREATE_SETTER_FUNCTION(EnableInverse, bool, aEnableInverse, "enableinverse")
        CREATE_GETTER_FUNCTION(EnableInverse, bool, "enableinverse")

        CREATE_SETTER_FUNCTION(MPIIO, bool, aMPIIO, "mpiio")
        CREATE_GETTER_FUNCTION(MPIIO, bool, "mpiio")

        CREATE_SETTER_FUNCTION(StageZero, bool, aIsEnabled, "stagezero")
        CREATE_GETTER_FUNCTION(StageZero, bool, "stagezero")

        CREATE_SETTER_FUNCTION(ForcingDataPath, const std::string&, aPath, "forcingdatapath")
        CREATE_GETTER_FUNCTION(ForcingDataPath, std::string, "forcingdatapath")

        CREATE_SETTER_FUNCTION(NetCDFDataPath, const std::string&, aPath, "netcdfdatapath")
        CREATE_GETTER_FUNCTION(NetCDFDataPath, std::string, "netcdfdatapath")

        CREATE_SETTER_FUNCTION(StartYear, int, aStartYear, "startyear")
        CREATE_GETTER_FUNCTION(StartYear, int, "startyear")

        CREATE_SETTER_FUNCTION(EndYear, int, aEndYear, "endyear")
        CREATE_GETTER_FUNCTION(EndYear, int, "endyear")

/** END OF THE HICMA-PARSEC SPECIFIC ARGUMENTS. **/
/** START OF THE DATA GENERATION MODULES. **/

        CREATE_SETTER_FUNCTION(Dimension, exageostat::common::Dimension, aDimension, "dimension")
        CREATE_GETTER_FUNCTION(Dimension, exageostat::common::Dimension, "dimension")

        CREATE_SETTER_FUNCTION(IsSynthetic, bool, aIsSynthetic, "issynthetic")
        CREATE_GETTER_FUNCTION(IsSynthetic, bool, "issynthetic")

        CREATE_SETTER_FUNCTION(DataPath, const std::string&, aDataPath, "datapath")
        CREATE_GETTER_FUNCTION(DataPath, std::string, "datapath")

        // Results output directory path (for CSVs and parameters)
        CREATE_SETTER_FUNCTION(ResultsPath, const std::string&, aResultsPath, "resultspath")
        CREATE_GETTER_FUNCTION(ResultsPath, std::string, "resultspath")

        // StageZero climate data grid parameters (required for Stage Zero)
        CREATE_SETTER_FUNCTION(LatitudeBand, int, aLatBand, "lat")
        CREATE_GETTER_FUNCTION(LatitudeBand, int, "lat")
        
        CREATE_SETTER_FUNCTION(LongitudeCount, int, aLonCount, "lon")
        CREATE_GETTER_FUNCTION(LongitudeCount, int, "lon")

/** END OF THE DATA GENERATION MODULES. **/
/** START OF THE DATA MODELING MODULES. **/

        CREATE_SETTER_FUNCTION(RecoveryFile, const std::string&, aRecoveryFile, "recoveryfile")
        CREATE_GETTER_FUNCTION(RecoveryFile, std::string, "recoveryfile")

        CREATE_SETTER_FUNCTION(FileLogPath, FILE *, apFileLogPath, "filelogpath")
        CREATE_GETTER_FUNCTION(FileLogPath, FILE *, "filelogpath")

        CREATE_SETTER_FUNCTION(FileLogName, const std::string&, aFileLogName, "filelogname")

        CREATE_SETTER_FUNCTION(DistanceMetric, common::DistanceMetric, aDistanceMetric, "distancemetric")
        CREATE_GETTER_FUNCTION(DistanceMetric, common::DistanceMetric, "distancemetric")

        CREATE_SETTER_FUNCTION(MaxMleIterations, int, aMaxMleIterations, "maxmleiterations")
        CREATE_GETTER_FUNCTION(MaxMleIterations, int, "maxmleiterations")

        CREATE_SETTER_FUNCTION(Accuracy, int, aAccuracy, "accuracy")
        CREATE_GETTER_FUNCTION(Accuracy, int, "accuracy")

        void SetTolerance(double aTolerance);
        CREATE_GETTER_FUNCTION(Tolerance, double, "tolerance")

        /** END OF THE DATA MODELING MODULES. **/
        /** START OF THE DATA PREDICTION MODULES. **/

        CREATE_SETTER_FUNCTION(UnknownObservationsNb, int, aUnknownObservationsNumber, "zmiss")
        CREATE_GETTER_FUNCTION(UnknownObservationsNb, int, "zmiss")

        CREATE_SETTER_FUNCTION(IsMSPE, bool, aIsMSPE, "mspe")
        CREATE_GETTER_FUNCTION(IsMSPE, bool, "mspe")

        CREATE_SETTER_FUNCTION(IsIDW, bool, aIsIDW, "idw")
        CREATE_GETTER_FUNCTION(IsIDW, bool, "idw")

        CREATE_SETTER_FUNCTION(IsMLOEMMOM, bool, aIsMLOEMMOM, "mloemmom")
        CREATE_GETTER_FUNCTION(IsMLOEMMOM, bool, "mloemmom")

        CREATE_SETTER_FUNCTION(IsFisher, bool, aIsFisher, "fisher")
        CREATE_GETTER_FUNCTION(IsFisher, bool, "fisher")

        CREATE_SETTER_FUNCTION(ObservationNumber, int, aObservationsNumber, "observationnumber")
        CREATE_GETTER_FUNCTION(ObservationNumber, int, "observationnumber")

        /** END OF THE DATA PREDICTION MODULES. **/

        /**
         * @brief Initialize a vector with a given size to contain zeros.
         * @param[in, out] aTheta A reference to the vector to initialize.
         * @param[in] aSize The size of the vector to initialize.
         * @return void.
         *
         */
        static void InitTheta(std::vector<double> &aTheta, const int &aSize);

        /**
         * @brief print the summary of MLE inputs.
         * @return void
         *
         */
        void PrintSummary();

        /**
         * @brief Calculates the number of observed measurements.
         * @return number of observed measurements.
         *
         */
        int CalculateZObsNumber();

    private:
        /// Used Dictionary
        std::unordered_map<std::string, std::any> mDictionary;
        /// Used Argument counter
        int mArgC = 0;
        /// Used Argument vectors
        char **mpArgV = nullptr;
        //// Used run mode
        static exageostat::common::Verbose mVerbosity;
        //// Used bool for init theta
        static bool mIsThetaInit;
        //// Used bool for R allocated memory on heap
        static bool mHeapAllocated;
        //// Used bool for indicating the first init of configurations for printing summary
        static bool mFirstInit;
    };
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
