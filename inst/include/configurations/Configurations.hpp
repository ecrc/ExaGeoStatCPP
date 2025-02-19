
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
         * @brief Initialize data generation arguments..
         * @return void
         *
         */
        void InitializeDataGenerationArguments();

        /**
         * @brief Initialize data Modeling arguments.
         * @return void
         *
         */
        void InitializeDataModelingArguments();

        /**
         * @brief Initialize data Prediction arguments.
         * @return void
         *
         */
        void InitializeDataPredictionArguments();

        /**
         * @brief Print the usage and accepted Arguments.
         * @return void
         *
         */
        static void PrintUsage();

        /** START OF THE COMMON ARGUMENTS BETWEEN ALL MODULES. **/

        CREATE_SETTER_FUNCTION(ProblemSize, int, aProblemSize, "ProblemSize")

        CREATE_GETTER_FUNCTION(ProblemSize, int, "ProblemSize")

        CREATE_SETTER_FUNCTION(KernelName, const std::string&, aKernel, "Kernel")

        CREATE_GETTER_FUNCTION(KernelName, const std::string&, "Kernel")

        CREATE_SETTER_FUNCTION(PGrid, int, aPGrid, "PGrid")

        CREATE_GETTER_FUNCTION(PGrid, int, "PGrid")

        CREATE_SETTER_FUNCTION(QGrid, int, aQGrid, "QGrid")

        CREATE_GETTER_FUNCTION(QGrid, int, "QGrid")

        CREATE_SETTER_FUNCTION(TimeSlot, int, aTimeSlot, "TimeSlot")

        CREATE_GETTER_FUNCTION(TimeSlot, int, "TimeSlot")

        CREATE_SETTER_FUNCTION(Computation, common::Computation, aComputation, "Computation")

        CREATE_GETTER_FUNCTION(Computation, common::Computation, "Computation")

        CREATE_SETTER_FUNCTION(Precision, common::Precision, aPrecision, "Precision")

        CREATE_GETTER_FUNCTION(Precision, common::Precision, "Precision")

        CREATE_SETTER_FUNCTION(CoresNumber, int, aCoresNumbers, "CoresNumbers")

        CREATE_GETTER_FUNCTION(CoresNumber, int, "CoresNumbers")

        CREATE_SETTER_FUNCTION(GPUsNumbers, int, aGPUsNumber, "GPUsNumbers")

        CREATE_GETTER_FUNCTION(GPUsNumbers, int, "GPUsNumbers")

        CREATE_SETTER_FUNCTION(DenseTileSize, int, aTileSize, "DTS")

        CREATE_GETTER_FUNCTION(DenseTileSize, int, "DTS")

        CREATE_SETTER_FUNCTION(LowTileSize, int, aTileSize, "LTS")

        CREATE_GETTER_FUNCTION(LowTileSize, int, "LTS")

        CREATE_SETTER_FUNCTION(Band, int, aBand, "Band")

        CREATE_GETTER_FUNCTION(Band, int, "Band")

        CREATE_SETTER_FUNCTION(MaxRank, int, aMaxRank, "MaxRank")

        CREATE_GETTER_FUNCTION(MaxRank, int, "MaxRank")

        CREATE_SETTER_FUNCTION(ActualObservationsFilePath, const std::string &, aActualObservationsFilePath,
                               "ActualObservationsFilePath")

        CREATE_GETTER_FUNCTION(ActualObservationsFilePath, std::string, "ActualObservationsFilePath")

        CREATE_SETTER_FUNCTION(Seed, int, aSeed, "Seed")

        CREATE_GETTER_FUNCTION(Seed, int, "Seed")

        CREATE_SETTER_FUNCTION(LoggerPath, const std::string&, aLoggerPath, "LoggerPath")

        CREATE_GETTER_FUNCTION(LoggerPath, std::string, "LoggerPath")

        CREATE_SETTER_FUNCTION(InitialTheta, const std::vector<double> &, apTheta, "InitialTheta")

        CREATE_GETTER_FUNCTION(InitialTheta, std::vector<double> &, "InitialTheta")

        CREATE_SETTER_FUNCTION(IsOOC, bool, aIsOOC, "OOC")

        CREATE_GETTER_FUNCTION(IsOOC, bool, "OOC")

        CREATE_SETTER_FUNCTION(ApproximationMode, int, aApproximationMode, "ApproximationMode")

        CREATE_GETTER_FUNCTION(ApproximationMode, int, "ApproximationMode")

        CREATE_SETTER_FUNCTION(Logger, bool, aLogger, "Logger")

        CREATE_GETTER_FUNCTION(Logger, bool, "Logger")

        CREATE_SETTER_FUNCTION(LowerBounds, const std::vector<double> &, apTheta, "LowerBounds")

        CREATE_GETTER_FUNCTION(LowerBounds, std::vector<double> &, "LowerBounds")

        CREATE_SETTER_FUNCTION(UpperBounds, const std::vector<double> &, apTheta, "UpperBounds")

        CREATE_GETTER_FUNCTION(UpperBounds, std::vector<double> &, "UpperBounds")

        CREATE_SETTER_FUNCTION(EstimatedTheta, const std::vector<double> &, apTheta, "EstimatedTheta")

        CREATE_GETTER_FUNCTION(EstimatedTheta, std::vector<double> &, "EstimatedTheta")

        CREATE_SETTER_FUNCTION(StartingTheta, const std::vector<double> &, apTheta, "StartingTheta")

        CREATE_GETTER_FUNCTION(StartingTheta, std::vector<double> &, "StartingTheta")

        CREATE_SETTER_FUNCTION(IsNonGaussian, bool, aIsNonGaussian, "IsNonGaussian")

        CREATE_GETTER_FUNCTION(IsNonGaussian, bool, "IsNonGaussian")

        /**
         * @brief Getter for the verbosity.
         * @return The verbosity mode.
         *
         */
        static exageostat::common::Verbose GetVerbosity();

        static void SetVerbosity(const common::Verbose &aVerbose);

        /** END OF THE COMMON ARGUMENTS BETWEEN ALL MODULES. **/

        /** START OF THE HICMA-PARSEC SPECIFIC ARGUEMNTS. **/

        CREATE_SETTER_FUNCTION(DenseBandDP, int, aDenseBandDP, "DenseBandDoublePrecision")

        CREATE_GETTER_FUNCTION(DenseBandDP, int, "DenseBandDoublePrecision")

        CREATE_SETTER_FUNCTION(ObjectsNumber, int, aObjectsNumber, "ObjectsNumber")

        CREATE_GETTER_FUNCTION(ObjectsNumber, int, "ObjectsNumber")

        CREATE_SETTER_FUNCTION(AdaptiveDecision, int, aAdaptiveDecision, "AdaptiveDecision")

        CREATE_GETTER_FUNCTION(AdaptiveDecision, int, "AdaptiveDecision")

        CREATE_SETTER_FUNCTION(DiagonalAddition, int, aDiagonalAddition, "DiagonalAddition")

        CREATE_GETTER_FUNCTION(DiagonalAddition, int, "DiagonalAddition")

        CREATE_SETTER_FUNCTION(TimeSlotPerFile, int, aTimeSlotPerFile, "TimeSlotPerFile")

        CREATE_GETTER_FUNCTION(TimeSlotPerFile, int, "TimeSlotPerFile")

        CREATE_SETTER_FUNCTION(FileNumber, int, aFileNumber, "FileNumber")

        CREATE_GETTER_FUNCTION(FileNumber, int, "FileNumber")

        CREATE_SETTER_FUNCTION(EnableInverse, bool , aEnableInverse, "EnableInverse")

        CREATE_GETTER_FUNCTION(EnableInverse, bool , "EnableInverse")

        CREATE_SETTER_FUNCTION(MPIIO, bool, aMPIIO, "MPIIO")

        CREATE_GETTER_FUNCTION(MPIIO, bool, "MPIIO")

        /** END OF THE HICMA-PARSEC SPECIFIC ARGUEMNTS. **/
        /** START OF THE DATA GENERATION MODULES. **/

        CREATE_SETTER_FUNCTION(Dimension, exageostat::common::Dimension, aDimension, "Dimension")

        CREATE_GETTER_FUNCTION(Dimension, exageostat::common::Dimension, "Dimension")

        CREATE_SETTER_FUNCTION(IsSynthetic, bool, aIsSynthetic, "IsSynthetic")

        CREATE_GETTER_FUNCTION(IsSynthetic, bool, "IsSynthetic")

        CREATE_SETTER_FUNCTION(DataPath, const std::string&, aDataPath, "DataPath")

        CREATE_GETTER_FUNCTION(DataPath, std::string, "DataPath")

        /** END OF THE DATA GENERATION MODULES. **/
        /** START OF THE DATA MODELING MODULES. **/

        CREATE_SETTER_FUNCTION(RecoveryFile, const std::string&, aRecoveryFile, "RecoveryFile")

        CREATE_GETTER_FUNCTION(RecoveryFile, std::string, "RecoveryFile")

        CREATE_SETTER_FUNCTION(FileLogPath, FILE *, apFileLogPath, "FileLogPath")

        CREATE_GETTER_FUNCTION(FileLogPath, FILE *, "FileLogPath")

        CREATE_SETTER_FUNCTION(FileLogName, const std::string&, aFileLogName, "FileLogName")

        CREATE_SETTER_FUNCTION(DistanceMetric, common::DistanceMetric, aDistanceMetric, "DistanceMetric")

        CREATE_GETTER_FUNCTION(DistanceMetric, common::DistanceMetric, "DistanceMetric")

        CREATE_SETTER_FUNCTION(MaxMleIterations, int, aMaxMleIterations, "MaxMleIterations")

        CREATE_GETTER_FUNCTION(MaxMleIterations, int, "MaxMleIterations")

        CREATE_SETTER_FUNCTION(Accuracy, int, aAccuracy, "Accuracy")

        CREATE_GETTER_FUNCTION(Accuracy, int, "Accuracy")

        void SetTolerance(double aTolerance);

        CREATE_GETTER_FUNCTION(Tolerance, double, "Tolerance")

        /** END OF THE DATA MODELING MODULES. **/
        /** START OF THE DATA PREDICTION MODULES. **/

        CREATE_SETTER_FUNCTION(UnknownObservationsNb, int, aUnknownObservationsNumber, "UnknownObservationsNb")

        CREATE_GETTER_FUNCTION(UnknownObservationsNb, int, "UnknownObservationsNb")

        CREATE_SETTER_FUNCTION(IsMSPE, bool, aIsMSPE, "IsMSPE")

        CREATE_GETTER_FUNCTION(IsMSPE, bool, "IsMSPE")

        CREATE_SETTER_FUNCTION(IsIDW, bool, aIsIDW, "IsIDW")

        CREATE_GETTER_FUNCTION(IsIDW, bool, "IsIDW")

        CREATE_SETTER_FUNCTION(IsMLOEMMOM, bool, aIsMLOEMMOM, "IsMLOEMMOM")

        CREATE_GETTER_FUNCTION(IsMLOEMMOM, bool, "IsMLOEMMOM")

        CREATE_SETTER_FUNCTION(IsFisher, bool, aIsFisher, "IsFisher")

        CREATE_GETTER_FUNCTION(IsFisher, bool, "IsFisher")

        CREATE_SETTER_FUNCTION(ObservationNumber, int, aObservationsNumber, "ObservationNumber")

        CREATE_GETTER_FUNCTION(ObservationNumber, int, "ObservationNumber")

        /** END OF THE DATA PREDICTION MODULES. **/

        /**
         * @brief Check if input value is numerical.
         * @param[in] aValue The input from the user side.
         * @return The int casted value.
         *
         */
        static int CheckNumericalValue(const std::string &aValue);

        /**
         * @brief Checks the value of the dimension parameter.
         * @param[in] aDimension A string represents the dimension.
         * @return The corresponding dimension value.
         *
         */
        static exageostat::common::Dimension CheckDimensionValue(const std::string &aDimension);

        /**
         * @brief Checks if the kernel value is valid.
         * @param[in] aKernel The kernel to check.
         * @return void
         *
         */
        void CheckKernelValue(const std::string &aKernel);

        /**
        * @brief Check input computation value.
        * @param[in] aValue The input from the user side.
        * @return Enum with the selected computation, Error if not exist.
        *
        */
        static common::Computation CheckComputationValue(const std::string &aValue);

        /**
         * @brief Check input precision value.
         * @param[in] aValue The input from the user side.
         * @return Enum with the selected Precision, Error if not exist.
         *
         */
        static common::Precision CheckPrecisionValue(const std::string &aValue);

        /**
         * @brief Checks the value of the unknown observations parameter.
         * @param[in] aValue A string represents the number of unknown observations.
         * @return The corresponding integer value.
         *
         */
        int CheckUnknownObservationsValue(const std::string &aValue);

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

        /**
         * @brief Parses a string of theta values and returns an array of doubles.
         * @param[in] aInputValues The input string of theta values.
         * @return A vector of parsed theta values.
         *
         */
        static std::vector<double> ParseTheta(const std::string &aInputValues);


        /**
         * @brief parse user's input to distance metric.
         * @param[in] aDistanceMetric string specifying the used distance metric.
         * @return void
         *
         */
        void ParseDistanceMetric(const std::string &aDistanceMetric);

    private:

        /**
         * @brief Checks the run mode and sets the verbosity level.
         * @param[in] aVerbosity A string represents the desired run mode ("verbose" or "standard").
         * @throws std::range_error if the input string is not "verbose" or "standard".
         * @return void
         *
         */
        static void ParseVerbose(const std::string &aVerbosity);

        /**
         * @brief Checks if a given string is in camel case format.
         * @param[in] aString The string to check.
         * @return true if the string is in camel case format, false otherwise.
         *
         */
        static bool IsCamelCase(const std::string &aString);

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
