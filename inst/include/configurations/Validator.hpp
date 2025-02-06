
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
* @file Validator.hpp
* @version 1.1.0
* @brief Contains the declaration of the Validator class and its member functions for configuration validation.
* @details Provides a set of static methods to validate the input configuration parsed from CLI or JSON for the ExaGeoStat software package.
* @author Mahmoud ElKarargy
* @date 2024-12-11
**/

#ifndef EXAGEOSTAT_CPP_VALIDATOR_HPP
#define EXAGEOSTAT_CPP_VALIDATOR_HPP

#include <functional>

namespace exageostat::configurations::validator {

    /**
     * @class Validator
     * @brief A class containing static methods for validating configuration parameters.
     */
    class Validator {
    public:

        /**
         * @brief Validates the configuration parameters parsed from the CLI or JSON input.
         * @param[in,out] apArgsMap A map containing configuration parameters and their values.
         * @return void
         *
         */
        static void Validate(std::unordered_map<std::string, std::any> &apArgsMap);

        /**
         * @brief Parses a string of theta values and converts them to a vector of doubles.
         * @param[in] aInputValues The input string containing theta values.
         * @return A vector of parsed theta values as doubles.
         *
         */
        static std::vector<double> CheckThetaValue(const std::string &aInputValues);

        /**
         * @brief Validates and converts a string representing a tolerance value into a double.
         * @param[in] aTolerance The input string for the tolerance value.
         * @return The validated and converted tolerance value as a double.
         *
         */
        static double CheckToleranceValue(const std::string& aTolerance);

        /**
         * @brief Validates a string representing a file name and returns file handler.
         * @param[in] aFileLogPath The input string for the file name.
         * @return The file handler of the validated file name.
         *
         */
        static FILE *CheckLogFileValue(const std::string &aFileLogPath);

        /**
         * @brief Validates a string representing a kernel name.
         * @param[in] aKernel The input string for the kernel name.
         * @return The validated kernel name as a string.
         *
         */
        static std::string CheckKernelValue(const std::string &aKernel);

        /**
         * @brief Checks if the input string is in camel case format.
         * @param[in] aString The input string to validate.
         * @return True if the string is in camel case; otherwise, false.
         *
         */
        bool static IsCamelCase(const std::string &aString);

        /**
         * @brief Validates and converts a verbosity level string into a Verbose enum.
         * @param[in] aVerbosity The input string for the verbosity level.
         * @return The validated Verbose enum value.
         *
         */
        static common::Verbose CheckVerboseValue(const std::string &aVerbosity);

        /**
         * @brief Validates and converts a string into a Precision enum.
         * @param[in] aValue The input string representing precision.
         * @return The validated Precision enum value.
         *
         */
        static common::Precision CheckPrecisionValue(const std::string &aValue);

        /**
         * @brief Validates and converts a string into a Computation enum.
         * @param[in] aValue The input string representing computation type.
         * @return The validated Computation enum value.
         *
         */
        static common::Computation CheckComputationValue(const std::string &aValue);

        /**
         * @brief Validates and converts a string into a boolean value.
         * @param[in] aBooleanValue The input string representing a boolean value.
         * @return The validated boolean value.
         *
         */
        static bool CheckBoolValue(const std::string& aBooleanValue);

        /**
         * @brief Validates and converts a string into a DistanceMetric enum.
         * @param[in] aDistanceMetric The input string representing the distance metric.
         * @return The validated DistanceMetric enum value.
         *
         */
        static common::DistanceMetric CheckDistanceMetricValue(const std::string &aDistanceMetric);

        /**
         * @brief Validates and converts a string into an integer.
         * @param[in] aValue The input string representing a numerical value.
         * @return The validated integer value.
         *
         */
        static int CheckNumericalValue(const std::string &aValue);

        /**
         * @brief Validates and converts a string into a Dimension enum.
         * @param[in] aDimension The input string representing the dimension.
         * @return The validated Dimension enum value.
         *
         */
        static common::Dimension CheckDimensionValue(const std::string &aDimension);

    private:
        /// A map of validation functions for different parameter types.
        static const std::unordered_map<std::string, std::function<std::any(const std::string&)>> mCheckersMap;
        /// A map linking arguments to their categories.
        static const std::unordered_map<std::string, std::string> mArgumentToCategoryMap;
    };
}// namespace exageostat::configurations::validator

#endif // EXAGEOSTAT_CPP_VALIDATOR_HPP
