
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataConfigurations.hpp
 * @version 1.0.0
 * @brief Contains the definition of the DataConfigurations class for configuring data settings in ExaGeoStat.
 * @author Sameh Abdulah
 * @date 2023-02-03
**/

#ifndef EXAGEOSTAT_CPP_DATACONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_DATACONFIGURATIONS_HPP

#include <iostream>

#include <configurations/Configurations.hpp>

namespace exageostat {
    namespace configurations {
        namespace data_configurations {

            /**
             * @class DataConfigurations
             * @brief A class for configuring data settings in ExaGeoStat.
             *
             */
            class DataConfigurations : public Configurations {

            public:
                /**
                 * @brief Setter for the kernel.
                 * @param[in] aKernel The kernel to set.
                 * @return void
                 *
                 */
                void SetKernel(const std::string &aKernel);

                /**
                 * @brief Getter for the kernel.
                 * @return The kernel.
                 *
                 */
                std::string GetKernel() const;

                /**
                 * @brief Setter for the data type.
                 * @param[in] aIsSynthetic The type of data to set.
                 * @return void
                 *
                 */
                void SetIsSynthetic(bool aIsSynthetic);

                /**
                 * @brief Getter for the data type.
                 * @return The data type.
                 *
                 */
                bool GetIsSynthetic() const;

                /**
                 * @brief Setter for the number of parameters.
                 * @param[in] aParameterNumbers The number of parameters to set.
                 * @return void
                 *
                 */
                void SetParametersNumber(int aParameterNumbers);

                /**
                 * @brief Getter for the number of parameters.
                 * @return The number of parameters.
                 *
                 */
                int GetParametersNumber() const;

                /**
                 * @brief Setter for the lower bounds.
                 * @param[in,out] apTheta A pointer to an array of lower bounds to set.
                 * @return void
                 *
                 */
                void SetLowerBounds(std::vector<double> &apTheta);

                /**
                 * @brief Getter for the lower bounds.
                 * @return A pointer to the array of lower bounds.
                 *
                 */
                std::vector<double> &GetLowerBounds();

                /**
                 * @brief Setter for the upper bounds.
                 * @param[in,out] apTheta A pointer to an array of upper bounds to set.
                 * @return void
                 *
                 */
                void SetUpperBounds(std::vector<double> &apTheta);

                /**
                 * @brief Getter for the upper bounds.
                 * @return A pointer to the array of upper bounds.
                 *
                 */
                std::vector<double> &GetUpperBounds();

                /**
                 * @brief Setter for the starting theta.
                 * @param[in,out] apTheta A pointer to an array of starting theta values to set.
                 * @return void
                 *
                 */
                void SetStartingTheta(std::vector<double> &apTheta);

                /**
                 * @brief Getter for the starting theta.
                 * @return A pointer to the array of starting theta values.
                 *
                 */
                std::vector<double> &GetStartingTheta();

                /**
                 * @brief Setter for the initial theta.
                 * @param[in,out] apTheta A pointer to an array of initial theta values to set.
                 * @return void
                 *
                 */
                void SetInitialTheta(std::vector<double> &apTheta);

                /**
                 * @brief Getter for the initial theta.
                 * @return A pointer to the array of initial theta values.
                 *
                 */
                std::vector<double> &GetInitialTheta();

                /**
                 * @brief Setter for the target theta.
                 * @param[in,out] apTheta A pointer to an array of target theta values to set.
                 * @return void
                 *
                 */
                void SetTargetTheta(std::vector<double> &apTheta);

                /**
                 * @brief Getter for the target theta.
                 * @return A pointer to the array of target theta values.
                 *
                 */
                std::vector<double> &GetTargetTheta();

                /**
                 * @brief Checks if the kernel value is valid.
                 * @param[in] aKernel The kernel to check.
                 * @return void
                 *
                 */
                void CheckKernelValue(const std::string &aKernel);

                /**
                 * @brief Checks if a given string is in camel case format.
                 * @param[in] aString The string to check.
                 * @return true if the string is in camel case format, false otherwise.
                 *
                 */
                static bool IsCamelCase(const std::string &aString);

                /**
                 * @brief Checks the run mode and sets the verbosity level.
                 * @param[in] aRunMode A string representing the desired run mode ("verbose" or "standard").
                 * @throws std::range_error if the input string isnot "verbose" or "standard".
                 * @return void
                 *
                 */
                static void ParseRunMode(const std::string &aRunMode);

                /**
                 * @brief Parses a string of theta values and returns an array of doubles.
                 * @param[in] aInputValues The input string of theta values.
                 * @return A vector of parsed theta values.
                 *
                 */
                static std::vector<double> ParseTheta(const std::string &aInputValues);

            protected:
                /// The kernel to use.
                std::string mKernel;
                /// The type of data to use.
                bool mIsSynthetic = true;
                /// The lower bounds to use.
                std::vector<double> mLowerBounds;
                /// The upper bounds to use.
                std::vector<double> mUpperBounds;
                /// The starting theta values to use.
                std::vector<double> mStartingTheta;
                /// The target theta values to use.
                std::vector<double> mTargetTheta;
                /// The initial theta values to use.
                std::vector<double> mInitialTheta;
                //// The number of parameters to use.
                int mParametersNumber = 0;
            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat
#endif //EXAGEOSTAT_CPP_DATACONFIGURATIONS_HPP