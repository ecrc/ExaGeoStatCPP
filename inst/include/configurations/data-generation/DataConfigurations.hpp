
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

            protected:
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
            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat
#endif //EXAGEOSTAT_CPP_DATACONFIGURATIONS_HPP