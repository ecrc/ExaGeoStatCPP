
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
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
         *
         * @brief A class for configuring data settings in ExaGeoStat.
         */
            class DataConfigurations : public Configurations {

            public:
                /**
                 * @brief Setter for the kernel.
                 *
                 * @param aKernel The kernel to set.
                 */
                void SetKernel(const std::string &aKernel);

                /**
                 * @brief Getter for the kernel.
                 *
                 * @return The kernel.
                 */
                std::string GetKernel();

                /**
                 * @brief Checks if the kernel value is valid.
                 *
                 * @param aKernel The kernel to check.
                 */
                void CheckKernelValue(const std::string& aKernel);
                static bool IsCamelCase(std::string aString);

                /**
                 * @brief Setter for the data type.
                 *
                 * @param aIsSynthetic The type of data to set.
                 */
                void SetIsSynthetic(bool aIsSynthetic);

                /**
                 * @brief Getter for the data type.
                 *
                 * @return The data type.
                 */
                bool GetIsSynthetic() const;

            protected:
                /// The kernel to use.
                std::string mKernel;
                /// The type of data to use.
                bool mIsSynthetic = false;
            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATACONFIGURATIONS_HPP
