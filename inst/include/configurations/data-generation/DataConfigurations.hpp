
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataConfigurations.hpp
 *
 * @version 1.0.0
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

            class DataConfigurations : public Configurations {

            public:
                /**
                 * @brief Kernel setter.
                 * @param aKernel
                 */
                void
                SetKernel(std::string aKernel);

                /**
                 * @brief Kernel getter.
                 * @return mKernel
                 */
                std::string
                GetKernel();

                /**
                 * @brief Check Kernel value.
                 * @param aKernel
                 */
                void CheckKernelValue(std::string aKernel);

                /**
                 * @brief Kernel setter.
                 * @param aKernel
                 */
                void
                SetIsSynthetic(bool aIsSynthetic);

                /**
                 * @brief Kernel getter.
                 * @return mKernel
                 */
                bool
                GetIsSynthetic();

            protected:
                /// Used Kernel.
                std::string mKernel;
                /// Used type of data.
                bool mIsSynthetic = false;
            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATACONFIGURATIONS_HPP
