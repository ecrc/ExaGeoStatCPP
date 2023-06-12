// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-30
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTAT_HPP
#define EXAGEOSTATCPP_EXAGEOSTAT_HPP

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>

namespace exageostat {
    namespace api {
        /**
         * @brief
         * High-Level Wrapper class containing the static API for ExaGeoStat operations.
         *
         * @tparam T
         * Data type: float or double
         */
        template<typename T>
        class ExaGeoStat {
        public:
            /**
            * @brief Initializes the ExaGeoStat hardware.
            *
            * @param apConfigurations Pointer to the configurations object.
            */
            static void ExaGeoStatInitializeHardware(configurations::Configurations *apConfigurations);

            /**
             * @brief Finalizes the ExaGeoStat hardware.
             *
             * @param apConfigurations Pointer to the configurations object.
             */
            static void ExaGeoStatFinalizeHardware(configurations::Configurations *apConfigurations);

            /**
             * @brief Generates Data whether it's synthetic data or real.
             *
             * @param apConfigurations Pointer to the configurations object.
             */
            static void ExaGeoStatGenerateData(configurations::Configurations *apConfigurations);

        private:
            /**
            * @brief
            * Prevent Class Instantiation for API Wrapper Class.
            */
            ExaGeoStat() = default;
        };

        /**
         * @brief Instantiates the ExaGeoStat class for float and double types.
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStat)
    }
}

#endif //EXAGEOSTATCPP_EXAGEOSTAT_HPP
