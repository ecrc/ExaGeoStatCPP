// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.hpp
 * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-30
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTAT_HPP
#define EXAGEOSTATCPP_EXAGEOSTAT_HPP

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <data-units/Locations.hpp>

namespace exageostat {
    namespace api {
        /**
         * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
         * @tparam T Data Type: float or double
         *
         */
        template<typename T>
        class ExaGeoStat {
        public:
            /**
             * @brief Initializes the ExaGeoStat hardware.
             * @param[in] apConfigurations Pointer to the configurations object.
             * @return void
             *
             */
            static void ExaGeoStatInitializeHardware();

            /**
             * @brief Finalizes the ExaGeoStat hardware.
             * @param[in] apConfigurations Pointer to the configurations object.
             * @return void
             *
             */
            static void ExaGeoStatFinalizeHardware();

            /**
             * @brief Generates Data whether it's synthetic data or real.
             * @param[in] apConfigurations Pointer to the configurations object.
             * @return void
             *
             */
            static dataunits::Locations * ExaGeoStatGenerateData();

            /**
             * @brief Models Data whether it's synthetic data or real.
             * @param[in] apDataLocations Pointer to Data Locations.
             * @return void
             *
             */
            static void ExaGeoStatDataModeling(dataunits::Locations *apDataLocations);

        private:
            /**
             * @brief
             * Prevent Class Instantiation for API Wrapper Class.
             *
             */
            ExaGeoStat() = default;
        };

        /**
         * @brief Instantiates the ExaGeoStat class for float and double types.
         * @tparam T Data Type: float or double
         *
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStat)
    }
}

#endif //EXAGEOSTATCPP_EXAGEOSTAT_HPP