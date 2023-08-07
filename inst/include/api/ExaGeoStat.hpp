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

#include <nlopt.h>

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

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
            static void ExaGeoStatInitializeHardware(common::Computation aComputation, int aCoreNumber, int aGpuNumber);

            /**
             * @brief Finalizes the ExaGeoStat hardware.
             * @param[in] apConfigurations Pointer to the configurations object.
             * @return void
             *
             */
            static void ExaGeoStatFinalizeHardware(common::Computation aComputation, dataunits::DescriptorData<T> *apDescriptorData);

            /**
             * @brief Generates Data whether it's synthetic data or real.
             * @param[in] apConfigurations Pointer to the configurations object.
             * @return void
             *
             */
            static exageostat::dataunits::ExaGeoStatData<T> *ExaGeoStatGenerateData(configurations::Configurations *apConfigurations);

            /**
             * @brief Models Data whether it's synthetic data or real.
             * @param[in] apDataLocations Pointer to Data Locations.
             * @return void
             *
             */
            static void ExaGeoStatDataModeling(configurations::Configurations *apConfigurations, exageostat::dataunits::ExaGeoStatData<T> *apData );

            /**
             * @brief Objective function used in optimization, and following the NLOPT objective function format.
             * @param aN The number of dimensions in the optimization problem.
             * @param apTheta An array of length n containing the current point in the parameter space.
             * @param apGrad  An array of length n where you can optionally return the gradient of the objective function.
             * @param apInfo pointers containing needed configurations and data.
             * @return double MLE results.
             */
            static double ExaGeoStatMleTileAPI(unsigned aN, const double *apTheta, double *apGrad, void *apInfo);

            /** @brief: Initialize the NLOPT optimizer
             * @param[out] apOpt: nlopt_opt object.
             * @param[in] apLowerBound: optimization lower bounds vector ( lb_1, lb_2, lb_3).
             * @param[in] apUpperBound: optimization upper bounds vector ( ub_1, ub_2, ub_3).
             * @param[in] aTolerance: a tolerance that is used for the purpose of stopping criteria only.
             * @return nlopt_opt object.
             */
            static void ExaGeoStatInitOptimizer(nlopt_opt *apOpt, double *apLowerBound, double *apUpperBound, double aTolerance);

        private:
            /**
             * @brief
             * Prevent Class Instantiation for API Wrapper Class.
             *
             */
            ExaGeoStat() = default;

            struct mModelingData{
                dataunits::ExaGeoStatData<T> *Data;
                configurations::Configurations *Configuration;
            };
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