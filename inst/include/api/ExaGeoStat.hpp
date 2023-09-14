
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.hpp
 * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-05-30
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTAT_HPP
#define EXAGEOSTATCPP_EXAGEOSTAT_HPP

#include <nlopt.hpp>

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

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
             * @brief Generates Data whether it's synthetic data or real.
             * @param[in] aHardware Reference to Hardware configuration for the ExaGeoStat solver.
             * @param[in] aConfigurations Reference to Configurations object containing user input data.
             * @param[out] aData Reference to an ExaGeoStatData<T> object where generated data will be stored
             * @return void
             *
             */
            static void ExaGeoStatGenerateData(const hardware::ExaGeoStatHardware &aHardware,
                                               configurations::Configurations &aConfigurations,
                                               dataunits::ExaGeoStatData<T> &aData);

            /**
             * @brief Models Data whether it's synthetic data or real.
             * @param[in] aHardware Reference to Hardware configuration for the ExaGeoStat solver.
             * @param[in] aConfigurations Reference to Configurations object containing user input data.
             * @param[in] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
             * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
             * @return void
             *
             */
            static T ExaGeoStatDataModeling(const hardware::ExaGeoStatHardware &aHardware,
                                            configurations::Configurations &aConfigurations,
                                            exageostat::dataunits::ExaGeoStatData<T> &aData,
                                            T *apMeasurementsMatrix = nullptr);


            /**
             * @brief Objective function used in optimization, and following the NLOPT objective function format.
             * @param[in] aTheta An array of length n containing the current point in the parameter space.
             * @param[in] aGrad  An array of length n where you can optionally return the gradient of the objective function.
             * @param[in] apInfo pointer containing needed configurations and data.
             * @return double MLE results.
             */
            static double
            ExaGeoStatMLETileAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo);

            /**
             * @brief Predict missing measurements values.
             * @param[in] aHardware Reference to Hardware configuration for the ExaGeoStat solver.
             * @param[in] aConfigurations Reference to Configurations object containing user input data.
             * @param[in, out] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
             * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
             * @return void
             */
            static void
            ExaGeoStatPrediction(const hardware::ExaGeoStatHardware &aHardware,
                                 configurations::Configurations &aConfigurations,
                                 exageostat::dataunits::ExaGeoStatData<T> &aData, T *apMeasurementsMatrix = nullptr);

            /**
             * @brief Predict missing measurements values.
             * @param[in] aHardware Reference to Hardware configuration for the ExaGeoStat solver.
             * @param[in] aConfigurations Reference to Configurations object containing user input data.
             * @param[in, out] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
             * @return void
             */
            static void
            ExaGeoStatPrediction(const hardware::ExaGeoStatHardware &aHardware,
                                 configurations::Configurations &aConfigurations,
                                 exageostat::dataunits::ExaGeoStatData<T> &aData);

        private:
            /**
             * @brief
             * Prevent Class Instantiation for API Wrapper Class.
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