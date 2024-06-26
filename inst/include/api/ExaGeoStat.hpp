
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.hpp
 * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTAT_HPP
#define EXAGEOSTATCPP_EXAGEOSTAT_HPP

#include <nlopt.hpp>

#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::api {
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
         * @param[in] aConfigurations Reference to Configurations object containing user input data.
         * @param[out] aData Reference to an ExaGeoStatData<T> object where generated data will be stored.
         * @return void
         *
         */
        static void ExaGeoStatLoadData(configurations::Configurations &aConfigurations,
                                       std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Models Data whether it's synthetic data or real.
         * @param[in] aConfigurations Reference to Configurations object containing user input data.
         * @param[in] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
         * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
         * @return the last optimum value of MLE.
         *
         */
        static T ExaGeoStatDataModeling(configurations::Configurations &aConfigurations,
                                        std::unique_ptr<ExaGeoStatData<T>> &aData,
                                        T *apMeasurementsMatrix = nullptr);


        /**
         * @brief Objective function used in optimization, and following the NLOPT objective function format.
         * @param[in] aTheta An array of length n containing the current point in the parameter space.
         * @param[in] aGrad  An array of length n where you can optionally return the gradient of the objective function.
         * @param[in] apInfo pointer containing needed configurations and data.
         * @return double MLE results.
         *
         */
        static double ExaGeoStatMLETileAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo);

        /**
         * @brief Predict missing measurements values.
         * @param[in] aConfigurations Reference to Configurations object containing user input data.
         * @param[in, out] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
         * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
         * @param[in] apTrainLocations (Optional) Pointer to Locations represents training locations. these are used in training phase.
         * @param[in] apTestLocations (Optional) Pointer to Locations represents test locations. These are used in prediction phase.
         * @return void
         *
         */
        static void
        ExaGeoStatPrediction(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData,
                             T *apMeasurementsMatrix = nullptr, dataunits::Locations<T> *apTrainLocations = nullptr,
                             dataunits::Locations<T> *apTestLocations = nullptr);

    };

    /**
     * @brief Instantiates the ExaGeoStat class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStat)
}

#endif //EXAGEOSTATCPP_EXAGEOSTAT_HPP