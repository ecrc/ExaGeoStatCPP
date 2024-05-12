
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Prediction.hpp
 * @brief Contains the definition of the Prediction class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

#ifndef EXAGEOSTATCPP_PREDICTION_HPP
#define EXAGEOSTATCPP_PREDICTION_HPP

namespace exageostat::prediction {

    /**
     * @Class Prediction
     * @brief Class to handle different Prediction Module calls.
     * @tparam T Data Type: float or double.
     *
     */
    template<typename T>
    class Prediction {
    public:

        /**
        * @brief Takes care of calling the MSPE function, and the appropriate auxiliary function.
        * @param[in, out] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
        * @param[in] aConfigurations Reference to Configurations object containing user input data.
        * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
        * @param[in] aKernel Reference to the kernel object to use.
        * @param[in] apTrainLocations (Optional) Pointer to Locations represents training locations. these are used in training phase.
        * @param[in] apTestLocations (Optional) Pointer to Locations represents test locations. These are used in prediction phase.
        * @return void
        *
        */
        static void
        PredictMissingData(std::unique_ptr<ExaGeoStatData<T>> &aData, configurations::Configurations &aConfigurations,
                           T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel,
                           dataunits::Locations<T> *apTrainLocations = nullptr,
                           dataunits::Locations<T> *apTestLocations = nullptr);

        /**
         * @brief Initializes needed pointers for prediction.
         * @param[in] aConfigurations Reference to Configurations object containing user input data.
         * @param[in, out] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
         * @param[in] aLinearAlgebraSolver linear algebra solver depending on implementation.
         * @param[out] apZObs Pointer to be filled with observation measurements
         * @param[out] apZActual Pointer to be filled with actual measurements
         * @param[out] aMissLocation Location object to be filled with missed locations.
         * @param[out] aObsLocation Location object to be filled with missed locations.
         * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
         * @param[in] aP the P value of the kernel multiplied by time slot.
         * @param[in] apTrainLocations (Optional) Pointer to Locations represents training locations. these are used in training phase.
         * @param[in] apTestLocations (Optional) Pointer to Locations represents test locations. These are used in prediction phase.
         * @return void
         *
         */
        static void
        InitializePredictionArguments(configurations::Configurations &aConfigurations,
                                      std::unique_ptr<ExaGeoStatData<T>> &aData,
                                      std::unique_ptr<exageostat::linearAlgebra::LinearAlgebraMethods<T>> &aLinearAlgebraSolver,
                                      T *apZObs, T *apZActual, exageostat::dataunits::Locations<T> &aMissLocation,
                                      exageostat::dataunits::Locations<T> &aObsLocation, T *apMeasurementsMatrix,
                                      const int &aP, dataunits::Locations<T> *apTrainLocations,
                                      dataunits::Locations<T> *apTestLocations);

    };

    /**
     * @brief Instantiates the Prediction class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(Prediction)
}
#endif //EXAGEOSTATCPP_PREDICTION_HPP