
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PredictionAuxiliaryFunctions.hpp
 * @brief Contains the definition of the PredictionAuxiliaryFunctions.hpp class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <data-units/Locations.hpp>

#ifndef EXAGEOSTATCPP_PREDICTION_AUXILIARY_FUNCTIONS_HPP
#define EXAGEOSTATCPP_PREDICTION_AUXILIARY_FUNCTIONS_HPP

namespace exageostat::prediction {

    /**
     * @Class PredictionAuxiliaryFunctions
     * @brief Class to define and implement different Prediction Module Auxiliary Functions.
     * @tparam T Data Type: float or double.
     */
    template<typename T>
    class PredictionAuxiliaryFunctions {
    public:

        /**
         * @brief Default constructor for PredictionAuxiliaryFunctions
         */
        PredictionAuxiliaryFunctions() = default;

        /**
         * @brief Default destructor for PredictionAuxiliaryFunctions
         */
        ~PredictionAuxiliaryFunctions() = default;

        /**
         * @brief  implements the Inverse Distance Weighting (IDW) interpolation method
         * for predicting missing values based on available observed values.
         * @param[in] apZMiss Pointer to the missed measurements.
         * @param[in] apZActual Pointer to the actual measurements.
         * @param[in] apZObs Pointer to the observed measurements.
         * @param[in] aZMissNumber Number of missed measurements.
         * @param[in] aZObsNumber Number of observed measurements.
         * @param[in] aMissLocation Reference to the missed locations.
         * @param[in] aObsLocation Reference to the observed locations.
         * @param[out] apMSPE Pointer to be filled with MSPE value.
         * @return T Array provides insight into the accuracy of the IDW-interpolated predictions for missing values
         */
        static void PredictIDW(T *apZMiss, T *apZActual, T *apZObs, const int &aZMissNumber, const int &aZObsNumber,
                               exageostat::dataunits::Locations<T> &aMissLocation,
                               exageostat::dataunits::Locations<T> &aObsLocation, T *apMSPE);

    };

    /**
      * @brief Instantiates the PredictionAuxiliaryFunctions class for float and double types.
      * @tparam T Data Type: float or double
      *
      */
    EXAGEOSTAT_INSTANTIATE_CLASS(PredictionAuxiliaryFunctions)
}

#endif //EXAGEOSTATCPP_PREDICTION_AUXILIARY_FUNCTIONS_HPP