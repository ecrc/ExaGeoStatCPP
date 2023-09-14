
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PredictionHelpers.hpp
 * @brief Contains the definition of the PredictionHelpers.hpp class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <data-units/Locations.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <configurations/Configurations.hpp>

#ifndef EXAGEOSTATCPP_PREDICTION_HELPERS_HPP
#define EXAGEOSTATCPP_PREDICTION_HELPERS_HPP

namespace exageostat {
    namespace prediction {
        /**
         * @Class PredictionHelpers
         * @brief Class to define and implement different Prediction Module helpers functions.
         * @tparam T Data Type: float or double.
         */

        template<typename T>
        class PredictionHelpers {
        public:
            /**
            * @brief Pick random Z points for prediction depending on p.
            * @param[in] aConfigurations Configurations object containing relevant settings.
            * @param[in,out] aData Reference ExaGeoStatData object populated with locations and descriptor data.
            * @param[out] apZObs Pointer to be filled with observed measurements.
            * @param[out] apZActual Pointer to be filled with actual measurements.
            * @param[in] apZ Pointer to a copy of the measurements matrix.
            * @param[out] aMissLocation Location object to be filled with missed locations.
            * @param[out] aObsLocation Location object to be filled with missed locations.
            * @return void
            */
            static void
            PickRandomPoints(exageostat::configurations::Configurations &aConfigurations,
                             exageostat::dataunits::ExaGeoStatData<T> &aData, T *apZObs, T *apZActual, T *apZ,
                             exageostat::dataunits::Locations<T> &aMissLocation,
                             exageostat::dataunits::Locations<T> &aObsLocation);

            /**
             * @brief Shuffle array.
             * @param[in, out] apArray Array to be shuffled.
             * @param[in, out] aLocations Locations to be shuffled.
             * @param[out] aSize Size of data.
             * @return void
             */
            static void
            Shuffle(T *apArray, exageostat::dataunits::Locations<T> &aLocations, int aSize);

            /**
             * @brief Shuffle array.
             * @param[in, out] apArray1 First Array to be shuffled.
             * @param[in, out] apArray2 Second Array to be shuffled.
             * @param[in, out] aLocations Locations to be shuffled.
             * @param[out] aSize Size of data.
             * @return void
             */
            static void
            Shuffle(T *apArray1, T *apArray2, exageostat::dataunits::Locations<T> &aLocations, int aSize);

            /**
             * @brief Shuffle array.
             * @param[in, out] apArray1 First Array to be shuffled.
             * @param[in, out] apArray2 Second Array to be shuffled.
             * @param[in, out] apArray3 Third Array to be shuffled.
             * @param[in, out] aLocations Locations to be shuffled.
             * @param[out] aSize Size of data.
             * @return void
             */
            static void
            Shuffle(T *apArray1, T *apArray2, T *apArray3, exageostat::dataunits::Locations<T> &aLocations, int aSize);

            /**
             * @brief Sorts the input data using the C++, the pre-defined standard library function sort().
             * @param aData[in] Pointer to the array of data to be sorted.
             * @param aCount[in] Number of elements in the input array.
             * @param aDimension[in] Dimension of the input data.
             * @return void
             */
            static void
            SortArray(uint32_t *aData, int aCount);

            /**
             *@brief Sorts location data and corresponding observation values in-place based on Locations coordinates.
             *
             * @param[in] aN Number of data points (input).
             * @param[in,out] aLocations Reference to the Locations object containing X and Y coordinates (input/output).
             * @param[in,out] apZ Pointer to the array containing observation values (input/output).
             * @return 0 if the sorting is successful.
             */
            static int
            SortInplace(int aN, exageostat::dataunits::Locations<T> &aLocations, T *apZ);

        };

        /**
         * @brief Instantiates the PredictionHelpers class for float and double types.
         * @tparam T Data Type: float or double
         *
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(PredictionHelpers)
    }
}

#endif //EXAGEOSTATCPP_PREDICTION_HELPERS_HPP