
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DistanceCalculationHelpers.hpp
 * @brief Contains the definition of the DistanceCalculationHelpers class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#ifndef EXAGEOSTATCPP_DistanceCalculationHelpers_HPP
#define EXAGEOSTATCPP_DistanceCalculationHelpers_HPP

#include <data-units/Locations.hpp>

namespace exageostat {
    namespace helpers {

        /**
         * @Class DistanceCalculationHelpers
         * @brief Class to calculate the distance between two points.
         * @tparam T Data Type: float or double.
         */

        template<typename T>
        class DistanceCalculationHelpers {
        public:
            /**
             * @brief Calculates the Euclidean distance between two points.
             * @param[in] aLocations1 Reference to the first set of locations.
             * @param[in] aLocations2 Reference to the second set of locations.
             * @param[in] aIdxLocation1 Index of the first location in the first set.
             * @param[in] aIdxLocation2 Index of the second location in the second set.
             * @param[in] aDistanceMetric Flag indicating the distance metric to use (1 for Manhattan distance, 2 for Euclidean distance).
             * @param[in] aFlagZ Flag indicating whether the points are in 2D or 3D space (0 for 2D, 1 for 3D).
             * @return The Euclidean distance between the two points.
             */
            static T CalculateDistance(exageostat::dataunits::Locations<T> &aLocations1,
                                       exageostat::dataunits::Locations<T> &aLocations2, const int &aIdxLocation1,
                                       const int &aIdxLocation2, const int &aDistanceMetric, const int &aFlagZ);

            /**
             * @brief Calculates the great-circle distance between two points on Earth using the Haversine formula.
             * @param[in] aLatitude1 Latitude of the first point in degrees.
             * @param[in] aLongitude1 Longitude of the first point in degrees.
             * @param[in] aLatitude2 Latitude of the second point in degrees.
             * @param[in] aLongitude2 Longitude of the second point in degrees.
             * @return The distance between the two points in kilometers.
             */
            static T DistanceEarth(T &aLatitude1, T &aLongitude1, T &aLatitude2, T &aLongitude2);

        };
        /**
          * @brief Instantiates the PredictionHelpers class for float and double types.
          * @tparam T Data Type: float or double
          *
          */
        EXAGEOSTAT_INSTANTIATE_CLASS(DistanceCalculationHelpers)
    }
}
#endif //EXAGEOSTATCPP_DistanceCalculationHelpers_HPP
