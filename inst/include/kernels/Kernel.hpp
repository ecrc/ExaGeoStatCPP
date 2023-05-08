
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Kernels.hpp
 * @brief Header file for the Kernels class, which contains the main kernel functions.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-05
**/

#ifndef EXAGEOSTAT_CPP_KERNELS_HPP
#define EXAGEOSTAT_CPP_KERNELS_HPP

#include <data-units/Locations.hpp>
#include <starpu.h>
#include <gsl/gsl_sf_bessel.h>
#include <common/PluginRegistry.hpp>
#include<cmath>

/*
* @def EARTH_RADIUS
* @brief The radius of the Earth in kilometers.
*
* This macro defines the radius of the Earth in kilometers, which is used by the Great Circle Distance (GCD) function.
*/
 #define EARTH_RADIUS 6371.0


namespace exageostat {
    namespace kernels {

        /**
         * @class Kernels
         * @brief A base class for kernel functions.
         *
         * This class provides a base class for kernel functions and contains several utility functions for computing
         * distance metrics and Bessel functions.
         *
         */
        class Kernel {
        public:

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @param[in] apMatrixA The output covariance matrix.
             * @param[in] aRowsNumber The number of rows in the output matrix.
             * @param[in] aColumnsNumber The number of columns in the output matrix.
             * @param[in] aRowOffset The row offset for the input locations.
             * @param[in] aColumnOffset The column offset for the input locations.
             * @param[in] apLocation1 The set of input locations 1.
             * @param[in] apLocation2 The set of input locations 2.
             * @param[in] apLocation3 The set of input locations 3.
             * @param[in] apLocalTheta An array of kernel parameters.
             * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
             */
            virtual void
            GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber, int aRowOffset,
                                     int aColumnOffset, dataunits::Locations *apLocation1,
                                     dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                     double *apLocalTheta, int aDistanceMetric) = 0;

            /**
             * @brief Calculates the Euclidean distance between two points.
             * @param apLocations1 Pointer to the first set of locations.
             * @param apLocations2 Pointer to the second set of locations.
             * @param aIdxLocation1 Index of the first location in the first set.
             * @param aIdxLocation2 Index of the second location in the second set.
             * @param aDistanceMetric Flag indicating the distance metric to use (1 for Manhattan distance, 2 for Euclidean distance).
             * @param aFlagZ Flag indicating whether the points are in 2D or 3D space (0 for 2D, 1 for 3D).
             * @return The Euclidean distance between the two points.
             */
            static double
            CalculateDistance(dataunits::Locations *apLocations1, dataunits::Locations *apLocations2, int aIdxLocationX,
                              int aIdxLocationY,
                              int aDistanceMetric, int aFlagZ);

            /**
             * @brief Calculates the great-circle distance between two points on Earth using the Haversine formula.
             * @param aLatitude1 Latitude of the first point in degrees.
             * @param aLongitude1 Longitude of the first point in degrees.
             * @param aLatitude2 Latitude of the second point in degrees.
             * @param aLongitude2 Longitude of the second point in degrees.
             * @return The distance between the two points in kilometers.
             */
            static double DistanceEarth(double aLatitude1, double aLongitude1, double aLatitude2, double aLongitude2);

            /**
             * @brief Calculates the derivative of the modified Bessel function of the second kind (K_nu)
             * with respect to its input, evaluated at input_value and order aOrder.
             * @param aOrder The order of the Bessel function.
             * @param input_value The input value at which to evaluate the derivative.
             * @return The value of the derivative of K_nu with respect to its input, evaluated at input_value and order aOrder.
             */
            static double CalculateDerivativeBesselInputNu(const double &aOrder, const double &aInputValue);

            /**
             * @brief Calculates the derivative of the modified Bessel function of the second kind (K_nu) with respect to its order,
             * evaluated at input_value and order aOrder.
             * @param aOrder The order of the Bessel function.
             * @param input_value The input value at which to evaluate the derivative.
             * @return The value of the derivative of K_nu with respect to its order, evaluated at input_value and order aOrder.
             */
            static double CalculateDerivativeBesselNu(const double &aOrder, const double &aInputValue);

            /**
             * @brief Calculates the second derivative of the modified Bessel function of the second kind (K_nu) with respect to its input,
             * evaluated at input_value and order aOrder.
             * @param aOrder The order of the Bessel function.
             * @param input_value The input value at which to evaluate the second derivative.
             * @return The value of the second derivative of K_nu with respect to its input, evaluated at input_value and order aOrder.
             */
            static double CalculateSecondDerivativeBesselNu(const double &aOrder, const double &aInputValue);

            /**
             * Calculates the second derivative of the modified Bessel function of the second kind (K_nu) with respect to its input, evaluated at input_value and order aOrder.
             * @param aOrder The order of the Bessel function.
             * @param input_value The input value at which to evaluate the derivative.
             * @return The value of the derivative of K_nu with respect to its input, evaluated at input_value and order aOrder.
             */
            static double CalculateSecondDerivativeBesselNuInput(const double &aOrder, const double &aInputValue);

            /**
             * @brief Returns the value of the parameter P used by the kernel function.
             * @return The value of P.
             */
            int GetPValue() const;

        protected:
            //// Used P.
            int mP = 1;

        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_KERNELS_HPP