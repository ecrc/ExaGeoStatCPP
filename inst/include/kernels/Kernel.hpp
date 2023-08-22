
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).


/**
 * @file Kernels.hpp
 * @brief Header file for the Kernels class, which contains the main kernel functions.
 * @version 1.0.0
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-03-05
**/

#ifndef EXAGEOSTAT_CPP_KERNELS_HPP
#define EXAGEOSTAT_CPP_KERNELS_HPP

#include<cmath>

extern "C" {
#include <starpu.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_psi.h>
}

#include <common/PluginRegistry.hpp>
#include <data-units/Locations.hpp>

/**
 * @def EARTH_RADIUS
 * @brief The radius of the Earth in kilometers.
 * @details This macro defines the radius of the Earth in kilometers, which is used by the Great Circle Distance (GCD) function.
 *
 */
#define EARTH_RADIUS 6371.0

/**
 * @def POW_e
 * @brief The value of e to the power of 1.
 * @details This macro defines the value of e to the power of 1, which is used in some kernel functions.
 *
 */
#define POW_e (2.71828182845904)


namespace exageostat {
    namespace kernels {

        /**
         * @class Kernels
         * @brief A base class for kernel functions.
         * @details This class provides a base class for kernel functions and contains several utility functions for computing distance metrics and Bessel functions.
         *
         */
        template<typename T>
        class Kernel {
        public:

            /**
             * Default virtual destructor to be overidden by the the suitable concrete kernel destructor.
             */
            virtual ~Kernel() = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @param[out] apMatrixA The output covariance matrix.
             * @param[in] aRowsNumber The number of rows in the output matrix.
             * @param[in] aColumnsNumber The number of columns in the output matrix.
             * @param[in] aRowOffset The row offset for the input locations.
             * @param[in] aColumnOffset The column offset for the input locations.
             * @param[in] apLocation1 The set of input locations 1.
             * @param[in] apLocation2 The set of input locations 2.
             * @param[in] apLocation3 The set of input locations 3.
             * @param[in] aLocalTheta An array of kernel parameters.
             * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
             * @return void
             */
            virtual void
            GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                     int &aColumnOffset, dataunits::Locations<T> *apLocation1,
                                     dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                     T *aLocalTheta, int &aDistanceMetric) = 0;

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
            T CalculateDistance(dataunits::Locations<T> &aLocations1, dataunits::Locations<T> &aLocations2,
                                int &aIdxLocation1, int &aIdxLocation2, int &aDistanceMetric, int &aFlagZ);

            /**
             * @brief Calculates the great-circle distance between two points on Earth using the Haversine formula.
             * @param[in] aLatitude1 Latitude of the first point in degrees.
             * @param[in] aLongitude1 Longitude of the first point in degrees.
             * @param[in] aLatitude2 Latitude of the second point in degrees.
             * @param[in] aLongitude2 Longitude of the second point in degrees.
             * @return The distance between the two points in kilometers.
             */
            static T DistanceEarth(T &aLatitude1, T &aLongitude1, T &aLatitude2, T &aLongitude2);

            /**
             * @brief Calculates the derivative of the modified Bessel function of the second kind (K_nu) with respect to its input, evaluated at input_value and order aOrder.
             * @param[in] aOrder The order of the Bessel function.
             * @param[in] aInputValue The input value at which to evaluate the derivative.
             * @return The value of the derivative of K_nu with respect to its input, evaluated at input_value and order aOrder.
             */
            static T CalculateDerivativeBesselInputNu(const T &aOrder, const T &aInputValue);

            /**
             * @brief Calculates the derivative of the modified Bessel function of the second kind (K_nu) with respect to its order, evaluated at input_value and order aOrder.
             * @param[in] aOrder The order of the Bessel function.
             * @param[in] aInputValue The input value at which to evaluate the derivative.
             * @return The value of the derivative of K_nu with respect to its order, evaluated at input_value and order aOrder.
             *
             */
            static T CalculateDerivativeBesselNu(const T &aOrder, const T &aInputValue);

            /**
             * @brief Calculates the second derivative of the modified Bessel function of the second kind (K_nu) with respect to its input, evaluated at input_value and order aOrder.
             * @param[in] aOrder The order of the Bessel function.
             * @param[in] aInputValue The input value at which to evaluate the second derivative.
             * @return The value of the second derivative of K_nu with respect to its input, evaluated at input_value and order aOrder.
             *
             */
            static T CalculateSecondDerivativeBesselNu(const T &aOrder, const T &aInputValue);

            /**
             * @brief Calculates the second derivative of the modified Bessel function of the second kind (K_nu) with respect to its input, evaluated at input_value and order aOrder.
             * @param[in] aOrder The order of the Bessel function.
             * @param[in] aInputValue The input value at which to evaluate the derivative.
             * @return The value of the derivative of K_nu with respect to its input, evaluated at input_value and order aOrder.
             *
             */
            static T CalculateSecondDerivativeBesselNuInput(const T &aOrder, const T &aInputValue);

            /**
             * @brief Returns the value of the parameter P used by the kernel function.
             * @return The value of P.
             *
             */
            int GetPValue() const;

            /**
             * @brief Sets the value of the parameter P used by the kernel function.
             * @param[in] aP Value to set `mP` with.
             * @return void
             *
             */
            void SetPValue(int aP);

            /**
             * @brief Returns the number of the parameters used by the kernel function.
             * @return The value of ParametersNumber.
             *
             */
            int GetParametersNumbers() const;

        protected:
            //// Used P.
            int mP = 1;
            //// Used number of paramters.
            int mParametersNumber = 3;
        };

        /**
         * @brief Instantiates the Data Generator class for float and double types.
         * @tparam T Data Type: float or double
         *
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(Kernel)

    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_KERNELS_HPP