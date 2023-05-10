
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonStat.hpp
 * @brief Defines the UnivariateMaternNonStat class, a Univariate Matern Non Stat kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNNONSTAT_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNNONSTAT_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternNonStat
         * @brief A class representing a Univariate Matern Non Stat kernel.
         *
         * This class represents a Univariate Matern Non Stat kernel, which is a subclass of the Kernel class. It provides
         * a method for generating a covariance matrix using a set of input locations and kernel parameters.
         */
        class UnivariateMaternNonStat : public Kernel {

        public:

            /**
             * @brief Constructs a new UnivariateMaternNonStat object.
             *
             * Initializes a new UnivariateMaternNonStat object with default values.
             */
            UnivariateMaternNonStat();

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
            void GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber, int aRowOffset,
                                          int aColumnOffset, dataunits::Locations *apLocation1,
                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                          double *apLocalTheta, int aDistanceMetric) override ;
            /**
             * @brief Creates a new UnivariateMaternNonStat object.
             * @return A pointer to the new UnivariateMaternNonStat object.
             *
             * This method creates a new UnivariateMaternNonStat object and returns a pointer to it.
             */
            static Kernel *Create();

            /**
             * Function for smoothness parameter
             * @param x x co-ordinate
             * @param y y co-ordinate
             * @param g parameter for function
             * @param h parameter for function
             * @param ti parameter for function
             * @return The function ge^(h(x+y)) + ti
             */
            static double Neu(double x, double y, double g, double h, double ti);

            /**
             * Function for partial sill
             * @param x x co-ordinate
             * @param y y co-ordinate
             * @param d parameter for function
             * @param e parameter for function
             * @param f parameter for function
             * @return The function de^(e(x+y)) + f
             */
            static double Sigma(double x, double y, double d, double e, double f);

            /**
             * Function for spatial range
             * @param x x co-ordinate
             * @param y y co-ordinate
             * @param a parameter for function
             * @param b parameter for function
             * @return The function ae^(sin bx + sin by)
             */
            static double Lambda(double x, double y, double a, double b);

            /**
             * Returns the Mahalanobis distance between two points to account for anisotropy
             * @param x1 x co-ordinate of first point
             * @param y1 y co-ordinate of first point
             * @param x2 x co-ordinate of second point
             * @param y2 y co-ordinate of second point
             * @param a11 First element of the positive definite matrix that defines the Mahalanobis Distance
             * @param a12 Second element of the positive definite matrix that defines the Mahalanobis Distance
             * @param a21 Third element of the positive definite matrix that defines the Mahalanobis Distance
             * @param a22 Fourth element of the positive definite matrix that defines the Mahalanobis Distance
             * @return The Mahalanobis Distance
             */
            static double CalculateMahalanobisDistanceSquared(double x1, double y1, double x2,
                                                       double y2, double a11, double a12,
                                                       double a21, double a22);

            /**
             * Utility function that evaluates the matern. Similiar to
             * (https://www.rdocumentation.org/packages/PrevMap/versions/1.5.3/topics/matern.kernel) in R
             * @param range Spatial Range parameter (Also known as rho)
             * @param smoothness Smoothness parameter (Also known as neu)
             * @param distance Distance between the two locations
             * @return Matern function evaluation
             */
            static double MaternUtil(double range, double smoothness, double distance);

        private:
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNNONSTAT_HPP
