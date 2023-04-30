
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternStationary.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-12
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNSTATIONARY_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNSTATIONARY_HPP

#include <kernels/Kernel.hpp>
#include <iostream>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>


namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternStationary
         * @brief A class representing a univariate stationary Matern kernel.
         */
        class UnivariateMaternStationary : public Kernel {

        public:

            UnivariateMaternStationary(){
                this->mP = 1;
            }
            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @param[in] apMatrixA The output covariance matrix.
             * @param[in] aRowsNumber The number of rows in the output matrix.
             * @param[in] aColumnsNumber The number of columns in the output matrix.
             * @param[in] aRowOffset The row offset for the input locations.
             * @param[in] aColumnOffset The column offset for the input locations.
             * @param[in] apLocation1 The set of input locations 1.
             * @param[in] apLocation2 The set of input locations 2.
             * @param[in] apLocalTheta An array of kernel parameters.
             * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
             */
            void GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber, int aRowOffset,
                                          int aColumnOffset, dataunits::Locations *apLocation1,
                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                          double *apLocalTheta, int aDistanceMetric) override {
                int i = 0, j = 0;
                int i0 = aRowOffset;
                int j0 = aColumnOffset;
                double x0, y0, z0;
                double expr = 0.0;
                double con = 0.0;
                double sigma_square = apLocalTheta[0];

                con = pow(2, (apLocalTheta[2] - 1)) * tgamma(apLocalTheta[2]);
                con = 1.0 / con;
                con = sigma_square * con;

                for (i = 0; i < aRowsNumber; i++) {
                    j0 = aColumnOffset;
                    for (j = 0; j < aColumnsNumber; j++) {
                        expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / apLocalTheta[1];
                        if (expr == 0) {
                            apMatrixA[i + j * aRowsNumber] = sigma_square /*+ 1e-4*/;
                        } else {
                            apMatrixA[i + j * aRowsNumber] = con * pow(expr, apLocalTheta[2])
                                                             * gsl_sf_bessel_Knu(apLocalTheta[2], expr); // Matern Function
                        }

                        j0++;
                    }
                    i0++;
                }
            };

            static Kernel* create() {
                return new UnivariateMaternStationary();
            }

        };
        EXAGEOSTAT_REGISTER_PLUGIN(UnivariateMaternStationary, UnivariateMaternStationary::create);
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNSTATIONARY_HPP
