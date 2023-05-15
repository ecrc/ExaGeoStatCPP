
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareBeta.hpp
 * @brief Defines the UnivariateMaternDdsigmaSquareBeta class, a Univariate Matern Ddsigma Square Beta kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUAREBETA_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUAREBETA_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternDdsigmaSquareBeta
         * @brief A class representing a Univariate Matern Ddsigma Square Beta kernel.
         *
         * This class represents a Univariate Matern Ddsigma Square Beta kernel, which is a subclass of the Kernel class. It provides
         * a method for generating a covariance matrix using a set of input locations and kernel parameters.
         */
        class UnivariateMaternDdsigmaSquareBeta : public Kernel {

        public:

            /**
             * @brief Constructs a new UnivariateMaternDdsigmaSquareBeta object.
             *
             * Initializes a new UnivariateMaternDdsigmaSquareBeta object with default values.
             */
            UnivariateMaternDdsigmaSquareBeta();

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
             * @param[in] aLocalTheta An array of kernel parameters.
             * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
             */
            void GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber, int aRowOffset,
                                          int aColumnOffset, dataunits::Locations *apLocation1,
                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                          double *aLocalTheta, int aDistanceMetric) override ;
            /**
             * @brief Creates a new UnivariateMaternDdsigmaSquareBeta object.
             * @return A pointer to the new UnivariateMaternDdsigmaSquareBeta object.
             *
             * This method creates a new UnivariateMaternDdsigmaSquareBeta object and returns a pointer to it.
             */
            static Kernel *Create();

        private:
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat
#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUAREBETA_HPP
