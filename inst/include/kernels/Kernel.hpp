
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

namespace exageostat {
    namespace Kernels {

        /**
         * @class Kernels
         * @brief A class containing the main kernel functions.
         */
        class Kernel{
        public:

            /**
             * @brief Generate covariance matrix A.
             *
             * @param[out]
             *
             * @param[in] m
             *          The number of rows in the tile A.
             *
             * @param[in] n
             *          The number of cols in the tile A.
             *
             * @param[in] m0
             *          Global row index of the tile A.
             *
             * @param[in] n0
             *          Global col index of the tile A.
             *
             * @param[in] l1
             *          Location struct of the first input.
             *
             * @param[in] l2
             *          Location struct of the second input.
             *
             * @param[in] localTheta
             *          Parameter vector that is used to generate the output covariance matrix.
             *
             * @param[in] distance_metric
             *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"             *
             * @return A
             *          The m-by-n matrix on which to compute the covariance matrix.
             */
            virtual void GenerateCovarianceMatrix() = 0;
            virtual void InitializeMatrixParameters() = 0;

        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_KERNELS_HPP