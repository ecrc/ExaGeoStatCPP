
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Definitions.hpp
 * @version 1.0.0
 * @brief This file contains common definitions used in ExaGeoStat software package.
 * @details These definitions include enums for dimension, computation, precision, and floating point arithmetic;
 * A macro for instantiating template classes with supported types; and a set of available kernels.
 * @author Sameh Abdulah
 * @date 2023-03-21
**/


#ifndef EXAGEOSTATCPP_DEFINITIONS_HPP
#define EXAGEOSTATCPP_DEFINITIONS_HPP

#include <iostream>
#include <set>

/**
 * @def EXAGEOSTAT_INSTANTIATE_CLASS
 * @brief Macro definition to instantiate the EXAGEOSTAT template classes with supported types.
**/
#define EXAGEOSTAT_INSTANTIATE_CLASS(TEMPLATE_CLASS)   template class TEMPLATE_CLASS<float>;  \
                                                    template class TEMPLATE_CLASS<double>;

// Variables sizes.
#define SIZE_OF_FLOAT 4
#define SIZE_OF_DOUBLE 8

namespace exageostat {
    namespace common {
        /**
         * @enum Dimension
         * @brief Enum denoting the dimension of generated data.
         */
        enum Dimension {
            Dimension2D = 0,
            Dimension3D = 1,
            DimensionST = 2,
        };

        /**
         * @enum Computation
         * @brief Enum denoting the types of computations that can be requested,
         * to use the required Linear Algebra solver library.
         */
        enum Computation {
            EXACT_DENSE = 0,
            DIAGONAL_APPROX = 1,
            TILE_LOW_RANK = 2,
        };

        /**
         * @enum Precision
         * @brief Enum denoting the precision of operations that are supported to be done on the matrix.
         */
        enum Precision {
            SINGLE = 0,
            DOUBLE = 1,
            MIXED = 2,
        };

        /**
         * @enum FloatPoint
         * @brief Enum denoting the floating point arithmetic of the matrix.
         */
        enum FloatPoint : int {
            EXAGEOSTAT_BYTE = 0,
            EXAGEOSTAT_INTEGER = 1,
            EXAGEOSTAT_REAL_FLOAT = 2,
            EXAGEOSTAT_REAL_DOUBLE = 3,
            EXAGEOSTAT_COMPLEX_FLOAT = 4,
            EXAGEOSTAT_COMPLEX_DOUBLE = 5,
        };


        /// TODO: Make it automatically generated from the kernels files names
        /**
         * @var availableKernels
         * @brief Set denoting the available kernels supported in matrix generation.
         * @details This set should be updated manually to add new kernels.
         */
        const std::set<std::string> availableKernels = {"univariate_matern_stationary",
                                                        "univariate_matern_non_stationary",
                                                        "bivariate_matern_flexible",
                                                        "bivariate_matern_parsimonious",
                                                        "bivariate_matern_parsimonious_profile",
                                                        "univariate_matern_nuggets_stationary",
                                                        "univariate_spacetime_matern_stationary",
                                                        "univariate_matern_dsigma_square",
                                                        "univariate_matern_dnu",
                                                        "univariate_matern_dbeta",
                                                        "univariate_matern_ddsigma_square",
                                                        "univariate_matern_ddsigma_square_beta",
                                                        "univariate_matern_ddsigma_square_nu",
                                                        "univariate_matern_ddbeta_beta",
                                                        "univariate_matern_ddbeta_nu",
                                                        "univariate_matern_ddnu_nu",
                                                        "bivariate_spacetime_matern_stationary",
                                                        "univariate_matern_non_gaussian",
                                                        "univariate_exp_non_gaussian",
                                                        "bivariate_spacetime_matern_stationary",
                                                        "trivariate_matern_parsimonious",
                                                        "trivariate_matern_parsimonious_profile",
                                                        "univariate_matern_non_stat"
        };
    }//namespace common
}//namespace exageostat

#endif //EXAGEOSTATCPP_DEFINITIONS_HPP
