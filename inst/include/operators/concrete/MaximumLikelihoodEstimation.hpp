// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file MaximumLikelihoodEstimation.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-30
**/

#ifndef EXAGEOSTATCPP_MAXIMUMLIKELIHOODESTIMATION_HPP
#define EXAGEOSTATCPP_MAXIMUMLIKELIHOODESTIMATION_HPP

#include <operators/OperatorMethods.hpp>
namespace exageostat {
    namespace operators {
        /**
         * @brief
         * ChameleonImplementationDense is a concrete implementation of LinearAlgebraMethods class for dense matrices.
         * @tparam T Type of matrix elements.
         */
        template<typename T>
        class MaximumLikelihoodEstimation : public OperatorMethods<T>{

        public:

        };

        EXAGEOSTAT_INSTANTIATE_CLASS(MaximumLikelihoodEstimation)

    }//namespace operators
}//namespace exageostat
#endif //EXAGEOSTATCPP_MAXIMUMLIKELIHOODESTIMATION_HPP
