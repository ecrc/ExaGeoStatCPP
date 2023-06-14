// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file OperatorFactory.cpp
 * @brief Implementation file for the OperatorFactory class, which creates operators based on the input operator type.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-30
 *
 * This file contains the implementation of the CreateOperator() method of the OperatorFactory class, which creates operators based on the input operator type.
 * The OperatorFactory class is templated on the data type of the operator, and provides a static CreateOperator() method for creating new operator instances.
 *
 * The implementation of the CreateOperator() method currently supports only the MLE operator type.
 *
**/

#include <operators/OperatorFactory.hpp>
#include <operators/concrete/MaximumLikelihoodEstimation.hpp>

using namespace exageostat::operators;
using namespace exageostat::common;

template<typename T>
std::unique_ptr<OperatorMethods<T>> OperatorFactory<T>::CreateOperator(Operators aOperators) {

    if (aOperators == MLE) {
        std::make_unique<MaximumLikelihoodEstimation<T>>();
    }
    throw std::runtime_error("You need to select operator");
}
