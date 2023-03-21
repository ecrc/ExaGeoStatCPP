
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

using namespace exageostat::linearAlgebra;
using namespace exageostat::configurations;

template<typename T>
void LinearAlgebraMethods<T>::SetConfigurations(Configurations *apConfigurations) {
    this->mpConfigurations = apConfigurations;
}