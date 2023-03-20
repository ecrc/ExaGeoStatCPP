
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors_main.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/dense/ChameleonAllocateDescriptors.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::linearAlgebra;
using namespace exageostat::configurations::data_configurations;
using namespace std;

int main(int argc, char **argv) {

    unique_ptr<LinearAlgebraFactory> chameleon;

    // Object has automatic storage duration (usually is on the stack)
    auto syntheticDataConfigurations = new SyntheticDataConfigurations(argc, argv);

//    chameleon->CreateLinearAlgebraSolver(syntheticDataConfigurations->GetComputation());
    chameleon->CreateLinearAlgebraSolver();

    return 0;
}