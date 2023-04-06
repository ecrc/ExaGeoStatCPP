
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
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace exageostat::linearAlgebra;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::common;
using namespace std;

int main(int argc, char **argv) {

    // Object has automatic storage duration (usually is on the stack)
    auto* syntheticDataConfigurations = new SyntheticDataConfigurations(argc, argv);
    if (syntheticDataConfigurations->GetPrecision() == SINGLE) {
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(syntheticDataConfigurations->GetComputation());
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
        linearAlgebraSolver->InitiateDescriptors();
    } else if (syntheticDataConfigurations->GetPrecision() == DOUBLE) {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(syntheticDataConfigurations->GetComputation());
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
        linearAlgebraSolver->InitiateDescriptors();

    } else if (syntheticDataConfigurations->GetPrecision() == MIXED) {

    }
    return 0;
}