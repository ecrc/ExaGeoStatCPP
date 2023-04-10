
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors_main.cpp
 * @brief This file demonstrates how to use the LinearAlgebraFactory class from the ExaGeoStat software package to create
 * and initialize linear algebra solvers for different precision types.
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


/**
 * @brief The main function of the program.
 *
 * This function demonstrates how to use the LinearAlgebraFactory class from the ExaGeoStat software package
 * to create and initialize linear algebra solvers for different precision types.
 *
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return The status code of the program.
 */
int main(int argc, char **argv) {

    // Create an instance of the SyntheticDataConfigurations class with user-defined configurations.
    auto syntheticDataConfigurations = new SyntheticDataConfigurations(argc, argv);

    // Create and initialize linear algebra solvers for different precision types.
    if (syntheticDataConfigurations->GetPrecision() == SINGLE) {
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(
                syntheticDataConfigurations->GetComputation());
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
        linearAlgebraSolver->InitiateDescriptors();
    } else if (syntheticDataConfigurations->GetPrecision() == DOUBLE) {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(
                syntheticDataConfigurations->GetComputation());
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
        linearAlgebraSolver->InitiateDescriptors();
    } else if (syntheticDataConfigurations->GetPrecision() == MIXED) {
        // TODO: Add implementation for mixed-precision linear algebra solver.
    }
    // Clean up allocated memory
    delete syntheticDataConfigurations;
    return 0;
}