// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestTrivariateMaternParsimonious.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <libraries/catch/catch.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <data-generators/DataGenerator.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::generators;
using namespace std;

void TEST_KERNEL_GENERATION_TrivariateMaternParsimonious() {

    // Create a unique pointer to a DataGenerator object
    std::unique_ptr<DataGenerator> synthetic_generator;

    // Create a new synthetic_data_configurations object with the provided command line arguments
    SyntheticDataConfigurations synthetic_data_configurations;

    synthetic_data_configurations.SetProblemSize(9);
    synthetic_data_configurations.SetKernel("TrivariateMaternParsimonious");
#ifdef EXAGEOSTAT_USE_CHAMELEON
    synthetic_data_configurations.SetDenseTileSize(5);
    synthetic_data_configurations.SetComputation(EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    synthetic_data_configurations.SetLowTileSize(5);
    synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
#endif
    synthetic_data_configurations.SetDimension(Dimension2D);
    synthetic_data_configurations.SetIsSynthetic(true);
    synthetic_data_configurations.SetPrecision(DOUBLE);

    vector<double> lb{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    synthetic_data_configurations.SetLowerBounds(lb);

    vector<double> ub{5, 5, 5, 5, 1, 1, 5, 5, 5, 5};
    synthetic_data_configurations.SetUpperBounds(ub);

    vector<double> initial_theta{1, 1, 1, 0.1, 0.5, 1, 1.5, 0.1, 0.1, 0};
    synthetic_data_configurations.SetInitialTheta(initial_theta);

    // Create the DataGenerator object
    synthetic_generator = synthetic_generator->CreateGenerator(&synthetic_data_configurations);

    // Initialize the locations of the generated data
    synthetic_generator->GenerateLocations();

    // Set the locations with these values.
    vector<double> x = {0.257389, 0.456062, 0.797269, 0.242161, 0.440742, 0.276432, 0.493965, 0.953933, 0.86952};
    vector<double> y = {0.138506, 0.238193, 0.170245, 0.579583, 0.514397, 0.752682, 0.867704, 0.610986, 0.891279};

    for (auto i = 0; i < x.size(); i++) {
        synthetic_generator->GetLocations()->GetLocationX()[i] = x[i];
        synthetic_generator->GetLocations()->GetLocationY()[i] = y[i];
    }

    synthetic_generator->GenerateDescriptors();

    auto descriptorC = synthetic_data_configurations.GetDescriptorC()[0];

    exageostat::dataunits::Locations *l1 = synthetic_generator->GetLocations();

    auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(
            synthetic_data_configurations.GetComputation());
    linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);
    linearAlgebraSolver->GenerateObservationsVector(descriptorC, l1, l1, nullptr, synthetic_data_configurations.GetInitialTheta(), 0,
                                                    synthetic_generator->GetKernel());
    auto *A = linearAlgebraSolver->GetMatrix();
    // Define the expected output
    double expected_output_data[] = {1.000000, 0.094281, 0.086603, 1.000000};


    size_t m = 4;
    size_t n = 1;
    for (size_t i = 0; i < m * n; i++) {
        double diff = A[i] - expected_output_data[i];
        REQUIRE(diff == Approx(0.0).margin(1e-6));
    }
}

TEST_CASE("TrivariateMaternParsimonious kernel test") {
    TEST_KERNEL_GENERATION_TrivariateMaternParsimonious();
}