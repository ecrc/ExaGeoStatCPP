// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestBivariateSpacetimeMaternStationary.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <libraries/catch/catch.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <data-generators/DataGenerator.hpp>
#include <vector>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::generators;
using namespace std;

void TEST_KERNEL_GENERATION_BivariateSpacetimeMaternStationary() {

    // Create a unique pointer to a DataGenerator object
    std::unique_ptr<DataGenerator> synthetic_generator;

    // Create a new synthetic_data_configurations object with the provided command line arguments
    SyntheticDataConfigurations synthetic_data_configurations;

    synthetic_data_configurations.SetProblemSize(9);
    synthetic_data_configurations.SetKernel("BivariateSpacetimeMaternStationary");
#ifdef EXAGEOSTAT_USE_CHAMELEON
    synthetic_data_configurations.SetDenseTileSize(5);
    synthetic_data_configurations.SetComputation(EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    synthetic_data_configurations.SetLowTileSize(5);
    synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
#endif
    synthetic_data_configurations.SetDimension(DimensionST);
    synthetic_data_configurations.SetIsSynthetic(true);
    synthetic_data_configurations.SetPrecision(DOUBLE);
    synthetic_data_configurations.SetTimeSlot(5);

    vector<double> target_theta{-1, -1, -1, -1, -1, -1, -1, -1, -1};
    synthetic_data_configurations.SetTargetTheta(target_theta);

    vector<double> lb{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    synthetic_data_configurations.SetLowerBounds(lb);

    vector<double> ub{5, 5, 5, 5, 5, 5, 5, 5};
    synthetic_data_configurations.SetUpperBounds(ub);

    vector<double> initial_theta{1, 1, 1, 0.1, 0.5, 1, 1.5, 0.1, 0.1};
    synthetic_data_configurations.SetInitialTheta(initial_theta);

    // Create the DataGenerator object
    synthetic_generator = synthetic_generator->CreateGenerator(&synthetic_data_configurations);

    // Initialize the locations of the generated data
    synthetic_generator->GenerateLocations();

    // Set the locations with these values.
    vector<double> x = {0.257383, 0.24216, 0.276439, 0.456059, 0.44074, 0.493973, 0.797276, 0.953928, 0.869523, 0.257383, 0.24216, 0.276439, 0.456059, 0.44074, 0.493973, 0.797276, 0.953928, 0.869523, 0.257383, 0.24216, 0.276439, 0.456059, 0.44074, 0.493973, 0.797276, 0.953928, 0.869523, 0.257383, 0.24216, 0.276439, 0.456059, 0.44074, 0.493973, 0.797276, 0.953928, 0.869523, 0.257383, 0.24216, 0.276439, 0.456059, 0.44074, 0.493973, 0.797276, 0.953928, 0.869523};
    vector<double> y = {0.138502, 0.579584, 0.75268, 0.238195, 0.514392, 0.867699, 0.17024, 0.610985, 0.891279, 0.138502, 0.579584, 0.75268, 0.238195, 0.514392, 0.867699, 0.17024, 0.610985, 0.891279, 0.138502, 0.579584, 0.75268, 0.238195, 0.514392, 0.867699, 0.17024, 0.610985, 0.891279, 0.138502, 0.579584, 0.75268, 0.238195, 0.514392, 0.867699, 0.17024, 0.610985, 0.891279, 0.138502, 0.579584, 0.75268, 0.238195, 0.514392, 0.867699, 0.17024, 0.610985, 0.891279};
    vector<int> z = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5};

    for (auto i = 0; i < x.size(); i++) {
        synthetic_generator->GetLocations()->GetLocationX()[i] = x[i];
        synthetic_generator->GetLocations()->GetLocationY()[i] = y[i];
        synthetic_generator->GetLocations()->GetLocationZ()[i] = z[i];
    }

    synthetic_generator->GenerateDescriptors();

    auto descriptorC = synthetic_data_configurations.GetDescriptorC()[0];

    exageostat::dataunits::Locations *l1 = synthetic_generator->GetLocations();

    auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(
            synthetic_data_configurations.GetComputation());
    linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);
    linearAlgebraSolver->CovarianceMatrixCodelet(descriptorC, EXAGEOSTAT_LOWER, l1, l1, nullptr,
                                                 synthetic_data_configurations.GetInitialTheta(), 0,
                                                 synthetic_generator->GetKernel());
    auto *A = linearAlgebraSolver->GetMatrix();

    // Define the expected output
    double expected_output_data[] = {1.000000, 0.745356, 0.788999, 0.738211, 0.000000, 0.745356, 1.000000, 0.738211,
                                     0.999532, 0.000000, 0.788999, 0.738211, 1.000000, 0.745356, 0.000000, 0.738211,
                                     0.999532, 0.745356, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                     0.000000};
    int m = 6;
    int n = 4;

    for (int i = 0; i < m * n; i++) {
        double diff = A[i] - expected_output_data[i];
        REQUIRE(diff == Approx(0.0).margin(1e-6));
    }
}

TEST_CASE("BivariateSpacetimeMaternStationary kernel test") {
    TEST_KERNEL_GENERATION_BivariateSpacetimeMaternStationary();
}