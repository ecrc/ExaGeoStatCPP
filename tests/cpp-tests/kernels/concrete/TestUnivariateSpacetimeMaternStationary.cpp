// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateSpacetimeMaternStationary.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-10
**/
#include <libraries/catch/catch.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <data-generators/DataGenerator.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::generators;
using namespace std;

void TEST_KERNEL_GENERATION_UnivariateSpacetimeMaternStationary() {

    SECTION("UnivariateSpacetimeMaternStationary") {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        SyntheticDataConfigurations synthetic_data_configurations;

        synthetic_data_configurations.SetProblemSize(9);
        synthetic_data_configurations.SetKernel("UnivariateSpacetimeMaternStationary");
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

        // Create a unique pointer to a DataGenerator object
        std::unique_ptr<DataGenerator<double>> synthetic_generator;

        vector<double> lb{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5, 5, 5, 5, 0};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{1, 1, 0.1, 0.5, 0.5, 0.1, 0};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        vector<double> target_theta{-1, -1, -1, -1, -1, -1, 0};
        synthetic_data_configurations.SetTargetTheta(target_theta);

        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        // Create the DataGenerator object
        synthetic_generator = synthetic_generator->CreateGenerator(&synthetic_data_configurations);

        // Initialize the seed manually with zero, to get the first generated seeded numbers.
        srand(0);
        // Generated locations data
        synthetic_generator->GenerateLocations();
        synthetic_generator->GenerateDescriptors();

        auto descriptorC = synthetic_data_configurations.GetDescriptorC()[0];
        exageostat::dataunits::Locations *l1 = synthetic_generator->GetLocations();

        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(
                synthetic_data_configurations.GetComputation());
        linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);
        linearAlgebraSolver->CovarianceMatrixCodelet(descriptorC, EXAGEOSTAT_LOWER, l1, l1, nullptr,
                                                     synthetic_data_configurations.GetInitialTheta().data(), 0,
                                                     synthetic_generator->GetKernel());

        auto *A = linearAlgebraSolver->GetMatrix();

        // Define the expected output.
        double expected_output_data[] = {1.000000, 0.216541, 0.171440, 0.164316, 0.152318,
                                         0.216541, 1.000000, 0.115157, 0.163143, 0.209362,
                                         0.171440, 0.115157, 1.000000, 0.179330, 0.121224,
                                         0.164316, 0.163143, 0.179330, 1.000000, 0.244122,
                                         0.152318, 0.209362, 0.121224, 0.244122, 1.000000};
        size_t m = 4;
        size_t n = 6;
        for (size_t i = 0; i < m * n; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Approx(0.0).margin(1e-6));
        }
        synthetic_generator->DestoryDescriptors();

        // Finalize ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);
        delete linearAlgebraSolver;

    }
}
TEST_CASE("UnivariateSpacetimeMaternStationary kernel test") {
    TEST_KERNEL_GENERATION_UnivariateSpacetimeMaternStationary();
}