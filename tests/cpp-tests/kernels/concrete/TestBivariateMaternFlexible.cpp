// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestBivariateMaternFlexible.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-09
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

void TEST_KERNEL_GENERATION_BivariateMaternFlexible() {

    SECTION("BivariateMaternFlexible")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        SyntheticDataConfigurations synthetic_data_configurations;

        synthetic_data_configurations.SetProblemSize(9);
        synthetic_data_configurations.SetKernel("BivariateMaternFlexible");
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

        vector<double> target_theta{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
        synthetic_data_configurations.SetTargetTheta(target_theta);

        vector<double> lb{0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{0.3, 0.6, 0.01, 0.3, 0.9, 0.9, 0.05, 0.3, 1.5, 0.9, 0.99};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        // Create the DataGenerator object
        auto synthetic_generator = DataGenerator<double>::CreateGenerator(&synthetic_data_configurations);

        // Initialize the seed manually with zero, to get the first generated seeded numbers.
        srand(0);
        // Generated locations data
        synthetic_generator->GenerateLocations();
        synthetic_generator->GenerateDescriptors();

        auto descriptorC = synthetic_data_configurations.GetDescriptorC()[0];
        exageostat::dataunits::Locations *l1 = synthetic_generator->GetLocations();
        int upper_lower = EXAGEOSTAT_LOWER;
        int distance_metric = 0;
        synthetic_generator->GetLinearAlgberaSolver()->CovarianceMatrixCodelet(descriptorC, upper_lower, l1, l1,
                                                                               nullptr,
                                                                               synthetic_data_configurations.GetInitialTheta().data(),
                                                                               distance_metric,
                                                                               synthetic_generator->GetKernel());

        auto *A = synthetic_generator->GetLinearAlgberaSolver()->GetMatrix();

        // Define the expected output
        double expected_output_data[] = {0.300000, 0.217899, 0.140362,
                                         0.148157, 0.600000, 0.148157,
                                         0.264357, 0.148157, 0.300000};
        int m = 3;
        int n = 3;
        for (int i = 0; i < m * n; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Approx(0.0).margin(1e-6));
        }

        // Finalize ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);
    }
}

TEST_CASE("Bivariate Matern Flexible kernel test") {
    TEST_KERNEL_GENERATION_BivariateMaternFlexible();
}