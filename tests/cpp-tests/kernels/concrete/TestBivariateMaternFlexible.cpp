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
#include <api/ExaGeoStat.hpp>

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::api;

using namespace std;

void TEST_KERNEL_GENERATION_BivariateMaternFlexible() {

    SECTION("BivariateMaternFlexible")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 9;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("BivariateMaternFlexible");
        synthetic_data_configurations.SetDimension(Dimension2D);

        vector<double> target_theta{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
        synthetic_data_configurations.SetTargetTheta(target_theta);

        vector<double> lb{0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{0.3, 0.6, 0.01, 0.3, 0.9, 0.9, 0.05, 0.3, 1.5, 0.9, 0.99};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

//#ifdef EXAGEOSTAT_USE_CHAMELEON
//        int dts = 5;
//        synthetic_data_configurations.SetDenseTileSize(dts);
//        synthetic_data_configurations.SetComputation(EXACT_DENSE);
//        // Initialise ExaGeoStat Hardware.
//        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(EXACT_DENSE, 3, 0);
//#endif
//#ifdef EXAGEOSTAT_USE_HiCMA
//        synthetic_data_configurations.SetLowTileSize(5);
//        synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
//#endif
//
//        srand(0);
//        auto *data = ExaGeoStat<double>::ExaGeoStatGenerateData(&synthetic_data_configurations);
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
//                                                                          exageostat::common::DESCRIPTOR_Z).chameleon_desc;
//        auto *A = (double *) CHAM_descriptorZ->mat;
//#endif
//#ifdef EXAGEOSTAT_USE_HiCMA
//        auto *A = data->GetDescriptorData()->GetDescriptor(HiCMA_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
//#endif
//
//        // Define the expected output
//        double expected_output_data[] = {0.300000, 0.217899, 0.140362,
//                                         0.148157, 0.600000, 0.148157,
//                                         0.264357, 0.148157, 0.300000};
//
//        for (int i = 0; i < N; i++) {
//            double diff = A[i] - expected_output_data[i];
//            REQUIRE(diff == Approx(0.0).margin(1e-6));
//        }
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//        // Finalize ExaGeoStat Hardware.
//        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(EXACT_DENSE, data->GetDescriptorData());
//#endif
//#ifdef EXAGEOSTAT_USE_HiCMA
//        // Finalize ExaGeoStat Hardware.
//        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(TILE_LOW_RANK, data->GetDescriptorData());
//#endif
    }
}

TEST_CASE("Bivariate Matern Flexible kernel test") {
    TEST_KERNEL_GENERATION_BivariateMaternFlexible();
}