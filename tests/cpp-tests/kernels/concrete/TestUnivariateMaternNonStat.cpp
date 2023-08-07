// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternNonStat.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <libraries/catch/catch.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::common;

using namespace std;

void TEST_KERNEL_GENERATION_UnivariateMaternNonStat() {

    SECTION("UnivariateMaternNonStat") {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 9;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternNonStat");
        synthetic_data_configurations.SetDimension(Dimension2D);

        vector<double> lb{0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5, 5, 5, 5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{0.04, 1.57, 0.33, -1, 0.8, 0.1, -0.5, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 5;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(EXACT_DENSE, 3, 0);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
        synthetic_data_configurations.SetLowTileSize(5);
        synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(TILE_LOW_RANK, 4, 0);
#endif
        // Initialize the seed manually with zero, to get the first generated seeded numbers.
        srand(0);

        auto *data = ExaGeoStat<double>::ExaGeoStatGenerateData(&synthetic_data_configurations);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                          exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
        auto *A = data->GetDescriptorData()->GetDescriptor(HiCMA_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
#endif

        // Define the expected output
        //// TODO: FIX VALUES
        double expected_output_data[] = {0.842571, 0.368249, 0.087037, 0.120736,
                                         0.368249, 0.782441, 0.165305, 0.265572,
                                         0.087037, 0.165305, 0.755169, 0.341821,
                                         0.120736, 0.265572, 0.341821, 0.733985};

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            cout << A[i] << endl;
//            REQUIRE(diff == Approx(0.0).margin(1e-6));
        }

#ifdef EXAGEOSTAT_USE_CHAMELEON
        // Finalize ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(EXACT_DENSE, data->GetDescriptorData());
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
        // Finalize ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(TILE_LOW_RANK, data->GetDescriptorData());
#endif

    }
}

TEST_CASE("UnivariateMaternNonStat kernel test") {
    TEST_KERNEL_GENERATION_UnivariateMaternNonStat();
}