// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternDdsigmaSquare.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <libraries/catch/catch.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::common;
using namespace exageostat::api;

using namespace std;

void TEST_KERNEL_GENERATION_UnivariateMaternDdsigmaSquare() {

    SECTION("UnivariateMaternDdsigmaSquare")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 9;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternDdsigmaSquare");
        synthetic_data_configurations.SetDimension(Dimension2D);

        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 5;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(EXACT_DENSE, 3, 0);

#endif
#ifdef EXAGEOSTAT_USE_HiCMA
        int lts = 5;
        synthetic_data_configurations.SetLowTileSize(lts);
        synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(TILE_LOW_RANK, 3, 0);
#endif


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

        for (int i = 0; i < N; i++) {
            double diff = A[i] - 0;
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

TEST_CASE("UnivariateMaternDdsigmaSquare kernel test") {
    TEST_KERNEL_GENERATION_UnivariateMaternDdsigmaSquare();
}