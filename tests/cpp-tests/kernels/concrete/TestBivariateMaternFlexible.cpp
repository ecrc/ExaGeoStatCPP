
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestBivariateMaternFlexible.cpp
 * @brief Unit tests for the BivariateMaternFlexible kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the BivariateMaternFlexible kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-09
**/

#include <libraries/catch/catch.hpp>
#include <api/ExaGeoStat.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::hardware;

void TEST_KERNEL_GENERATION_BivariateMaternFlexible() {

    SECTION("BivariateMaternFlexible")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 27;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("BivariateMaternFlexible");
        synthetic_data_configurations.SetDimension(Dimension2D);
        synthetic_data_configurations.SetCoresNumber(4);

        vector<double> target_theta{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
        synthetic_data_configurations.SetTargetTheta(target_theta);

        vector<double> lb{0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{0.3, 0.6, 0.01, 0.3, 0.9, 0.9, 0.05, 0.3, 1.5, 0.9, 0.99};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 16;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
                exageostat::dataunits::ExaGeoStatData<double> data(synthetic_data_configurations.GetProblemSize(), synthetic_data_configurations.GetDimension(), hardware);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, synthetic_data_configurations, data);
        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {-0.696887, -2.069051, -0.417382, -1.001165, -0.235721, -1.726871, -0.285059,
                                         -0.516932, -1.100053, -1.760795, -1.430539, -1.356296, -0.274180, -0.753934,
                                         -1.039373, -1.076400, -0.099143, -1.549572, 0.778055, -0.730216, 0.706544,
                                         -0.804701, 0.315835, -0.136930, -0.089561, -1.776487, -0.101322, -1.548683,
                                         0.245666, -2.050699, -1.080789, -2.173827, -0.816287, -2.922996, -0.933723,
                                         -2.220099, -0.380451, -1.029619, -0.376941, -1.691071, -0.200486, -2.772512,
                                         -0.234095, -1.313128, -0.647894, -0.937577, 0.012360, -2.657097, -0.401142,
                                         -1.074726, -0.245235, -1.372391, 0.096746, -1.470435};

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Approx(0.0).margin(1e-6));
        }

#endif
    }
}

TEST_CASE("Bivariate Matern Flexible kernel test") {
TEST_KERNEL_GENERATION_BivariateMaternFlexible();

}