
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestTrendModel.cpp
 * @brief Unit tests for the TestTrendModel kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestTrendModel kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-11-11
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

using namespace std;

using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::configurations;

void TEST_KERNEL_GENERATION_TrendModel() {

    //TODO: check for correct values to generate positive matrix
    SECTION("TrendModel")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;

        int N = 2;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("TrendModel");
        synthetic_data_configurations.SetDimension(Dimension2D);

        vector<double> lb{0.1,0.1,0.1};
        vector<double> ub{5,5,5};
        vector<double> initial_theta{1,0.1,0.5};
        synthetic_data_configurations.SetLowerBounds(lb);
        synthetic_data_configurations.SetUpperBounds(ub);
        synthetic_data_configurations.SetInitialTheta(initial_theta);
        synthetic_data_configurations.SetEstimatedTheta({1,0.1,0.5});

        int dts = 1;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        std::unique_ptr<ExaGeoStatData<double>> data;
//        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(synthetic_data_configurations, data);
//        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
//                                                                          exageostat::common::DESCRIPTOR_Z).chameleon_desc;
//        auto *A = (double *) CHAM_descriptorZ->mat;
//         Define the expected output
        double expected_output_data[] = {-1.272336, -2.475473, 0.545850, -0.120985, 0.242569, -1.544215, 0.098647,
                                         0.779835, -1.481391};

        for (size_t i = 0; i < N; i++) {
//            double diff = A[i] - expected_output_data[i];
//            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("TrendModel kernel test") {
    TEST_KERNEL_GENERATION_TrendModel();
}