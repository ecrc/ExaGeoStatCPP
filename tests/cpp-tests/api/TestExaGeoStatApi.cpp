
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatApi.cpp
 * @brief Test suite for the ExaGeoStat API's data generation functionality.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-08-07
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;
using namespace exageostat::hardware;

void TEST_GENERATE_DATA() {
    SECTION("Data generation - Observations")
    {
        int seed = 0;
        srand(seed);

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 16;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternStationary");
#ifdef EXAGEOSTAT_USE_CHAMELEON
        synthetic_data_configurations.SetDenseTileSize(9);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
        synthetic_data_configurations.SetLowTileSize(5);
        synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
#endif

        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{1, 0.1, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 4, 0); // Or you could use configurations.GetComputation().
        exageostat::dataunits::ExaGeoStatData<double> data(synthetic_data_configurations.GetProblemSize(),
                                                           synthetic_data_configurations.GetDimension(), hardware);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, synthetic_data_configurations, data);

        // Define the expected output for desk Z
        double expected_output_data[] = {-1.272336, -2.590700, 0.512143, -0.163880, 0.313504, -1.474411, 0.161705,
                                         0.623389, -1.341858, -1.054282, -1.669383, 0.219171, 0.971214, 0.538973,
                                         -0.752828, 0.290822};

        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                         DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        double diff;

        for (int i = 0; i < N; i++) {
            diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

void TEST_MODEL_DATA() {
    int seed = 0;
    srand(seed);

    // Create a new configurations object with the provided command line arguments
    Configurations configurations;
    int N = 16;
    configurations.SetProblemSize(N);
    configurations.SetKernelName("UnivariateMaternStationary");
    int dts = 8;

    configurations.SetDenseTileSize(dts);
    configurations.SetComputation(EXACT_DENSE);

    configurations.SetPrecision(DOUBLE);

    configurations.SetMaxMleIterations(3);
    configurations.SetTolerance(pow(10, -4));

    vector<double> lb{0.1, 0.1, 0.1};
    configurations.SetLowerBounds(lb);
    configurations.SetStartingTheta(lb);

    vector<double> ub{5, 5, 5};
    configurations.SetUpperBounds(ub);

    vector<double> initial_theta{1, 0.1, 0.5};
    configurations.SetInitialTheta(initial_theta);


    SECTION("Data Modeling")
    {

        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 4, 0); // Or you could use configurations.GetComputation().
        exageostat::dataunits::ExaGeoStatData<double> data(configurations.GetProblemSize(),
                                                           configurations.GetDimension(), hardware);

        FloatPoint float_point = EXAGEOSTAT_REAL_DOUBLE;

        //initiating the matrix of the CHAMELEON Descriptor Z.
        auto *z_matrix = (double *) malloc(sizeof(double) * N);
        // Assign z_values to z_matrix
        double z_values[] = {-1.272336140360187606, -2.590699695867695773, 0.512142584178685967, -0.163880452049749520,
                             0.313503633252489700, -1.474410682226017677, 0.161705025505231914, 0.623389205185149065,
                             -1.341858445399783495, -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                             0.971213790000161170, 0.538973474182433021, -0.752828466476077041, 0.290822066007430102};

        for (int i = 0; i < N; ++i) {
            z_matrix[i] = z_values[i];
        }
        data.GetDescriptorData()->SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z, configurations.GetIsOOC(), z_matrix,
                                                float_point, dts, dts,
                                                dts * dts, N, 1, 0, 0, N, 1, configurations.GetPGrid(),
                                                configurations.GetQGrid());

        //initiating the matrix of the CHAMELEON Descriptor C.
        auto *c_matrix = (double *) malloc(sizeof(double) * N * N);
        for (int i = 0; i < N * N; i++) {
            c_matrix[i] = 0;
        }
        data.GetDescriptorData()->SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C, configurations.GetIsOOC(), c_matrix,
                                                float_point, dts, dts,
                                                dts * dts, N, N, 0, 0, N, N, configurations.GetPGrid(),
                                                configurations.GetQGrid());

        //Initiating the descriptors needed by data modeling.
        data.GetDescriptorData()->SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, configurations.GetIsOOC(),
                                                nullptr, float_point, dts,
                                                dts, dts * dts, N, 1, 0, 0, N, 1, configurations.GetPGrid(),
                                                configurations.GetQGrid());

        data.GetDescriptorData()->SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT, configurations.GetIsOOC(),
                                                nullptr,
                                                float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
                                                configurations.GetPGrid(), configurations.GetQGrid());

        data.GetDescriptorData()->SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT, configurations.GetIsOOC(),
                                                nullptr, float_point,
                                                dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, configurations.GetPGrid(),
                                                configurations.GetQGrid());

        //creating locations x and y.
        double location_x[] = {0.193041886015106440, 0.330556191348134576, 0.181612878614480805, 0.370473792629892440,
                               0.652140077821011688, 0.806332494087129037, 0.553322652018005678, 0.800961318379491916,
                               0.207324330510414295, 0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                               0.528267338063630132, 0.974792095826657490, 0.552452887769893985, 0.877592126344701295};
        for (int i = 0; i < N; ++i) {
            data.GetLocations()->GetLocationX()[i] = location_x[i];
        }


        double location_y[] = {0.103883421072709245, 0.135790035858701447, 0.434683756771190977, 0.400778210116731537,
                               0.168459601739528508, 0.105195696955825133, 0.396398870832379624, 0.296757457846952011,
                               0.564507515068284116, 0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                               0.573571374074921758, 0.568657969024185528, 0.935835812924391552, 0.942824444953078489};
        for (int i = 0; i < N; ++i) {
            data.GetLocations()->GetLocationY()[i] = location_y[i];
        }

        //TODO: the log_likelihood value needs to be saved somewhere for it to be tested.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data);
        //TODO: check for thr log_likelihood values.

        free(z_matrix);
        free(c_matrix);

    }
    SECTION("Data Generation and Modeling")
    {
        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 4, 0); // Or you could use configurations.GetComputation().
        exageostat::dataunits::ExaGeoStatData<double> data(configurations.GetProblemSize(),
                                                           configurations.GetDimension(), hardware);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, configurations, data);
        //TODO: the log_likelihood value needs to be saved somewhere for it to be tested.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data);
        //TODO: check for thr log_likelihood values
    }

}

TEST_CASE("ExaGeoStat API tests") {
TEST_GENERATE_DATA();

TEST_MODEL_DATA();

}