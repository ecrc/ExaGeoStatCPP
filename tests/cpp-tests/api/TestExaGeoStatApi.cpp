
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatApi.cpp
 * @brief Test suite for the ExaGeoStat APIs data generation functionality.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-08-07
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::configurations;
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
        synthetic_data_configurations.SetDenseTileSize(9);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        Configurations::SetVerbosity(QUIET_MODE);
        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{1, 0.1, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 4, 0); // Or you could use configurations.GetComputation().
        std::unique_ptr<exageostat::dataunits::ExaGeoStatData<double>> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(hardware, synthetic_data_configurations,
                                                                data);

        // Define the expected output for desk Z
        double expected_output_data[] = {-1.272336, -2.590700, 0.512143, -0.163880, 0.313504, -1.474411, 0.161705,
                                         0.623389, -1.341858, -1.054282, -1.669383, 0.219171, 0.971214, 0.538973,
                                         -0.752828, 0.290822};
        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                         DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        double diff;

        for (int i = 0; i < N; i++) {
            diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

void TEST_MODEL_DATA(Computation aComputation) {
    int seed = 0;
    srand(seed);

    // Create a new configurations object with the provided command line arguments
    Configurations configurations;
    int N = 16;
    configurations.SetProblemSize(N);
    configurations.SetKernelName("UnivariateMaternStationary");
    int dts = 8;
    configurations.SetDenseTileSize(dts);
    configurations.SetComputation(aComputation);
    configurations.SetMaxMleIterations(3);
    configurations.SetTolerance(pow(10, -4));

    vector<double> lb{0.1, 0.1, 0.1};
    configurations.SetLowerBounds(lb);
    configurations.SetStartingTheta(lb);

    vector<double> ub{5, 5, 5};
    configurations.SetUpperBounds(ub);

    vector<double> initial_theta{1, 0.1, 0.5};
    configurations.SetInitialTheta(initial_theta);

    double expected = 0;
    if (aComputation == EXACT_DENSE) {
        expected = -24.026000;
    } else if (aComputation == DIAGONAL_APPROX) {
        expected = -24.028197;
        configurations.SetBand(1);
    } else if (aComputation == TILE_LOW_RANK) {
        expected = -24.0049327;
        configurations.SetLowTileSize(dts);
        configurations.SetMaxRank(500);
    }

    SECTION("Data Modeling")
    {
        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(aComputation, 4, 0); // Or you could use configurations.GetComputation().
        std::unique_ptr<exageostat::dataunits::ExaGeoStatData<double>> data = std::make_unique<exageostat::dataunits::ExaGeoStatData<double>>(configurations.GetProblemSize(), configurations.GetDimension());


        //initiating the matrix of the CHAMELEON Descriptor Z.
        auto *z_matrix = new double[N]{-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
                                       -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
                                       0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
                                       -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                                       0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
                                       0.290822066007430102};
        //creating locations x and y.
        auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                         0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                         0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                         0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                         0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                         0.877592126344701295};

        auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                         0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                         0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                         0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                         0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                         0.942824444953078489};

        data->GetLocations()->SetLocationX(*location_x, N);
        data->GetLocations()->SetLocationY(*location_y, N);

        double log_likelihood = exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations,
                                                                                            data, z_matrix);
        REQUIRE((log_likelihood - expected) == Catch::Approx(0.0).margin(1e-6));

        delete[] location_x;
        delete[] location_y;
        delete[] z_matrix;
    }SECTION("Data Generation and Modeling")
    {
        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(aComputation, 4, 0); // Or you could use configurations.GetComputation().
        std::unique_ptr<exageostat::dataunits::ExaGeoStatData<double>> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(hardware, configurations,
                                                                data);
        double log_likelihood = exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations,
                                                                                            data);
        REQUIRE((log_likelihood - expected) == Catch::Approx(0.0).margin(1e-6));
    }
}

void TEST_PREDICTION() {
    Configurations configurations;
    configurations.SetUnknownObservationsNb(4);
    int N = 16;
    configurations.SetProblemSize(N);
    configurations.SetKernelName("UnivariateMaternStationary");
    int dts = 8;

    configurations.SetDenseTileSize(dts);
    configurations.SetComputation(EXACT_DENSE);
    configurations.SetMaxMleIterations(3);
    configurations.SetTolerance(pow(10, -4));

    vector<double> lb{0.1, 0.1, 0.1};
    configurations.SetLowerBounds(lb);
    configurations.SetStartingTheta(lb);
    vector<double> ub{5, 5, 5};
    configurations.SetUpperBounds(ub);
    vector<double> initial_theta{1, 0.1, 0.5};
    configurations.SetInitialTheta(initial_theta);
    vector<double> estimated_theta{-1, -1, -1};
    configurations.SetEstimatedTheta(estimated_theta);
    Configurations::SetVerbosity(QUIET_MODE);

    auto hardware = ExaGeoStatHardware(EXACT_DENSE, configurations.GetCoresNumber(), configurations.GetGPUsNumbers());
    std::unique_ptr<exageostat::dataunits::ExaGeoStatData<double>> data = std::make_unique<exageostat::dataunits::ExaGeoStatData<double>>(configurations.GetProblemSize(), configurations.GetDimension());

    auto *z_matrix = new double[N]{-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
                                   -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
                                   0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
                                   -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                                   0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
                                   0.290822066007430102};
    //creating locations x and y.
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                     0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                     0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                     0.877592126344701295};
    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                     0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                     0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                     0.942824444953078489};

    data->GetLocations()->SetLocationX(*location_x, N);
    data->GetLocations()->SetLocationY(*location_y, N);

    SECTION("Test Prediction - MSPE ONLY")
    {
        vector<double> new_estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(new_estimated_theta);
        configurations.SetIsMSPE(true);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data, z_matrix);

    }SECTION("Test Prediction - IDW ONLY") {
        configurations.SetIsMSPE(false);
        configurations.SetIsIDW(true);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data, z_matrix);
    }SECTION("Test Prediction - MLOE_MMOM ONLY") {
        configurations.SetIsMSPE(false);
        configurations.SetIsIDW(false);
        vector<double> new_estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(new_estimated_theta);
        configurations.SetIsMLOEMMOM(true);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data, z_matrix);
    }SECTION("Test Prediction - FISHER\n") {
        configurations.SetIsMLOEMMOM(false);
        configurations.SetIsFisher(true);
        vector<double> new_estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(new_estimated_theta);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data, z_matrix);
    }SECTION("Test Prediction - ALL OPERATIONS") {
        vector<double> new_estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(new_estimated_theta);
        configurations.SetIsMSPE(true);
        configurations.SetIsIDW(true);
        configurations.SetIsMLOEMMOM(true);
        configurations.SetIsFisher(true);
        // Setting Estimated with initial theta will require mloe_mmom to be zero
        configurations.SetEstimatedTheta(initial_theta);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data, z_matrix);
    }SECTION("Test Prediction - ALL MODULES") {
        configurations.SetIsMSPE(true);
        configurations.SetIsIDW(true);
        configurations.SetIsFisher(true);
        configurations.SetIsMLOEMMOM(true);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(hardware, configurations,
                                                                data);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data, z_matrix);
    }
    delete[] location_x;
    delete[] location_y;
    delete[] z_matrix;
}

TEST_CASE("ExaGeoStat API tests") {
    TEST_GENERATE_DATA();
    TEST_MODEL_DATA(EXACT_DENSE);
    TEST_MODEL_DATA(DIAGONAL_APPROX);
    TEST_MODEL_DATA(TILE_LOW_RANK);
    TEST_PREDICTION();
}