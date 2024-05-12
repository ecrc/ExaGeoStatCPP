// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestAllRFunctions.cpp
 * @brief Test suite for R/Rcpp adapters in C++.
 * @details This file contains tests for R methods adapted for use in C++ within the ExaGeoStat software package.
 * It includes tests for initializing arguments, loading data, modeling data, and predicting data using Rcpp adapters.
 * The test suite specifically verifies the integration and functionality of R methods through the Rcpp interface
 * within a C++ environment, ensuring that statistical models and algorithms are correctly initialized,
 * data is accurately loaded and processed, and predictions are properly executed.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-09
**/

#include <catch2/catch_all.hpp>

#include <Rcpp-adapters/FunctionsAdapter.hpp>

using namespace std;

using namespace exageostat::adapters;
using namespace exageostat::common;

void TEST_DATA_GENERATION() {

    const string computation = "exact";
    const string dimension = "2D";
    const string kernel = "univariate_matern_stationary";
    const string distance_matrix = "eg";
    const int cores_number = 1;
    const int gpus_number = 0;
    const int problem_size = 9;
    const int dts = 3;
    const int lts = 0;
    const vector<double> initial_theta = {1, 0.1, 0.5};
    const vector<double> lower_bound = {0.1, 0.1, 0.1};
    const vector<double> upper_bound = {5, 5, 5};
    const vector<double> estimated_theta = {0.9, 0.2, 0.5};
    srand(0);

    auto hardware = ExaGeoStatHardware(computation, cores_number, gpus_number);
    auto exageostat_data = R_ExaGeoStatLoadData(kernel, initial_theta, distance_matrix, problem_size, 0, dts, lts,
                                                dimension, "", "", "", "");

    vector<double> x = {0.257389, 0.456062, 0.797269, 0.242161, 0.440742, 0.276432, 0.493965, 0.953933, 0.86952};
    vector<double> y = {0.138506, 0.238193, 0.170245, 0.579583, 0.514397, 0.752682, 0.867704, 0.610986, 0.891279};
    vector<double> z = {-1.27234, -2.47547, 0.54585, -0.120985, 0.242569, -1.54421, 0.0986468, 0.779835, -1.48139};

    for (int i = 0; i < problem_size; i++) {
        REQUIRE((exageostat_data->GetLocations()->GetLocationX()[i] - x[i]) == Catch::Approx(0.0).margin(1e-6));
        REQUIRE((exageostat_data->GetLocations()->GetLocationY()[i] - y[i]) == Catch::Approx(0.0).margin(1e-6));
        REQUIRE((exageostat_data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z)[i] -
                 z[i]) == Catch::Approx(0.0).margin(1e-5));
    }
    delete exageostat_data;
}

void TEST_ALL_R_METHODS() {

    const string computation = "exact";
    const int cores_number = 1;
    const int gpus_number = 0;
    const string dimension = "2D";
    const string kernel = "univariate_matern_stationary";
    const string distance_matrix = "eg";
    const int dts = 8;
    const int lts = 0;

    const vector<double> initial_theta = {1, 0.1, 0.5};
    const vector<double> estimated_theta = {0.9, 0.2, 0.5};
    const vector<double> lower_bound = {0.1, 0.1, 0.1};
    const vector<double> upper_bound = {5, 5, 5};

    vector<vector<double>> train_data = {
            {0.193041886015106440,  0.330556191348134576,  0.181612878614480805, 0.370473792629892440,
                    0.652140077821011688, 0.806332494087129037,  0.553322652018005678, 0.800961318379491916,
                    0.207324330510414295,  0.347951476310368490,  0.092042420080872822,  0.465445944914930965,
                    0.528267338063630132, 0.974792095826657490, 0.552452887769893985,  0.877592126344701295},
            {0.103883421072709245,  0.135790035858701447,  0.434683756771190977, 0.400778210116731537,
                    0.168459601739528508, 0.105195696955825133,  0.396398870832379624, 0.296757457846952011,
                    0.564507515068284116,  0.627679865720607300,  0.928648813611047563,  0.958236057068741931,
                    0.573571374074921758, 0.568657969024185528, 0.935835812924391552,  0.942824444953078489},
            {-1.272336140360187606, -2.590699695867695773, 0.512142584178685967, -0.163880452049749520,
                    0.313503633252489700, -1.474410682226017677, 0.161705025505231914, 0.623389205185149065,
                    -1.341858445399783495, -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                    0.971213790000161170, 0.538973474182433021, -0.752828466476077041, 0.290822066007430102}
    };
    vector<vector<double>> test_data = {
            {0.193041886015106440, 0.330556191348134576},
            {0.103883421072709245, 0.135790035858701447}
    };

    auto hardware = ExaGeoStatHardware(computation, cores_number, gpus_number);

    SECTION("Prediction Method") {

        // Since we're testing with the same x and y that we sent to the train data, then we expect the same output.
        auto result = R_ExaGeoStatPredictData(kernel, distance_matrix, estimated_theta, dts, lts, dimension, train_data,
                                              test_data);
        for (int i = 0; i < result.size(); i++) {
            REQUIRE((result[i] - train_data[2][i]) == Catch::Approx(0.0).margin(1e-6));
        }

        // Now let's test with new values
        test_data = {
                {0.2, 0.3},
                {0.2, 0.2}
        };
        result = R_ExaGeoStatPredictData(kernel, distance_matrix, estimated_theta, dts, lts, dimension, train_data,
                                         test_data);
        vector<double> expected_output = {-1.04437, -1.63133};
        for (int i = 0; i < result.size(); i++) {
            REQUIRE((result[i] - expected_output[i]) == Catch::Approx(0.0).margin(1e-5));
        }

    }SECTION("PREDICTION - Fisher") {

        auto result = R_ExaGeoStatFisher(kernel, distance_matrix, estimated_theta, dts, lts, dimension, train_data,
                                         test_data);

        vector<double> expected_output = {0.143552740378975224, 0.022642203822768981, 0.013749788660229965,
                                          0.022642203822769023, 0.143052423579629107, -0.371722296732991009,
                                          0.013749788660229986, -0.371722296732990953, 1.101996735797240667};

        for (int i = 0; i < result.size(); i++) {
            REQUIRE((result[i] - expected_output[i]) == Catch::Approx(0.0).margin(1e-6));
        }
    }

    SECTION("PREDICTION - MLOE-MMOM") {
        test_data = {
                {0.347951, 0.62768},
                {0.806332, 0.105196},
        };
        auto result = R_ExaGeoStatMLOE_MMOM(kernel, distance_matrix, estimated_theta, initial_theta, dts, lts,
                                            dimension, train_data, test_data);
        vector<double> expected_output = {0.0529825, -0.408526};
        for (int i = 0; i < result.size(); i++) {
            REQUIRE((result[i] - expected_output[i]) == Catch::Approx(0.0).margin(1e-6));
        }

        result = R_ExaGeoStatMLOE_MMOM("univariate_matern_stationary", "eg", initial_theta, initial_theta, 8, 0, "2D",
                                       train_data, test_data);
        for (double i: result) {
            REQUIRE((i - 0) == Catch::Approx(0.0).margin(1e-6));
        }
    }

    SECTION("PREDICTION - IDW") {

        vector<double> test_measurements = {
                -1.05428, -1.47441
        };
        auto result = R_ExaGeoStatIDW(kernel, distance_matrix, estimated_theta, dts, lts, dimension, train_data,
                                      test_data, test_measurements);
        vector<double> expected_output = {0.447903, 0.100866, 0.79494};
        for (int i = 0; i < result.size(); i++) {
            REQUIRE((result[i] - expected_output[i]) == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("Test R/Rcpp adapters in C++") {
    TEST_DATA_GENERATION();
    TEST_ALL_R_METHODS();
}