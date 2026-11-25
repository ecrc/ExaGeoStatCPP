// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestRPredictProblemSize.cpp
 * @brief Regression test for the problem-size configured by the R prediction adapter.
 * @details PredictionSetupHelper (used by predict_data / R_ExaGeoStatPredictData) configured the problem
 * size from the training locations only:
 *
 *     aConfigurations.SetProblemSize(data->GetLocations()->GetSize());   // == number of TRAINING points
 *
 * The prediction pipeline, however, works over the combined set of observed (training) and missing (test)
 * locations, and sizes its internal Z descriptor from the configured problem size. When the number of test
 * (missing) points is comparable to or larger than the number of training points, the problem size is too
 * small and the prediction fails / produces incorrect values.
 *
 * The existing R-adapter test never triggers this because it uses 16 training points and only 2 test
 * points (test << train). This test deliberately uses MORE test points than training points and checks the
 * fundamental kriging interpolation property: predicting at an observed location must return the observed
 * value (for an exact, nugget-free model). It uses the non-nugget univariate_matern_stationary kernel so it
 * is independent of any other kernel-specific behaviour.
 * @version 1.1.0
**/

#include <vector>

#include <catch2/catch_all.hpp>

#include <Rcpp-adapters/FunctionsAdapter.hpp>

using namespace std;

using namespace exageostat::adapters;
using namespace exageostat::common;

void TEST_R_PREDICT_PROBLEM_SIZE() {

    const string kernel = "univariate_matern_stationary";
    const string distance_matrix = "eg";
    const string dimension = "2D";
    const int dts = 3;
    const int lts = 0;
    const vector<double> estimated_theta = {1.0, 0.1, 0.5};

    // 9 training points with distinct measurements.
    vector<double> train_x = {0.10, 0.22, 0.31, 0.43, 0.55, 0.62, 0.74, 0.81, 0.93};
    vector<double> train_y = {0.12, 0.27, 0.39, 0.41, 0.58, 0.63, 0.71, 0.86, 0.95};
    vector<double> train_z = {-1.0, 2.0, 0.5, -0.5, 1.5, -1.5, 0.8, -0.8, 1.1};
    vector<vector<double>> train_data = {train_x, train_y, train_z};

    auto hardware = ExaGeoStatHardware("exact", 1, 0);

    SECTION("Prediction with more test points than training points") {

        // 12 test points = the 9 training locations + 3 extra. N_test (12) > N_train (9), so a problem size
        // computed from the training set alone is too small for the prediction descriptors.
        vector<double> test_x = train_x;
        vector<double> test_y = train_y;
        test_x.push_back(0.18); test_y.push_back(0.20);
        test_x.push_back(0.50); test_y.push_back(0.50);
        test_x.push_back(0.88); test_y.push_back(0.90);
        vector<vector<double>> test_data = {test_x, test_y};

        REQUIRE(test_x.size() > train_x.size());

        auto result = R_ExaGeoStatPredictData(kernel, distance_matrix, estimated_theta, dts, lts, dimension,
                                              train_data, test_data);

        // One prediction per test point.
        REQUIRE(result.size() == test_x.size());

        // Interpolation property: the prediction at each training location (the first 9 test points, which
        // coincide with the training locations) must recover the observed value for an exact, nugget-free
        // model. With the under-sized problem size the prediction is wrong; with the correct size it holds.
        for (size_t i = 0; i < train_x.size(); i++) {
            INFO("training location " << i << " at (" << train_x[i] << ", " << train_y[i] << ")");
            REQUIRE(result[i] == Catch::Approx(train_z[i]).margin(1e-4));
        }
    }
}

TEST_CASE("R predict_data problem-size regression") {
    TEST_R_PREDICT_PROBLEM_SIZE();
}
