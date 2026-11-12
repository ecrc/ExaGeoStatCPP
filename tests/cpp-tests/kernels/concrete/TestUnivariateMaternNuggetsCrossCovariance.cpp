
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternNuggetsCrossCovariance.cpp
 * @brief Regression test for the location-index order of the UnivariateMaternNuggetsStationary kernel.
 * @details The existing kernel test exercises only the symmetric (train-train) covariance matrix, where
 * the matrix entry (i, j) is computed from the distance between two points drawn from the SAME location
 * set. Because the Euclidean distance is symmetric, an accidental swap of the two location indices is
 * invisible there: d(loc[i], loc[j]) == d(loc[j], loc[i]).
 *
 * This test instead builds a rectangular CROSS-covariance matrix between two DIFFERENT location sets,
 * exactly as the prediction path does when assembling C12 between observed and missing points. For two
 * distinct location sets, entry (i, j) must use d(location1[i], location2[j]); swapping the indices to
 * d(location1[j], location2[i]) produces a different (effectively transposed) matrix and therefore wrong
 * predictions.
 *
 * The expected values are computed independently of the kernel, using the closed-form Matern covariance
 * for nu = 0.5, which reduces to the exponential covariance sigma^2 * exp(-d / beta). This makes the
 * reference values self-evident and independent of any internal index convention.
 * @version 1.1.0
**/

#include <cmath>
#include <vector>

#include <catch2/catch_all.hpp>

#include <kernels/Kernel.hpp>
#include <data-units/Locations.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::plugins;

void TEST_NUGGETS_CROSS_COVARIANCE() {

    SECTION("UnivariateMaternNuggetsStationary cross-covariance respects location-index order") {

        // Matern parameters: sigma^2 = 1, beta (range) = 1, nu = 0.5, nugget = 0.5.
        // With nu = 0.5 the Matern covariance reduces to the exponential covariance:
        //     C(d) = sigma^2 * exp(-d / beta)        for d > 0
        //     C(0) = sigma^2 + nugget                for d = 0
        double theta[4] = {1.0, 1.0, 0.5, 0.5};
        const double sigma_square = theta[0];
        const double beta = theta[1];
        const double nugget = theta[3];

        const int rows = 3;    // number of "row" points, indexed from location set 1
        const int cols = 3;    // number of "column" points, indexed from location set 2

        // Two DISTINCT location sets so the cross-covariance is asymmetric.
        // Row points lie on the y-axis, column points lie on the x-axis.
        vector<double> row_x = {0.0, 0.0, 0.0};
        vector<double> row_y = {0.0, 1.0, 2.0};
        vector<double> col_x = {0.0, 3.0, 6.0};
        vector<double> col_y = {0.0, 0.0, 0.0};

        Locations<double> location1(rows, Dimension2D);
        location1.SetLocationX(*row_x.data(), rows);
        location1.SetLocationY(*row_y.data(), rows);

        Locations<double> location2(cols, Dimension2D);
        location2.SetLocationX(*col_x.data(), cols);
        location2.SetLocationY(*col_y.data(), cols);

        auto *kernel = PluginRegistry<Kernel<double>>::Create("UnivariateMaternNuggetsStationary", 1);
        REQUIRE(kernel != nullptr);

        // Output matrix in column-major layout: A[i + j * rows].
        vector<double> matrix(rows * cols, 0.0);

        // Distance metric 0 -> Euclidean (metric 1 is reserved for great-circle/Haversine).
        const int distance_metric = 0;
        kernel->GenerateCovarianceMatrix(matrix.data(), rows, cols, /*row offset*/ 0, /*column offset*/ 0,
                                         location1, location2, location1, theta, distance_metric);

        // Build the reference cross-covariance independently of the kernel.
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                const double dx = row_x[i] - col_x[j];
                const double dy = row_y[i] - col_y[j];
                const double distance = sqrt(dx * dx + dy * dy);
                const double expected = (distance == 0.0) ? (sigma_square + nugget)
                                                          : sigma_square * exp(-distance / beta);
                const double actual = matrix[i + j * rows];
                INFO("entry (" << i << ", " << j << "): distance = " << distance);
                REQUIRE(actual == Catch::Approx(expected).margin(1e-9));
            }
        }

        // The single most diagnostic entry (0, 1): row point (0,0), column point (3,0).
        //   Correct : d(location1[0]=(0,0), location2[1]=(3,0)) = 3 -> exp(-3) ~= 0.0498
        //   Swapped : d(location1[1]=(0,1), location2[0]=(0,0)) = 1 -> exp(-1) ~= 0.3679  (the bug)
        REQUIRE(matrix[0 + 1 * rows] == Catch::Approx(exp(-3.0)).margin(1e-9));

        delete kernel;
    }
}

TEST_CASE("UnivariateMaternNuggetsStationary cross-covariance index-order regression") {
    TEST_NUGGETS_CROSS_COVARIANCE();
}
