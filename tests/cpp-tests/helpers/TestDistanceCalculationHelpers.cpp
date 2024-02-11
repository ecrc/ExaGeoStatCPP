
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestDistanceCalculationHelpers.cpp
 * @brief Unit tests for the DistanceCalculationHelpers.cpp in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the class DistanceCalculationHelpers.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2024-01-24
**/

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>

#include <common/Definitions.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::common;
using namespace exageostat::helpers;
using namespace exageostat::dataunits;

void TEST_RADIAN_CONVERSION() {
    exageostat::helpers::DistanceCalculationHelpers<double> distance_helper;
    REQUIRE(distance_helper.DegreeToRadian(0) == Catch::Approx(0.0));
    REQUIRE(distance_helper.DegreeToRadian(90) == Catch::Approx(PI / 2));
    REQUIRE(distance_helper.DegreeToRadian(180) == Catch::Approx(PI));
}

void TEST_CALCULATE_DISTANCE() {

    SECTION("2D - Euclidean & Haversine Distance") {
        int size = 2;

        // 2D - Locations initializations
        auto *location_x = new double[size]{10.0, 4.0};
        auto *location_y = new double[size]{12.0, 4.0};

        Locations<double> location1(size, Dimension2D);
        Locations<double> location2(size, Dimension2D);

        location1.SetLocationX(*location_x, size);
        location1.SetLocationY(*location_y, size);
        location2.SetLocationX(*location_x, size);
        location2.SetLocationY(*location_y, size);

        int idx1 = 0; // index of x-coordinate
        int idx2 = 1; // index of y-coordinate
        int flagZ = 0; // 2D case

        int distanceMetric = 0; // Default - Use Euclidean distance
        REQUIRE(DistanceCalculationHelpers<double>::CalculateDistance(location1, location2, idx1, idx2, distanceMetric,
                                                                      flagZ) == 10);

        distanceMetric = 1; // Use Haversine distance
        auto distance_earth = DistanceCalculationHelpers<double>::DistanceEarth(location_x[idx1], location_y[idx1],
                                                                                location_x[idx2], location_y[idx2]);
        REQUIRE(DistanceCalculationHelpers<double>::CalculateDistance(location1, location2, idx1, idx2, distanceMetric,
                                                                      flagZ) == distance_earth);

        delete[] location_x;
        delete[] location_y;
    }SECTION("3D - Euclidean & Haversine Distance") {
        int size = 3;

        // 3D - Locations initializations
        auto *location_x = new double[size]{10.0, 4.0};
        auto *location_y = new double[size]{12.0, 4.0};
        auto *location_z = new double[size]{1, 2};

        Locations<double> location1(size, Dimension3D);
        Locations<double> location2(size, Dimension3D);

        location1.SetLocationX(*location_x, size);
        location1.SetLocationY(*location_y, size);
        location1.SetLocationZ(*location_z, size);

        location2.SetLocationX(*location_x, size);
        location2.SetLocationY(*location_y, size);
        location2.SetLocationZ(*location_z, size);

        int idx1 = 0; // index of x-coordinate
        int idx2 = 1; // index of y-coordinate
        int flagZ = 1; // 2D case

        int distanceMetric = 0; // Use Euclidean distance
        REQUIRE(DistanceCalculationHelpers<double>::CalculateDistance(location1, location2, idx1, idx2, distanceMetric,
                                                                      flagZ) == Catch::Approx(10.05).epsilon(0.01));

        distanceMetric = 1; // Use Haversine distance
        REQUIRE_THROWS(
                DistanceCalculationHelpers<double>::CalculateDistance(location1, location2, idx1, idx2, distanceMetric,
                                                                      flagZ));

        delete[] location_x;
        delete[] location_y;
        delete[] location_z;
    }
}

void TEST_DISTANCE_EARTH() {
    SECTION("Distance between two identical points") {

        double lat = 0.1, lon = 0.2;
        REQUIRE(DistanceCalculationHelpers<double>::DistanceEarth(lat, lon, lat, lon) == Catch::Approx(0.0));

    }SECTION("Distance between North Pole and South Pole") {

        double north_pole_Lat = 90.0, pole_Lon = 0.0;
        double south_pole_Lat = -90.0;
        REQUIRE(DistanceCalculationHelpers<double>::DistanceEarth(north_pole_Lat, pole_Lon, south_pole_Lat, pole_Lon) ==
                Catch::Approx(EARTH_RADIUS * M_PI));

    }SECTION("Distance between two points on the equator 90 degrees apart") {

        double equator_Lat = 0.0;
        double lon1 = 0.0, lon2 = 90.0;
        REQUIRE(DistanceCalculationHelpers<double>::DistanceEarth(equator_Lat, lon1, equator_Lat, lon2) ==
                Catch::Approx(EARTH_RADIUS * M_PI / 2));

    }SECTION("Distance between two arbitrary points") {

        // Arbitrary location coordinates
        double location1_lat = 40.7128;
        double location1_lon = -74.0060;
        double location2_lat = 34.0522;
        double location2_lon = -118.2437;

        // The expected distance is approximately 3940 kilometers
        double expectedDistance = 3940.0;
        double calculatedDistance = DistanceCalculationHelpers<double>::DistanceEarth(location1_lat, location1_lon,
                                                                                      location2_lat, location2_lon);
        REQUIRE(calculatedDistance == Catch::Approx(expectedDistance).epsilon(0.01));

    }
}

TEST_CASE("Degree to Radian Conversion") {
    TEST_RADIAN_CONVERSION();
    TEST_CALCULATE_DISTANCE();
    TEST_DISTANCE_EARTH();
}