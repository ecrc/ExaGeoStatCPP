
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestSyntheticGenerator.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-08
**/

#include <libraries/catch/catch.hpp>
#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <iostream>

using namespace std;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::generators::Synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::common;

void TEST_SPREAD_REVERSED_BITS() {

    SyntheticDataConfigurations synthetic_data_configurations;

    synthetic_data_configurations.SetProblemSize(16);
    synthetic_data_configurations.SetKernel("UnivariateMaternStationary");
#ifdef EXAGEOSTAT_USE_CHAMELEON
    synthetic_data_configurations.SetComputation(exageostat::common::EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    synthetic_data_configurations.SetComputation(exageostat::common::TILE_LOW_RANK);
#endif

    SyntheticGenerator syntheticGenerator = SyntheticGenerator<double>(&synthetic_data_configurations);

    SECTION("Spread Bytes")
    {
        uint16_t randomByte = INT16_MAX;
        REQUIRE(randomByte == 0x7FFF);
        uint64_t returnedByte = SyntheticGenerator<double>::SpreadBits(randomByte);
        // This because 7FFF will first be 16 hex = 64 bits
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 7FFF
        // 7FFF to bits is 0111111111111111
        // So, at the end it will be
        // ---- ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1
        REQUIRE(returnedByte == 0x0111111111111111);
    }SECTION("Reverse Spread Bytes")
    {
        uint64_t randomByte = 0x0111111111111111;
        uint16_t returnedByte = SyntheticGenerator<double>::ReverseSpreadBits(randomByte);
        REQUIRE(returnedByte == 0x7FFF);
    }SECTION("Spread & reverse 3D")
    {
        // Test spreading and shifting 3D and getting back values correctly.
        uint16_t x = INT16_MAX;
        uint16_t y = INT16_MAX;
        uint16_t z = INT16_MAX;
        uint64_t vectorZ;

        vectorZ = (SyntheticGenerator<double>::SpreadBits(z) << 2);
        // vector Z will be
        // ---- ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1
        // After shifting by 2
        // ---- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1--
        // Since 0100 is equal to 4 in hex then expected to be 0x0444444444444444
        REQUIRE(vectorZ == 0x0444444444444444);

        // Do the same for Y
        vectorZ += (SyntheticGenerator<double>::SpreadBits(y) << 1);
        // If vector Z was empty it will be
        // ---- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1-
        // But sine vectorZ is already contains Z. then by adding both we get
        // ---- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11-
        // Since 0110 is equal to 6 in hex then expected to be 0x0666666666666666
        REQUIRE(vectorZ == 0x0666666666666666);

        // Lastly, Adding X
        vectorZ += SyntheticGenerator<double>::SpreadBits(x);
        // Adding X without shifting will result in
        // ---- -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111
        // Since 0111 is equal to 7 in hex then expected to be 0x0777777777777777
        REQUIRE(vectorZ == 0x0777777777777777);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 2);

        // What we reversed is what we send.
        REQUIRE(reversed_x == INT16_MAX);
        REQUIRE(reversed_y == INT16_MAX);
        REQUIRE(reversed_z == INT16_MAX);

        REQUIRE(x == reversed_x);
        REQUIRE(y == reversed_y);
        REQUIRE(z == reversed_z);

        // --------------------------------------------------------------------------------------------

        // Task with different values
        uint16_t x_random = 32007;
        uint16_t y_random = 37;
        uint16_t z_random = 22222;

        //Spreading
        vectorZ = (SyntheticGenerator<double>::SpreadBits(z_random) << 2) +
                  (SyntheticGenerator<double>::SpreadBits(y_random) << 1) +
                  SyntheticGenerator<double>::SpreadBits(x_random);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x_random = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y_random = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z_random = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 2);

        REQUIRE(x_random == reversed_x_random);
        REQUIRE(y_random == reversed_y_random);
        REQUIRE(z_random == reversed_z_random);
    }

    SECTION("Spread & reverse 3D")
    {
        // Test spreading and shifting 3D and getting back values correctly.
        uint16_t x = INT16_MAX;
        uint16_t y = INT16_MAX;
        uint16_t z = 0;
        uint64_t vectorZ;

        vectorZ = (SyntheticGenerator<double>::SpreadBits(z) << 2);
        // vector Z will be zeros
        REQUIRE(vectorZ == 0x0000000000000000);

        // Do the same for Y
        vectorZ += (SyntheticGenerator<double>::SpreadBits(y) << 1);
        // vector Z after shift by one will be
        // ---- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1-
        // Since 0010 is equal to 2 in hex then expected to be 0x022222222222222
        REQUIRE(vectorZ == 0x0222222222222222);

        // Lastly, Adding X
        vectorZ += SyntheticGenerator<double>::SpreadBits(x);
        // Adding X without shifting will result in
        // ---- --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11
        // Since 0011 is equal to 3 in hex then expected to be 0x0333333333333333
        REQUIRE(vectorZ == 0x0333333333333333);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 2);

        // What we reversed is what we send.
        REQUIRE(reversed_x == INT16_MAX);
        REQUIRE(reversed_y == INT16_MAX);
        REQUIRE(reversed_z == 0);

        REQUIRE(x == reversed_x);
        REQUIRE(y == reversed_y);
        REQUIRE(z == reversed_z);

        // --------------------------------------------------------------------------------------------

        // Task with different values
        uint16_t x_random = 32007;
        uint16_t y_random = 37;
        uint16_t z_random = 0;

        //Spreading
        vectorZ = (SyntheticGenerator<double>::SpreadBits(z_random) << 2) +
                  (SyntheticGenerator<double>::SpreadBits(y_random) << 1) +
                  SyntheticGenerator<double>::SpreadBits(x_random);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x_random = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y_random = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z_random = SyntheticGenerator<double>::ReverseSpreadBits(vectorZ >> 2);

        REQUIRE(x_random == reversed_x_random);
        REQUIRE(y_random == reversed_y_random);
        REQUIRE(z_random == reversed_z_random);
    }
}

void TEST_GENERATE_LOCATIONS() {

    SyntheticDataConfigurations synthetic_data_configurations;
    synthetic_data_configurations.SetProblemSize(16);
    synthetic_data_configurations.SetKernel("UnivariateMaternStationary");
#ifdef EXAGEOSTAT_USE_CHAMELEON
    synthetic_data_configurations.SetComputation(exageostat::common::EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    synthetic_data_configurations.SetComputation(exageostat::common::TILE_LOW_RANK);
#endif

    SyntheticGenerator syntheticGenerator = SyntheticGenerator<double>(&synthetic_data_configurations);

    Locations locations;

    SECTION("2D Generation") {
        synthetic_data_configurations.SetDimension(Dimension2D);
        syntheticGenerator.GenerateLocations();

        double *x = syntheticGenerator.GetLocations()->GetLocationX();
        double *y = syntheticGenerator.GetLocations()->GetLocationY();
        REQUIRE(syntheticGenerator.GetLocations()->GetLocationZ() == nullptr);

        for (auto i = 0; i < synthetic_data_configurations.GetProblemSize(); i++) {
            REQUIRE(x[i] != 0);
            REQUIRE(y[i] != 0);
        }
    }

    SECTION("3D Generation") {
        synthetic_data_configurations.SetDimension(Dimension3D);
        syntheticGenerator.GenerateLocations();

        double *x = syntheticGenerator.GetLocations()->GetLocationX();
        double *y = syntheticGenerator.GetLocations()->GetLocationY();
        double *z = syntheticGenerator.GetLocations()->GetLocationZ();

        for (auto i = 0; i < synthetic_data_configurations.GetProblemSize(); i++) {
            REQUIRE(x[i] != 0);
            REQUIRE(y[i] != 0);
            REQUIRE(z[i] != 0);
        }
    }SECTION("ST Generation") {
        synthetic_data_configurations.SetDimension(DimensionST);
        synthetic_data_configurations.SetTimeSlot(3);
        syntheticGenerator.GenerateLocations();

        double *x = syntheticGenerator.GetLocations()->GetLocationX();
        double *y = syntheticGenerator.GetLocations()->GetLocationY();
        double *z = syntheticGenerator.GetLocations()->GetLocationZ();

        for (auto i = 0;
             i < synthetic_data_configurations.GetProblemSize() * synthetic_data_configurations.GetTimeSlot(); i++) {
            REQUIRE(x[i] != 0);
            REQUIRE(y[i] != 0);
            REQUIRE(z[i] != 0);
        }
    }
}

void TEST_HELPERS_FUNCTIONS() {

    SyntheticDataConfigurations synthetic_data_configurations;
    synthetic_data_configurations.SetProblemSize(16);
    synthetic_data_configurations.SetKernel("UnivariateMaternStationary");
#ifdef EXAGEOSTAT_USE_CHAMELEON
    synthetic_data_configurations.SetComputation(exageostat::common::EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    synthetic_data_configurations.SetComputation(exageostat::common::TILE_LOW_RANK);
#endif

    SyntheticGenerator syntheticGenerator = SyntheticGenerator<double>(&synthetic_data_configurations);

    SECTION("Uniform distribution") {
        double lowerRange = -0.4;
        double higherRange = 0.4;
        double uniformed_num = SyntheticGenerator<double>::UniformDistribution(lowerRange, higherRange);
        REQUIRE(uniformed_num > lowerRange);
        REQUIRE(uniformed_num < 1);
    }

    SECTION("Compare Uint32") {

        uint32_t num1 = 16;
        REQUIRE(syntheticGenerator.CompareUint64(num1, num1) == false);
        REQUIRE(syntheticGenerator.CompareUint64(num1, num1 + num1) == true);
        REQUIRE(syntheticGenerator.CompareUint64(num1 + num1, num1) == false);
    }
}

void TEST_GENERATION() {

    SECTION("test Generated location") {

        SyntheticDataConfigurations synthetic_data_configurations;
        synthetic_data_configurations.SetDimension(Dimension2D);
        int N = 9;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernel("UnivariateMaternStationary");
        SyntheticGenerator syntheticGenerator = SyntheticGenerator<double>(&synthetic_data_configurations);

        // Initialize the seed manually with zero, to get the first generated seeded numbers.
        srand(0);
        syntheticGenerator.GenerateLocations();

        // The expected output of the locations.
        vector<double> x = {0.257389, 0.456062, 0.797269, 0.242161, 0.440742, 0.276432, 0.493965, 0.953933, 0.86952};
        vector<double> y = {0.138506, 0.238193, 0.170245, 0.579583, 0.514397, 0.752682, 0.867704, 0.610986, 0.891279};

        for (int i = 0; i < N; i++) {
            REQUIRE((syntheticGenerator.GetLocations()->GetLocationX()[i] - x[i]) == Approx(0.0).margin(1e-6));
            REQUIRE((syntheticGenerator.GetLocations()->GetLocationY()[i] - y[i]) == Approx(0.0).margin(1e-6));
        }

        // Now test re-generating locations again, but without modifying seed manually which will results in completely new locations values
        syntheticGenerator.GenerateLocations();
        for (int i = 0; i < N; i++) {
            REQUIRE((syntheticGenerator.GetLocations()->GetLocationX()[i] - x[i]) != Approx(0.0).margin(1e-6));
            REQUIRE((syntheticGenerator.GetLocations()->GetLocationY()[i] - y[i]) != Approx(0.0).margin(1e-6));
        }

        // Now if we modified seed again, we will get the first generated locations again.
        srand(0);
        syntheticGenerator.GenerateLocations();
        for (int i = 0; i < N; i++) {
            REQUIRE((syntheticGenerator.GetLocations()->GetLocationX()[i] - x[i]) == Approx(0.0).margin(1e-6));
            REQUIRE((syntheticGenerator.GetLocations()->GetLocationY()[i] - y[i]) == Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("Synthetic Data Generation tests") {
    TEST_SPREAD_REVERSED_BITS();
    TEST_GENERATE_LOCATIONS();
    TEST_HELPERS_FUNCTIONS();
    TEST_GENERATION();
}
