
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

    // Init using unique ptr as when unique_ptr is destroyed/get out of scope, the resource is automatically claimed.
    auto syntheticDataConfigurations = new SyntheticDataConfigurations();
    SyntheticGenerator syntheticGenerator = SyntheticGenerator(syntheticDataConfigurations);

    SECTION("Spread Bytes")
    {
        uint16_t randomByte = INT16_MAX;
        REQUIRE(randomByte == 0x7FFF);
        uint64_t returnedByte = syntheticGenerator.SpreadBits(randomByte);
        // This because 7FFF will first be 16 hex = 64 bits
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 7FFF
        // 7FFF to bits is 0111111111111111
        // So, at the end it will be
        // ---- ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1
        REQUIRE(returnedByte == 0x0111111111111111);
    }
    SECTION("Reverse Spread Bytes")
    {
        uint64_t randomByte = 0x0111111111111111;
        uint16_t returnedByte = syntheticGenerator.ReverseSpreadBits(randomByte);
        REQUIRE(returnedByte == 0x7FFF);
    }
    SECTION("Spread & reverse 3D")
    {
        // Test spreading and shifting 3D and getting back values correctly.
        uint16_t x = INT16_MAX;
        uint16_t y = INT16_MAX;
        uint16_t z = INT16_MAX;
        uint64_t vectorZ;

        vectorZ = (syntheticGenerator.SpreadBits(z) << 2);
        // vector Z will be
        // ---- ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1 ---1
        // After shifting by 2
        // ---- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1-- -1--
        // Since 0100 is equal to 4 in hex then expected to be 0x0444444444444444
        REQUIRE(vectorZ == 0x0444444444444444);

        // Do the same for Y
        vectorZ += (syntheticGenerator.SpreadBits(y) << 1);
        // If vector Z was empty it will be
        // ---- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1-
        // But sine vectorZ is already contains Z. then by adding both we get
        // ---- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11- -11-
        // Since 0110 is equal to 6 in hex then expected to be 0x0666666666666666
        REQUIRE(vectorZ == 0x0666666666666666);

        // Lastly, Adding X
        vectorZ += syntheticGenerator.SpreadBits(x);
        // Adding X without shifting will result in
        // ---- -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111 -111
        // Since 0111 is equal to 7 in hex then expected to be 0x0777777777777777
        REQUIRE(vectorZ == 0x0777777777777777);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 2);

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
        vectorZ = (syntheticGenerator.SpreadBits(z_random) << 2) + (syntheticGenerator.SpreadBits(y_random) << 1) + syntheticGenerator.SpreadBits(x_random);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x_random =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y_random =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z_random =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 2);

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

        vectorZ = (syntheticGenerator.SpreadBits(z) << 2);
        // vector Z will be zeros
        REQUIRE(vectorZ == 0x0000000000000000);

        // Do the same for Y
        vectorZ += (syntheticGenerator.SpreadBits(y) << 1);
        // vector Z after shift by one will be
        // ---- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1- --1-
        // Since 0010 is equal to 2 in hex then expected to be 0x022222222222222
        REQUIRE(vectorZ == 0x0222222222222222);

        // Lastly, Adding X
        vectorZ += syntheticGenerator.SpreadBits(x);
        // Adding X without shifting will result in
        // ---- --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11 --11
        // Since 0011 is equal to 3 in hex then expected to be 0x0333333333333333
        REQUIRE(vectorZ == 0x0333333333333333);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 2);

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
        vectorZ = (syntheticGenerator.SpreadBits(z_random) << 2) + (syntheticGenerator.SpreadBits(y_random) << 1) + syntheticGenerator.SpreadBits(x_random);

        // Spreading is Done, Now reversing.
        uint16_t reversed_x_random =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 0);
        uint16_t reversed_y_random =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 1);
        uint16_t reversed_z_random =  syntheticGenerator.ReverseSpreadBits(vectorZ >> 2);

        REQUIRE(x_random == reversed_x_random);
        REQUIRE(y_random == reversed_y_random);
        REQUIRE(z_random == reversed_z_random);
    }
}

void TEST_GENERATE_LOCATIONS(){

    // Init using unique ptr as when unique_ptr is destroyed/get out of scope, the resource is automatically claimed.
    auto syntheticDataConfigurations = new SyntheticDataConfigurations();
    syntheticDataConfigurations->SetProblemSize(16);
    syntheticDataConfigurations->SetKernel("UnivariateMaternStationary");
    SyntheticGenerator syntheticGenerator = SyntheticGenerator(syntheticDataConfigurations);

    Locations locations;

    SECTION("2D Generation"){
        syntheticDataConfigurations->SetDimension(Dimension2D);
        syntheticGenerator.GenerateLocations();

        double* x = syntheticGenerator.GetLocations()->GetLocationX();
        double* y = syntheticGenerator.GetLocations()->GetLocationY();
        REQUIRE(syntheticGenerator.GetLocations()->GetLocationZ() == nullptr);

        for (auto i = 0; i < syntheticDataConfigurations->GetProblemSize(); i ++){
            REQUIRE( x[i] != 0 );
            REQUIRE( y[i] != 0 );
        }
    }

    SECTION("3D Generation"){
        syntheticDataConfigurations->SetDimension(Dimension3D);
        syntheticGenerator.GenerateLocations();

        double* x = syntheticGenerator.GetLocations()->GetLocationX();
        double* y = syntheticGenerator.GetLocations()->GetLocationY();
        double* z = syntheticGenerator.GetLocations()->GetLocationZ();

        for (auto i = 0; i < syntheticDataConfigurations->GetProblemSize(); i ++){
            REQUIRE( x[i] != 0 );
            REQUIRE( y[i] != 0 );
            REQUIRE( z[i] != 0 );
        }
    }
    SECTION("ST Generation"){
        syntheticDataConfigurations->SetDimension(DimensionST);
        syntheticDataConfigurations->SetTimeSlot(3);
        syntheticGenerator.GenerateLocations();

        double* x = syntheticGenerator.GetLocations()->GetLocationX();
        double* y = syntheticGenerator.GetLocations()->GetLocationY();
        double* z = syntheticGenerator.GetLocations()->GetLocationZ();

        for (auto i = 0; i < syntheticDataConfigurations->GetProblemSize() * syntheticDataConfigurations->GetTimeSlot(); i ++){
            REQUIRE( x[i] != 0 );
            REQUIRE( y[i] != 0 );
            REQUIRE( z[i] != 0 );
        }
    }
}

void TEST_HELPERS_FUNCTIONS(){
    auto syntheticDataConfigurations = new SyntheticDataConfigurations();
    SyntheticGenerator syntheticGenerator = SyntheticGenerator(syntheticDataConfigurations);

    SECTION("Uniform distribution"){
        double lowerRange = -0.4;
        double higherRange = 0.4;
        double uniformed_num = syntheticGenerator.UniformDistribution(lowerRange, higherRange);
        REQUIRE(uniformed_num > lowerRange);
        REQUIRE(uniformed_num < 1);
    }

    SECTION("Compare Uint32"){

        uint32_t num1 = 16;
        REQUIRE(syntheticGenerator.CompareUint64( num1, num1) == false);
        REQUIRE(syntheticGenerator.CompareUint64( num1, num1 + num1) == true);
        REQUIRE(syntheticGenerator.CompareUint64( num1 + num1, num1) == false);

    }
}

void TEST_GENERATION(){

        auto syntheticDataConfigurations = new SyntheticDataConfigurations();
        syntheticDataConfigurations->SetDimension(Dimension2D);
        syntheticDataConfigurations->SetProblemSize(2);
        syntheticDataConfigurations->SetKernel("UnivariateMaternStationary");
        SyntheticGenerator syntheticGenerator = SyntheticGenerator(syntheticDataConfigurations);
        syntheticGenerator.GenerateLocations();

        if (fabs(syntheticGenerator.GetLocations()->GetLocationX()[0] - 0.386069) >= 1e-6) {
            REQUIRE(false);
        }
        if (fabs(syntheticGenerator.GetLocations()->GetLocationY()[0] - 0.207752) >= 1e-6) {
            REQUIRE(false);
        }
        if (fabs(syntheticGenerator.GetLocations()->GetLocationX()[1] - 0.363241) >= 1e-6) {
            REQUIRE(false);
        }
        if (fabs(syntheticGenerator.GetLocations()->GetLocationY()[1] - 0.869383) >= 1e-6) {
            REQUIRE(false);
        }
        free(syntheticDataConfigurations);
}


TEST_CASE("Synthetic Data Generation values tests") {
    TEST_GENERATION();
}

TEST_CASE("Synthetic Data Generation tests") {
    TEST_SPREAD_REVERSED_BITS();
    TEST_GENERATE_LOCATIONS();
    TEST_HELPERS_FUNCTIONS();
}
