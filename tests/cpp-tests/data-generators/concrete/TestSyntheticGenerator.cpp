
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

void TEST_GENERATE_LOCATIONS() {

    // Init using unique ptr as when unique_ptr is destroyed/get out of scope, the resource is automatically claimed.
    auto syntheticDataConfigurations = new SyntheticDataConfigurations();
    syntheticDataConfigurations->SetProblemSize(16);

    SyntheticGenerator syntheticGenerator = SyntheticGenerator(syntheticDataConfigurations);

    SECTION("SpreadBits")
    {
        uint32_t randomByte = 10;
        cout << syntheticGenerator.SpreadBits(randomByte) << endl;

    }
}

TEST_CASE("Synthetic Data Generation") {
    TEST_GENERATE_LOCATIONS();
}
