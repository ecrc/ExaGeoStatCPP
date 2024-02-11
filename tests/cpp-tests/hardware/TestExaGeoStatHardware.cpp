
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatHardware.cpp
 * @brief Unit tests for the ExaGeoStatHardware class in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the ExaGeoStatHardware class
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2024-01-24
**/

#include <catch2/catch_all.hpp>

#include <hardware/ExaGeoStatHardware.hpp>

using namespace exageostat::hardware;
using namespace exageostat::common;
using namespace exageostat::hardware;

void TEST_HARDWARE_CONSTRUCTION() {

    //Initialize multiple instances of the hardware, and verify chameleon context is not null.
    auto hardware_1 = ExaGeoStatHardware(EXACT_DENSE, 4, 0);
    REQUIRE(hardware_1.GetChameleonContext() != nullptr);

    auto hardware_2 = ExaGeoStatHardware(DIAGONAL_APPROX, 4, 0);
    REQUIRE(hardware_2.GetChameleonContext() != nullptr);
    REQUIRE(hardware_1.GetChameleonContext() == hardware_2.GetChameleonContext());

    //In case of HICMA initialize a TLR hardware and verify non-null context,
    // otherwise exception is raised in case of initialization.
#ifdef USE_HICMA
    auto hardware_3 = ExaGeoStatHardware(TILE_LOW_RANK, 4, 0);
    REQUIRE(hardware_3.GetHicmaContext() != nullptr);

    auto hardware_4 = ExaGeoStatHardware(TILE_LOW_RANK, 4, 0);
    REQUIRE(hardware_4.GetHicmaContext() != nullptr);
    REQUIRE(hardware_3.GetHicmaContext() == hardware_4.GetHicmaContext());
#else
    REQUIRE_THROWS(ExaGeoStatHardware(TILE_LOW_RANK, 4, 0));
#endif

}

void TEST_STATIC_CONTEXT() {

    //Initialize multiple instances of the hardware: EXACT_DENSE,DIAGONAL_APPROX ,and verify they all have same static Chameleon context.
    auto hardware_1 = ExaGeoStatHardware(EXACT_DENSE, 4, 0);
    auto context_hardware_1 = ExaGeoStatHardware::GetContext(EXACT_DENSE);

    auto hardware_2 = ExaGeoStatHardware(EXACT_DENSE, 1, 0);
    auto context_hardware_2 = ExaGeoStatHardware::GetContext(EXACT_DENSE);
    REQUIRE(context_hardware_1 == context_hardware_2);

    auto hardware_3 = ExaGeoStatHardware(DIAGONAL_APPROX, 7, 0);
    auto context_hardware_3 = ExaGeoStatHardware::GetContext(DIAGONAL_APPROX);
    REQUIRE(context_hardware_2 == context_hardware_3);

    auto hardware_4 = ExaGeoStatHardware(DIAGONAL_APPROX, 3, 0);
    auto context_hardware_4 = ExaGeoStatHardware::GetContext(DIAGONAL_APPROX);
    REQUIRE(context_hardware_3 == context_hardware_4);

    //In case of HICMA initialize multiple instances of the TLR hardware.
#ifdef USE_HICMA
    // Verify they have same static HICMA context
    auto hardware_5 = ExaGeoStatHardware(TILE_LOW_RANK, 4, 0);
    REQUIRE(hardware_5.GetContext(TILE_LOW_RANK) != nullptr);

    auto hardware_6 = ExaGeoStatHardware(TILE_LOW_RANK, 4, 0);
    REQUIRE(hardware_6.GetContext(TILE_LOW_RANK) != nullptr);
    REQUIRE(hardware_5.GetContext(TILE_LOW_RANK) == hardware_6.GetContext(TILE_LOW_RANK));

    // Verify they have same Chameleon context as the EXACT_DENSE and DIAGONAL_APPROX hardware
    REQUIRE(hardware_5.GetContext(EXACT_DENSE) == hardware_6.GetContext(EXACT_DENSE));
    REQUIRE(hardware_5.GetContext(EXACT_DENSE) == hardware_1.GetContext(EXACT_DENSE));
#endif

}

void TEST_CONTEXT_GETTER() {

    // Initialize multiple instances of hardware, and verify context getter results.
    auto hardware_1 = ExaGeoStatHardware(EXACT_DENSE, 1, 0);
    auto hardware_2 = ExaGeoStatHardware(DIAGONAL_APPROX, 1, 0);

    // Chameleon context is always non-null
    REQUIRE(hardware_1.GetContext(EXACT_DENSE) != nullptr);
    REQUIRE(hardware_2.GetContext(DIAGONAL_APPROX) != nullptr);

    //In case of HICMA, initialize TLR hardware.
#ifdef USE_HICMA
    auto hardware_3 = ExaGeoStatHardware(TILE_LOW_RANK, 4, 0);

    // Hicma context and Chameleon context are always non-null
    REQUIRE(hardware_3.GetContext(TILE_LOW_RANK) != nullptr);
    REQUIRE(hardware_3.GetContext(EXACT_DENSE) != nullptr);
    REQUIRE(hardware_3.GetContext(DIAGONAL_APPROX) != nullptr);
#else
    // Otherwise an exception is raised
    REQUIRE_THROWS(hardware_1.GetContext(TILE_LOW_RANK));
    REQUIRE_THROWS(hardware_2.GetContext(TILE_LOW_RANK));
#endif
}

TEST_CASE("ExaGeoStat Hardware Tests") {
    TEST_HARDWARE_CONSTRUCTION();
    TEST_CONTEXT_GETTER();
    TEST_STATIC_CONTEXT();
}

