
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestChameleonImplmentationDense.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-06
**/
#include <libraries/catch/catch.hpp>
#include <iostream>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include "control/context.h"

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;


void INIT_HARDWARE(){
    ChameleonImplementationDense<double> chameleonImpl;
    chameleonImpl.ExaGeoStatInitContext(4, 0);
    CHAM_context_t *chameleonContext = chameleon_context_self();
    REQUIRE(chameleonContext != nullptr);
}
// Test that the function initializes all the required descriptors without errors.
void TEST_INITIALIZETION(){
    auto* syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("Single"){
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(1000);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 2);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct().size() == 1);

        REQUIRE(syntheticDataConfigurations->GetDescriptorC()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }

    SECTION("Double"){
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(1024);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 4);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 3);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct().size() == 3);

        for(auto & descriptorC : syntheticDataConfigurations->GetDescriptorC()){
            REQUIRE(descriptorC != nullptr);
        }
        for (auto & i : syntheticDataConfigurations->GetDescriptorZ()){
            REQUIRE(i != nullptr);
        }
        for (auto & i : syntheticDataConfigurations->GetDescriptorProduct()){
            REQUIRE(i != nullptr);
        }
        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }
}
//Test that the function initializes the CHAM_descriptorC descriptor correctly.
void TEST_CHAMELEON_DESCRIPTORS_C() {
    auto* syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("SINGLE"){
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetDenseTileSize(512);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

    }
}

TEST_CASE("Chameleon Implementation Dense"){
    INIT_HARDWARE();
    TEST_INITIALIZETION();
    TEST_CHAMELEON_DESCRIPTORS_C();
}

/*

// Test case 2: Test that the function initializes the CHAM_descriptorC descriptor correctly.
TEST_CASE("InitiateDescriptorsTest2", "[ChameleonImplementationDense]") {
exageostat::linearAlgebra::dense::ChameleonImplementationDense<double> solver;
solver.mpConfigurations->SetProblemSize(100);
solver.mpConfigurations->SetDenseTileSize(10);
solver.InitiateDescriptors();

auto* CHAM_descriptorC = (CHAM_desc_t*) solver.mpConfigurations->GetDescriptorC()[0];
int N = solver.mpConfigurations->GetProblemSize() * solver.mpConfigurations->GetP();
int dts = solver.mpConfigurations->GetDenseTileSize();
int pGrid = solver.mpConfigurations->GetPGrid();
int qGrid = solver.mpConfigurations->GetQGrid();

REQUIRE(CHAM_descriptorC->m == N);
REQUIRE(CHAM_descriptorC->n == N);
REQUIRE(CHAM_descriptorC->mb == dts);
REQUIRE(CHAM_descriptorC->nb == dts);
REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
REQUIRE(CHAM_descriptorC->i == 0);
REQUIRE(CHAM_descriptorC->j == 0);
REQUIRE(CHAM_descriptorC->mt == N / dts);
REQUIRE(CHAM_descriptorC->nt == N / dts);
REQUIRE(CHAM_descriptorC->m0 == 0);
REQUIRE(CHAM_descriptorC->n0 == 0);
REQUIRE(CHAM_descriptorC->lm == N);
REQUIRE(CHAM_descriptorC->ln == N);
REQUIRE(CHAM_descriptorC->p == pGrid);
REQUIRE(CHAM_descriptorC->q == qGrid);

for (int i = 0; i < CHAM_descriptorC->mt * CHAM_descriptorC->nt * CHAM_descriptorC->bsiz; i++) {
REQUIRE(CHAM_descriptorC->mat[i] == 0);
}
}

// Test case 3: Test that the function initializes the CHAM_descriptorZ descriptor correctly.
TEST_CASE("InitiateDescriptorsTest3", "[ChameleonImplementationDense]") {
exageostat::linearAlgebra::dense::ChameleonImplementationDense<double> solver;
solver.mpConfigurations->SetProblemSize(100);
solver.mpConfigurations->SetDenseTileSize(10);
solver.InitiateDescriptors();

auto* CHAM_descriptorZ = (CHAM_desc_t*) solver.mpConfigurations->GetDescriptorZ()[0];
int N = solver.mpConfigurations->GetProblemSize() * solver.mpConfigurations->GetP();
int dts = solver.mpConfigurations->GetDenseTileSize();
int pGrid = solver.mpConfigurations->GetPGrid();
int qGrid = solver.mpConfigurations->GetQGrid();

REQUIRE(CHAM_descriptorZ->m == N);
REQUIRE(CHAM_descriptorZ->n == 1);
REQUIRE(CHAM_descriptorZ->mb == dts);
REQUIRE(CHAM_descriptorZ->nb == dts);
REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
REQUIRE(CHAM_descriptorZ->i == 0);
REQUIRE(CHAM_descriptorZ->j == 0);
REQUIRE(CHAM_descriptorZ->mt == N / dts);
REQUIRE(CHAM_descriptorZ->nt == 1);
REQUIRE(CHAM_descriptorZ->m0 == 0);
REQUIRE(CHAM_descriptorZ->n0 == 0);
REQUIRE(CHAM_descriptorZ->lm == N);
REQUIRE(CHAM_descriptorZ->ln == 1);
REQUIRE(CHAM_descriptorZ->p == pGrid);
REQUIRE(CHAM_descriptorZ->q == qGrid);




 // Test that the function initializes all the required descriptors without errors.
void TEST_INITIALIZETION(){

    auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);
    auto* syntheticDataConfigurations = new SyntheticDataConfigurations();
    linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

    SECTION("Single"){
        syntheticDataConfigurations->SetProblemSize(1024);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);

        // Check the values set in the CHAM_desc_s struct for the C descriptor
        auto* Cdesc = (CHAM_desc_t*) syntheticDataConfigurations->GetDescriptorC()[0];
        REQUIRE(Cdesc->m == 1024);
        REQUIRE(Cdesc->n == 1024);
        REQUIRE(Cdesc->mb == 64);
        REQUIRE(Cdesc->nb == 64);
        REQUIRE(Cdesc->bsiz == 4096);
        REQUIRE(Cdesc->lmt*64 == 1024);
        REQUIRE(Cdesc->lnt*64 == 1024);
        REQUIRE(Cdesc->p == 1);
        REQUIRE(Cdesc->q == 1);
        REQUIRE(Cdesc->i == 1);
        REQUIRE(Cdesc->j == 1);
        REQUIRE(Cdesc->m0 == 0);
        REQUIRE(Cdesc->n0 == 0);
        REQUIRE(Cdesc->mt == 16);
        REQUIRE(Cdesc->nt == 16);
        REQUIRE(Cdesc->m1 == 1024);
        REQUIRE(Cdesc->n1 == 1024);

        // Check the values set in the CHAM_desc_s struct for the Z descriptor
        auto* Zdesc = (CHAM_desc_t*) syntheticDataConfigurations->GetDescriptorZ()[0];
        REQUIRE(Zdesc->m == 1024);
        REQUIRE(Zdesc->n == 1);
        REQUIRE(Zdesc->mb == 64);
        REQUIRE(Zdesc->nb == 1);
        REQUIRE(Zdesc->bsiz == 64);
        REQUIRE(Zdesc->lmt*64 == 1024);
        REQUIRE(Zdesc->lnt == 1);
        REQUIRE(Zdesc->p == 1);
        REQUIRE(Zdesc->q == 1);
        REQUIRE(Zdesc->i == 1);
        REQUIRE(Zdesc->j == 1);
        REQUIRE(Zdesc->m0 == 0);
        REQUIRE(Zdesc->n0 == 0);
        REQUIRE(Zdesc->mt == 16);
        REQUIRE(Zdesc->nt == 1);
        REQUIRE(Zdesc->m1 == 1024);
        REQUIRE(Zdesc->n1 == 1);

        // Check the values set in the CHAM_desc_s struct for the Product descriptor
        auto* Pdesc = (CHAM_desc_t*) syntheticDataConfigurations->GetDescriptorProduct()[0];
        REQUIRE(Pdesc->m == 1024);
        REQUIRE(Pdesc->n == 1);
        REQUIRE(Pdesc->mb == 64);
        REQUIRE(Pdesc->nb == 1);
        REQUIRE(Pdesc->bsiz == 64);
        REQUIRE(Pdesc->lmt*64 == 1024);
        REQUIRE(Pdesc->lnt == 1);
        REQUIRE(Pdesc->p == 1);
        REQUIRE(Pdesc->q == 1);
        REQUIRE(Pdesc->i == 1);
        REQUIRE(Pdesc->j == 1);
        REQUIRE(Pdesc->m0 == 0);
        REQUIRE(Pdesc->n0 == 0);
        REQUIRE(Pdesc->mt == 16);
        REQUIRE(Pdesc->nt == 1);
        REQUIRE(Pdesc->m1 == 1024);
        REQUIRE(Pdesc->n1 == 1);
    }

    SECTION("Multiple"){
        syntheticDataConfigurations->SetProblemSize(2048);
        syntheticDataConfigurations->SetDenseTileSize(128);
        syntheticDataConfigurations->SetNumRHS(2);
        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorC()[1] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ()[0] != nullptr);
        REQUIRE(syntheticDataConfig
 */

