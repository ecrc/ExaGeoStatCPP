// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestDescriptorData.cpp
 * @brief Tests for CHAMELEON to HICMA descriptor conversion in ExaGeoStat.
 * @details This test case verifies the conversion of matrix descriptors from the CHAMELEON format to the HICMA format.
 * It ensures that key properties of the matrix descriptor, such as dimensions, block sizes, and grid distribution parameters, are preserved during the conversion process.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <catch2/catch_all.hpp>

#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace std;

using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;

void TEST_CHAM_TO_HICMA_CONV() {

    // Initialize Configuration
    Configurations synthetic_data_configurations;
    synthetic_data_configurations.SetProblemSize(4);
    synthetic_data_configurations.SetDenseTileSize(1);

    // Initialize linear algebra solver
    int p = 1;
    auto hardware = ExaGeoStatHardware(EXACT_DENSE, 1, 0);
    auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);
    auto *data = new DescriptorData<float>();
    linearAlgebraSolver->InitiateDescriptors(synthetic_data_configurations, *data, p);

    // Create CHAM descriptor and convert it to HICMA descriptor
    auto *CHAM_descriptorC = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
    auto *HICMA_descriptor = data->ConvertChameleonToHicma(CHAM_descriptorC);

    // Verify common attributes are of same value
    REQUIRE(CHAM_descriptorC->m == HICMA_descriptor->m);
    REQUIRE(CHAM_descriptorC->n == HICMA_descriptor->n);
    REQUIRE(CHAM_descriptorC->mb == HICMA_descriptor->mb);
    REQUIRE(CHAM_descriptorC->nb == HICMA_descriptor->nb);
    REQUIRE(CHAM_descriptorC->bsiz == HICMA_descriptor->bsiz);
    REQUIRE(CHAM_descriptorC->i == HICMA_descriptor->i);
    REQUIRE(CHAM_descriptorC->j == HICMA_descriptor->j);
    REQUIRE(CHAM_descriptorC->mt == HICMA_descriptor->mt);
    REQUIRE(CHAM_descriptorC->nt == HICMA_descriptor->nt);
    REQUIRE(CHAM_descriptorC->lm == HICMA_descriptor->lm);
    REQUIRE(CHAM_descriptorC->ln == HICMA_descriptor->ln);
    REQUIRE(CHAM_descriptorC->p == HICMA_descriptor->p);
    REQUIRE(CHAM_descriptorC->q == HICMA_descriptor->q);

    delete data;
}

TEST_CASE(" CHAMELEON To HICMA Converter") {
    TEST_CHAM_TO_HICMA_CONV();
}
