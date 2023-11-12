
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestChameleonImplementationDST.cpp
 * @brief Unit tests for the Diagonal Super Tile computation in the ExaGeoStat software package.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-04-09
**/

extern "C" {
#include <control/context.h>
}

#include <catch2/catch_all.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <configurations/Configurations.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <data-units/DescriptorData.hpp>

using namespace std;

using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;
using namespace exageostat::hardware;


//Test that the function initializes the CHAM_descriptorC descriptor correctly.
void TEST_CHAMELEON_DESCRIPTORS_VALUES_DST() {

    Configurations synthetic_data_configurations;

    SECTION("SINGLE")
    {
        // initialize Hardware.
        auto hardware = ExaGeoStatHardware(DIAGONAL_APPROX, 1, 0);

        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(DIAGONAL_APPROX);
        linearAlgebraSolver->SetContext(hardware.GetChameleonContext());

        synthetic_data_configurations.SetProblemSize(64);
        synthetic_data_configurations.SetDenseTileSize(16);

        auto *data = new DescriptorData<float>();
        linearAlgebraSolver->InitiateDescriptors(synthetic_data_configurations, *data, nullptr);

        auto *CHAM_descriptorC = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        auto *CHAM_descriptorZ = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
        auto *CHAM_descriptorZcpy = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
        auto *CHAM_descriptorDeterminant = data->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                               DESCRIPTOR_DETERMINANT).chameleon_desc;
        auto *CHAM_descriptorProduct = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT).chameleon_desc;
        int N = synthetic_data_configurations.GetProblemSize() * synthetic_data_configurations.GetP();
        int dts = synthetic_data_configurations.GetDenseTileSize();
        int pGrid = synthetic_data_configurations.GetPGrid();
        int qGrid = synthetic_data_configurations.GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        REQUIRE(CHAM_descriptorProduct->m == 1);

        REQUIRE(CHAM_descriptorC->n == N);
        REQUIRE(CHAM_descriptorZ->n == 1);
        REQUIRE(CHAM_descriptorZcpy->n == 1);
        REQUIRE(CHAM_descriptorDeterminant->n == 1);
        REQUIRE(CHAM_descriptorProduct->n == 1);

        REQUIRE(CHAM_descriptorC->mb == dts);
        REQUIRE(CHAM_descriptorZ->mb == dts);
        REQUIRE(CHAM_descriptorZcpy->mb == dts);
        REQUIRE(CHAM_descriptorDeterminant->mb == dts);
        REQUIRE(CHAM_descriptorProduct->mb == dts);

        REQUIRE(CHAM_descriptorC->nb == dts);
        REQUIRE(CHAM_descriptorZ->nb == dts);
        REQUIRE(CHAM_descriptorZcpy->nb == dts);
        REQUIRE(CHAM_descriptorDeterminant->nb == dts);
        REQUIRE(CHAM_descriptorProduct->nb == dts);

        REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZcpy->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorDeterminant->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct->bsiz == dts * dts);

        REQUIRE(CHAM_descriptorC->i == 0);
        REQUIRE(CHAM_descriptorZ->i == 0);
        REQUIRE(CHAM_descriptorZcpy->i == 0);
        REQUIRE(CHAM_descriptorDeterminant->i == 0);
        REQUIRE(CHAM_descriptorProduct->i == 0);

        REQUIRE(CHAM_descriptorC->j == 0);
        REQUIRE(CHAM_descriptorZ->j == 0);
        REQUIRE(CHAM_descriptorZcpy->j == 0);
        REQUIRE(CHAM_descriptorDeterminant->j == 0);
        REQUIRE(CHAM_descriptorProduct->j == 0);

        REQUIRE(CHAM_descriptorC->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZcpy->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorDeterminant->mt == 1);
        REQUIRE(CHAM_descriptorProduct->mt == 1);

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        REQUIRE(CHAM_descriptorProduct->nt == 1);

        REQUIRE(CHAM_descriptorC->lm == N);
        REQUIRE(CHAM_descriptorZ->lm == N);
        REQUIRE(CHAM_descriptorZcpy->lm == N);
        REQUIRE(CHAM_descriptorDeterminant->lm == 1);
        REQUIRE(CHAM_descriptorProduct->lm == 1);

        REQUIRE(CHAM_descriptorC->ln == N);
        REQUIRE(CHAM_descriptorZ->ln == 1);
        REQUIRE(CHAM_descriptorZcpy->ln == 1);
        REQUIRE(CHAM_descriptorDeterminant->ln == 1);
        REQUIRE(CHAM_descriptorProduct->ln == 1);

        REQUIRE(CHAM_descriptorC->p == pGrid);
        REQUIRE(CHAM_descriptorZ->p == pGrid);
        REQUIRE(CHAM_descriptorZcpy->p == pGrid);
        REQUIRE(CHAM_descriptorDeterminant->p == pGrid);
        REQUIRE(CHAM_descriptorProduct->p == pGrid);

        REQUIRE(CHAM_descriptorC->q == qGrid);
        REQUIRE(CHAM_descriptorZ->q == qGrid);
        REQUIRE(CHAM_descriptorZcpy->q == qGrid);
        REQUIRE(CHAM_descriptorDeterminant->q == qGrid);
        REQUIRE(CHAM_descriptorProduct->q == qGrid);


        auto *mat = (float *) CHAM_descriptorZ->mat;
        auto *matZcpy = (float *) CHAM_descriptorZcpy->mat;
        for (auto i = 0;
             i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) * (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matZcpy[i] == 0.0f);
        }

        mat = (float *) CHAM_descriptorDeterminant->mat;
        auto *matProduct = (float *) CHAM_descriptorProduct->mat;
        for (auto i = 0; i < (CHAM_descriptorDeterminant->mt - 1) * (CHAM_descriptorDeterminant->nt - 1) *
                             (CHAM_descriptorDeterminant->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matProduct[i] == 0.0f);
        }

        delete data;
    }

    SECTION("DOUBLE")
    {

        // initialize Hardware.
        auto hardware = ExaGeoStatHardware(DIAGONAL_APPROX, 2, 0);

        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(DIAGONAL_APPROX);
        linearAlgebraSolver->SetContext(hardware.GetChameleonContext());

        synthetic_data_configurations.SetProblemSize(32);
        synthetic_data_configurations.SetDenseTileSize(16);

        auto *data = new DescriptorData<double>();
        linearAlgebraSolver->InitiateDescriptors(synthetic_data_configurations, *data, nullptr);

        auto *CHAM_descriptorC = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        auto *CHAM_descriptorZ = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
        auto *CHAM_descriptorZcpy = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
        auto *CHAM_descriptorDeterminant = data->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                               DESCRIPTOR_DETERMINANT).chameleon_desc;
        auto *CHAM_descriptorProduct = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT).chameleon_desc;
        auto *CHAM_descriptorProduct_1 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1).chameleon_desc;
        auto *CHAM_descriptorProduct_2 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2).chameleon_desc;
        int N = synthetic_data_configurations.GetProblemSize() * synthetic_data_configurations.GetP();
        int dts = synthetic_data_configurations.GetDenseTileSize();
        int pGrid = synthetic_data_configurations.GetPGrid();
        int qGrid = synthetic_data_configurations.GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        REQUIRE(CHAM_descriptorProduct->m == 1);
        REQUIRE(CHAM_descriptorProduct_1->m == 1);
        REQUIRE(CHAM_descriptorProduct_2->m == 1);

        REQUIRE(CHAM_descriptorC->n == N);
        REQUIRE(CHAM_descriptorZ->n == 1);
        REQUIRE(CHAM_descriptorZcpy->n == 1);
        REQUIRE(CHAM_descriptorDeterminant->n == 1);
        REQUIRE(CHAM_descriptorProduct->n == 1);
        REQUIRE(CHAM_descriptorProduct_1->n == 1);
        REQUIRE(CHAM_descriptorProduct_2->n == 1);

        REQUIRE(CHAM_descriptorC->mb == dts);
        REQUIRE(CHAM_descriptorZ->mb == dts);
        REQUIRE(CHAM_descriptorZcpy->mb == dts);
        REQUIRE(CHAM_descriptorDeterminant->mb == dts);
        REQUIRE(CHAM_descriptorProduct->mb == dts);
        REQUIRE(CHAM_descriptorProduct_1->mb == dts);
        REQUIRE(CHAM_descriptorProduct_2->mb == dts);

        REQUIRE(CHAM_descriptorC->nb == dts);
        REQUIRE(CHAM_descriptorZ->nb == dts);
        REQUIRE(CHAM_descriptorZcpy->nb == dts);
        REQUIRE(CHAM_descriptorDeterminant->nb == dts);
        REQUIRE(CHAM_descriptorProduct->nb == dts);
        REQUIRE(CHAM_descriptorProduct_1->nb == dts);
        REQUIRE(CHAM_descriptorProduct_2->nb == dts);

        REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZcpy->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorDeterminant->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct_1->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct_2->bsiz == dts * dts);

        REQUIRE(CHAM_descriptorC->i == 0);
        REQUIRE(CHAM_descriptorZ->i == 0);
        REQUIRE(CHAM_descriptorZcpy->i == 0);
        REQUIRE(CHAM_descriptorDeterminant->i == 0);
        REQUIRE(CHAM_descriptorProduct->i == 0);
        REQUIRE(CHAM_descriptorProduct_1->i == 0);
        REQUIRE(CHAM_descriptorProduct_2->i == 0);


        REQUIRE(CHAM_descriptorC->j == 0);
        REQUIRE(CHAM_descriptorZ->j == 0);
        REQUIRE(CHAM_descriptorZcpy->j == 0);
        REQUIRE(CHAM_descriptorDeterminant->j == 0);
        REQUIRE(CHAM_descriptorProduct->j == 0);
        REQUIRE(CHAM_descriptorProduct_1->j == 0);
        REQUIRE(CHAM_descriptorProduct_2->j == 0);

        REQUIRE(CHAM_descriptorC->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZcpy->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorDeterminant->mt == 1);
        REQUIRE(CHAM_descriptorProduct->mt == 1);
        REQUIRE(CHAM_descriptorProduct_1->mt == 1);
        REQUIRE(CHAM_descriptorProduct_2->mt == 1);

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        REQUIRE(CHAM_descriptorProduct->nt == 1);
        REQUIRE(CHAM_descriptorProduct_1->nt == 1);
        REQUIRE(CHAM_descriptorProduct_2->nt == 1);

        REQUIRE(CHAM_descriptorC->lm == N);
        REQUIRE(CHAM_descriptorZ->lm == N);
        REQUIRE(CHAM_descriptorZcpy->lm == N);
        REQUIRE(CHAM_descriptorDeterminant->lm == 1);
        REQUIRE(CHAM_descriptorProduct->lm == 1);
        REQUIRE(CHAM_descriptorProduct_1->lm == 1);
        REQUIRE(CHAM_descriptorProduct_2->lm == 1);

        REQUIRE(CHAM_descriptorC->ln == N);
        REQUIRE(CHAM_descriptorZ->ln == 1);
        REQUIRE(CHAM_descriptorZcpy->ln == 1);
        REQUIRE(CHAM_descriptorDeterminant->ln == 1);
        REQUIRE(CHAM_descriptorProduct->ln == 1);
        REQUIRE(CHAM_descriptorProduct_1->ln == 1);
        REQUIRE(CHAM_descriptorProduct_2->ln == 1);

        REQUIRE(CHAM_descriptorC->p == pGrid);
        REQUIRE(CHAM_descriptorZ->p == pGrid);
        REQUIRE(CHAM_descriptorZcpy->p == pGrid);
        REQUIRE(CHAM_descriptorDeterminant->p == pGrid);
        REQUIRE(CHAM_descriptorProduct->p == pGrid);
        REQUIRE(CHAM_descriptorProduct_1->p == pGrid);
        REQUIRE(CHAM_descriptorProduct_2->p == pGrid);

        REQUIRE(CHAM_descriptorC->q == qGrid);
        REQUIRE(CHAM_descriptorZ->q == qGrid);
        REQUIRE(CHAM_descriptorZcpy->q == qGrid);
        REQUIRE(CHAM_descriptorDeterminant->q == qGrid);
        REQUIRE(CHAM_descriptorProduct->q == qGrid);
        REQUIRE(CHAM_descriptorProduct_1->q == qGrid);
        REQUIRE(CHAM_descriptorProduct_2->q == qGrid);

        auto *mat = (double *) CHAM_descriptorZ->mat;
        auto *matZcpy = (double *) CHAM_descriptorZcpy->mat;
        for (auto i = 0;
             i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) * (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matZcpy[i] == 0.0f);
        }
        mat = (double *) CHAM_descriptorDeterminant->mat;
        for (auto i = 0; i < (CHAM_descriptorDeterminant->mt - 1) * (CHAM_descriptorDeterminant->nt - 1) *
                             (CHAM_descriptorDeterminant->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
        }

        delete data;
    }
}

TEST_CASE("Chameleon Implementation DST") {
    TEST_CHAMELEON_DESCRIPTORS_VALUES_DST();

}