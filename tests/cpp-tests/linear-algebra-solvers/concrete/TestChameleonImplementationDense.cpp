
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestChameleonImplmentationDense.cpp
 * @brief Unit tests for the Dense computation in the ExaGeoStat software package.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-04-06
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
using namespace exageostat::dataunits;

//Test that the function initializes the CHAM_descriptorC descriptor correctly.
void TEST_CHAMELEON_DESCRIPTORS_VALUES() {

    Configurations synthetic_data_configurations;
    int p = 1;
    SECTION("SINGLE")
    {

        // initialize Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 1, 0);

        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);
        linearAlgebraSolver->SetContext(hardware.GetChameleonContext());

        synthetic_data_configurations.SetProblemSize(4);
        synthetic_data_configurations.SetDenseTileSize(1);

        auto *data = new DescriptorData<float>();
        linearAlgebraSolver->InitiateDescriptors(synthetic_data_configurations, *data, p, nullptr);

        auto *CHAM_descriptorC = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        auto *CHAM_descriptorZ = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
        auto *CHAM_descriptorZcpy = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
        auto *CHAM_descriptorDeterminant = data->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                               DESCRIPTOR_DETERMINANT).chameleon_desc;
        auto *CHAM_descriptorProduct = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT).chameleon_desc;
        int N = synthetic_data_configurations.GetProblemSize();
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
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 4, 0);

        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);
        linearAlgebraSolver->SetContext(hardware.GetChameleonContext());

        synthetic_data_configurations.SetProblemSize(64);
        synthetic_data_configurations.SetDenseTileSize(8);

        auto *data = new DescriptorData<double>();
        linearAlgebraSolver->InitiateDescriptors(synthetic_data_configurations, *data, p, nullptr);

        auto *CHAM_descriptorC = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        auto *CHAM_descsubC11 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C11).chameleon_desc;
        auto *CHAM_descsubC12 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12).chameleon_desc;
        auto *CHAM_descsubC22 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22).chameleon_desc;
        auto *CHAM_descriptorZ = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
        auto *CHAM_descriptorZ_1 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_1).chameleon_desc;
        auto *CHAM_descriptorZ_2 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_2).chameleon_desc;
        auto *CHAM_descriptorZcpy = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
        auto *CHAM_descriptorDeterminant = data->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                               DESCRIPTOR_DETERMINANT).chameleon_desc;
        auto *CHAM_descriptorProduct = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT).chameleon_desc;
        auto *CHAM_descriptorProduct_1 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1).chameleon_desc;
        auto *CHAM_descriptorProduct_2 = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2).chameleon_desc;

        int N = synthetic_data_configurations.GetProblemSize();
        int dts = synthetic_data_configurations.GetDenseTileSize();
        int pGrid = synthetic_data_configurations.GetPGrid();
        int qGrid = synthetic_data_configurations.GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZ_1->m == N / 2);
        REQUIRE(CHAM_descriptorZ_2->m == N / 2);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        REQUIRE(CHAM_descsubC11->m == N / 2);
        REQUIRE(CHAM_descsubC12->m == N / 2);
        REQUIRE(CHAM_descsubC22->m == N / 2);
        REQUIRE(CHAM_descriptorProduct->m == 1);
        REQUIRE(CHAM_descriptorProduct_1->m == 1);
        REQUIRE(CHAM_descriptorProduct_2->m == 1);

        REQUIRE(CHAM_descriptorC->n == N);
        REQUIRE(CHAM_descriptorZ->n == 1);
        REQUIRE(CHAM_descriptorZ_1->n == 1);
        REQUIRE(CHAM_descriptorZ_2->n == 1);
        REQUIRE(CHAM_descriptorZcpy->n == 1);
        REQUIRE(CHAM_descriptorDeterminant->n == 1);
        REQUIRE(CHAM_descsubC11->n == N / 2);
        REQUIRE(CHAM_descsubC12->n == N / 2);
        REQUIRE(CHAM_descsubC22->n == N / 2);
        REQUIRE(CHAM_descriptorProduct->n == 1);
        REQUIRE(CHAM_descriptorProduct_1->n == 1);
        REQUIRE(CHAM_descriptorProduct_2->n == 1);

        REQUIRE(CHAM_descriptorC->mb == dts);
        REQUIRE(CHAM_descriptorZ->mb == dts);
        REQUIRE(CHAM_descriptorZ_1->mb == dts);
        REQUIRE(CHAM_descriptorZ_2->mb == dts);
        REQUIRE(CHAM_descriptorZcpy->mb == dts);
        REQUIRE(CHAM_descriptorDeterminant->mb == dts);
        REQUIRE(CHAM_descsubC11->mb == dts);
        REQUIRE(CHAM_descsubC12->mb == dts);
        REQUIRE(CHAM_descsubC22->mb == dts);
        REQUIRE(CHAM_descriptorProduct->mb == dts);
        REQUIRE(CHAM_descriptorProduct_1->mb == dts);
        REQUIRE(CHAM_descriptorProduct_2->mb == dts);

        REQUIRE(CHAM_descriptorC->nb == dts);
        REQUIRE(CHAM_descriptorZ->nb == dts);
        REQUIRE(CHAM_descriptorZ_1->nb == dts);
        REQUIRE(CHAM_descriptorZ_2->nb == dts);
        REQUIRE(CHAM_descriptorZcpy->nb == dts);
        REQUIRE(CHAM_descriptorDeterminant->nb == dts);
        REQUIRE(CHAM_descsubC11->nb == dts);
        REQUIRE(CHAM_descsubC12->nb == dts);
        REQUIRE(CHAM_descsubC22->nb == dts);
        REQUIRE(CHAM_descriptorProduct->nb == dts);
        REQUIRE(CHAM_descriptorProduct_1->nb == dts);
        REQUIRE(CHAM_descriptorProduct_2->nb == dts);

        REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ_1->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ_2->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZcpy->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorDeterminant->bsiz == dts * dts);
        REQUIRE(CHAM_descsubC11->bsiz == dts * dts);
        REQUIRE(CHAM_descsubC12->bsiz == dts * dts);
        REQUIRE(CHAM_descsubC22->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct_1->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct_2->bsiz == dts * dts);

        REQUIRE(CHAM_descriptorC->i == 0);
        REQUIRE(CHAM_descriptorZ->i == 0);
        REQUIRE(CHAM_descriptorZ_1->i == 0);
        REQUIRE(CHAM_descriptorZ_2->i == 0);
        REQUIRE(CHAM_descriptorZcpy->i == 0);
        REQUIRE(CHAM_descriptorDeterminant->i == 0);
        REQUIRE(CHAM_descsubC11->i == 0);
        REQUIRE(CHAM_descsubC12->i == N / 2);
        REQUIRE(CHAM_descsubC22->i == N / 2);
        REQUIRE(CHAM_descriptorProduct->i == 0);
        REQUIRE(CHAM_descriptorProduct_1->i == 0);
        REQUIRE(CHAM_descriptorProduct_2->i == 0);

        REQUIRE(CHAM_descriptorC->j == 0);
        REQUIRE(CHAM_descriptorZ->j == 0);
        REQUIRE(CHAM_descriptorZ_1->j == 0);
        REQUIRE(CHAM_descriptorZ_2->j == 0);
        REQUIRE(CHAM_descriptorZcpy->j == 0);
        REQUIRE(CHAM_descriptorDeterminant->j == 0);
        REQUIRE(CHAM_descsubC11->j == 0);
        REQUIRE(CHAM_descsubC12->j == 0);
        REQUIRE(CHAM_descsubC22->j == N / 2);
        REQUIRE(CHAM_descriptorProduct->j == 0);
        REQUIRE(CHAM_descriptorProduct_1->j == 0);
        REQUIRE(CHAM_descriptorProduct_2->j == 0);

        REQUIRE(CHAM_descriptorC->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ_1->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descriptorZ_2->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descriptorZcpy->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorDeterminant->mt == 1);
        REQUIRE(CHAM_descsubC11->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC12->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC22->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descriptorProduct->mt == 1);
        REQUIRE(CHAM_descriptorProduct_1->mt == 1);
        REQUIRE(CHAM_descriptorProduct_2->mt == 1);

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZ_1->nt == 1);
        REQUIRE(CHAM_descriptorZ_2->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        REQUIRE(CHAM_descsubC11->nt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC12->nt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC22->nt == ((N / 2) / dts));
        REQUIRE(CHAM_descriptorProduct->nt == 1);
        REQUIRE(CHAM_descriptorProduct_1->nt == 1);
        REQUIRE(CHAM_descriptorProduct_2->nt == 1);


        REQUIRE(CHAM_descriptorC->lm == N);
        REQUIRE(CHAM_descriptorZ->lm == N);
        REQUIRE(CHAM_descriptorZ_1->lm == N / 2);
        REQUIRE(CHAM_descriptorZ_2->lm == N / 2);
        REQUIRE(CHAM_descriptorZcpy->lm == N);
        REQUIRE(CHAM_descriptorDeterminant->lm == 1);
        REQUIRE(CHAM_descsubC11->lm == N);
        REQUIRE(CHAM_descsubC12->lm == N);
        REQUIRE(CHAM_descsubC22->lm == N);
        REQUIRE(CHAM_descriptorProduct->lm == 1);
        REQUIRE(CHAM_descriptorProduct_1->lm == 1);
        REQUIRE(CHAM_descriptorProduct_2->lm == 1);


        REQUIRE(CHAM_descriptorC->ln == N);
        REQUIRE(CHAM_descriptorZ->ln == 1);
        REQUIRE(CHAM_descriptorZ_1->ln == 1);
        REQUIRE(CHAM_descriptorZ_2->ln == 1);
        REQUIRE(CHAM_descriptorZcpy->ln == 1);
        REQUIRE(CHAM_descriptorDeterminant->ln == 1);
        REQUIRE(CHAM_descsubC11->ln == N);
        REQUIRE(CHAM_descsubC12->ln == N);
        REQUIRE(CHAM_descsubC22->ln == N);
        REQUIRE(CHAM_descriptorProduct->ln == 1);
        REQUIRE(CHAM_descriptorProduct_1->ln == 1);
        REQUIRE(CHAM_descriptorProduct_2->ln == 1);

        REQUIRE(CHAM_descriptorC->p == pGrid);
        REQUIRE(CHAM_descriptorZ->p == pGrid);
        REQUIRE(CHAM_descriptorZ_1->p == pGrid);
        REQUIRE(CHAM_descriptorZ_2->p == pGrid);
        REQUIRE(CHAM_descriptorZcpy->p == pGrid);
        REQUIRE(CHAM_descriptorDeterminant->p == pGrid);
        REQUIRE(CHAM_descsubC11->p == pGrid);
        REQUIRE(CHAM_descsubC12->p == pGrid);
        REQUIRE(CHAM_descsubC22->p == pGrid);
        REQUIRE(CHAM_descriptorProduct->p == pGrid);
        REQUIRE(CHAM_descriptorProduct_1->p == pGrid);
        REQUIRE(CHAM_descriptorProduct_2->p == pGrid);

        REQUIRE(CHAM_descriptorC->q == qGrid);
        REQUIRE(CHAM_descriptorZ->q == qGrid);
        REQUIRE(CHAM_descriptorZ_1->q == qGrid);
        REQUIRE(CHAM_descriptorZ_2->q == qGrid);
        REQUIRE(CHAM_descriptorZcpy->q == qGrid);
        REQUIRE(CHAM_descriptorDeterminant->q == qGrid);
        REQUIRE(CHAM_descsubC11->q == qGrid);
        REQUIRE(CHAM_descsubC12->q == qGrid);
        REQUIRE(CHAM_descsubC22->q == qGrid);
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
        mat = (double *) CHAM_descriptorZ->mat;
        for (auto i = 0; i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) *
                             (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
        }

        mat = (double *) CHAM_descriptorDeterminant->mat;
        for (auto i = 0; i < (CHAM_descriptorDeterminant->mt - 1) * (CHAM_descriptorDeterminant->nt - 1) *
                             (CHAM_descriptorDeterminant->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);

        }

        delete data;
    }
}

TEST_CASE("Chameleon Implementation Dense") {
    TEST_CHAMELEON_DESCRIPTORS_VALUES();

}
