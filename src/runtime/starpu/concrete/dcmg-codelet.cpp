
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dcmg-codelet.cpp
 * @brief A class for starpu codelet dcmg.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-19
**/

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/dcmg-codelet.hpp>

using namespace exageostat::runtime;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;

template<typename T>
struct starpu_codelet DCMGCodelet<T>::cl_dcmg = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dcmg_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dcmg_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 1,
        .modes        = {STARPU_W},
        .name         = "dcmg"
};

template<typename T>
void DCMGCodelet<T>::InsertTask(void *apDescriptor, const int &aTriangularPart, Locations<T> *apLocation1,
                                Locations<T> *apLocation2, Locations<T> *apLocation3, T *apLocalTheta,
                                const int &aDistanceMetric, const Kernel<T> *apKernel) {
    int rows_num, cols_num, row, col, tile_row = 0, tile_col = 0;
    auto *CHAM_apDescriptor = (CHAM_desc_t *) apDescriptor;

    for (col = 0; col < CHAM_apDescriptor->nt; col++) {
        cols_num = col == CHAM_apDescriptor->nt - 1 ? CHAM_apDescriptor->n - col * CHAM_apDescriptor->nb
                                                : CHAM_apDescriptor->nb;
        if (aTriangularPart == ChamUpperLower) {
            row = 0;
        } else {
            row = CHAM_apDescriptor->m == CHAM_apDescriptor->n ? col : 0;
        }
        for (; row < CHAM_apDescriptor->mt; row++) {
            rows_num = row == CHAM_apDescriptor->mt - 1 ? CHAM_apDescriptor->m - row * CHAM_apDescriptor->mb
                                                  : CHAM_apDescriptor->mb;
            tile_row = row * CHAM_apDescriptor->mb;
            tile_col = col * CHAM_apDescriptor->nb;
            starpu_insert_task(&this->cl_dcmg,
                               STARPU_VALUE, &rows_num, sizeof(int),
                               STARPU_VALUE, &cols_num, sizeof(int),
                               STARPU_VALUE, &tile_row, sizeof(int),
                               STARPU_VALUE, &tile_col, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, row, col),
                               STARPU_VALUE, &apLocation1, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocation2, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocation3, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(kernels::Kernel<T> *),
                               0);
        }
    }
}

template<typename T>
void DCMGCodelet<T>::cl_dcmg_function(void *apBuffers[], void *apCodeletArguments) {
    int rows_num, cols_num, tile_row, tile_col, distance_metric;
    Locations<T> *pLocation1, *pLocation2, *pLocation3;
    T *pLocal_theta, *pDescriptor_A;
    Kernel<T> *pKernel;

    pDescriptor_A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    starpu_codelet_unpack_args(apCodeletArguments, &rows_num, &cols_num, &tile_row, &tile_col, &pLocation1, &pLocation2, &pLocation3, &pLocal_theta,
                               &distance_metric, &pKernel);
    pKernel->GenerateCovarianceMatrix(pDescriptor_A, rows_num, cols_num, tile_row, tile_col, *pLocation1, *pLocation2, *pLocation3, pLocal_theta, distance_metric);
}
