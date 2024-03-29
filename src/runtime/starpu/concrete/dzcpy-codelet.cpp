
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dzcpy-codelet.cpp
 * @brief A class for starpu codelet dzcpy.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/dzcpy-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet DZCPYCodelet<T>::cl_dzcpy = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dzcpy_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dzcpy_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 1,
        .modes        = {STARPU_W},
        .name         = "dzcpy"
};

template<typename T>
void DZCPYCodelet<T>::InsertTask(void *apDescriptor, void *apDoubleVector) {
    int row, tile_row, rows_num;
    auto pDescriptor_A = (CHAM_desc_t *) apDescriptor;

    for (row = 0; row < pDescriptor_A->mt; row++) {
        rows_num = row == pDescriptor_A->mt - 1 ? pDescriptor_A->m - row * pDescriptor_A->mb : pDescriptor_A->mb;
        tile_row = row * pDescriptor_A->mb;
        starpu_insert_task(&this->cl_dzcpy,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_VALUE, &tile_row, sizeof(int),
                           STARPU_VALUE, &apDoubleVector, sizeof(double),
                           STARPU_W, RUNTIME_data_getaddr(pDescriptor_A, row, 0),
                           0);
    }
}

template<typename T>
void DZCPYCodelet<T>::cl_dzcpy_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num, tile_row;
    T *pDescriptor_A, *pR;

    pDescriptor_A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    starpu_codelet_unpack_args(apCodeletArguments, &rows_num, &tile_row, &pR);
    memcpy(pDescriptor_A, &pR[tile_row], rows_num * sizeof(T));
}
