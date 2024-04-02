
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file stride-vec-codelet.cpp
 * @brief A class for starpu codelet stride-vec.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/stride-vec-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet STRIDEVECCodelet<T>::cl_stride_vec = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_stride_vec_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_stride_vec_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 3,
        .modes        = {STARPU_R, STARPU_W, STARPU_W},
        .name         = "stride_vec"
};

template<typename T>
void STRIDEVECCodelet<T>::InsertTask(const void *apDescA, void *apDescB, void *apDescC) {
    int row, rows_num;
    const auto desc_A = (CHAM_desc_t *) apDescA;
    auto desc_B = (CHAM_desc_t *) apDescB;
    auto desc_C = (CHAM_desc_t *) apDescC;

    auto desc_mt = desc_A->mt;
    auto desc_m = desc_A->m;
    auto desc_mb = desc_A->mb;

    for (row = 0; row < desc_mt; row++) {
        rows_num = row == desc_mt - 1 ? desc_m - row * desc_mb : desc_mb;
        starpu_insert_task(&this->cl_stride_vec,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr(desc_A, row, 0),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(desc_B, (int) floor(row / 2.0), 0),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(desc_C, (int) floor(row / 2.0), 0),
                           0);
    }
}

template<typename T>
void STRIDEVECCodelet<T>::cl_stride_vec_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num;
    T *pDescriptor_A, *pDescriptor_B, *pDescriptor_C;

    pDescriptor_A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pDescriptor_B = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
    pDescriptor_C = (T *) STARPU_MATRIX_GET_PTR(apBuffers[2]);
    starpu_codelet_unpack_args(apCodeletArguments, &rows_num);

    for (int i = 0, j = 0; i < rows_num - 1; i += 2, j++) {
        pDescriptor_B[j] = pDescriptor_A[i];
        pDescriptor_C[j] = pDescriptor_A[i + 1];
    }
}
