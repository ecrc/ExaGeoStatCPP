
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmse-codelet.cpp
 * @brief A class for starpu codelet dmse.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-21
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/dmse-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet DMSECodelet<T>::cl_dmse = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dmse_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dmse_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers    = 3,
        .modes        = {STARPU_W, STARPU_W, STARPU_W},
        .name        = "dmse"
};

template<typename T>
void DMSECodelet<T>::InsertTask(void *apDescError, void *apDescZPredict, void *apDescZMiss) {
    int row, rows_num;
    auto pDesc_Z_predict = (CHAM_desc_t *) apDescZPredict;

    for (row = 0; row < pDesc_Z_predict->mt; row++) {
        rows_num =
                row == pDesc_Z_predict->mt - 1 ? pDesc_Z_predict->m - row * pDesc_Z_predict->mb : pDesc_Z_predict->mb;
        starpu_insert_task(&this->cl_dmse,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescError, 0, 0),
                           STARPU_W,
                           (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZPredict, row, 0),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZMiss, row, 0),
                           0);
    }
}

template<typename T>
void DMSECodelet<T>::cl_dmse_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num;
    T *pZPredict, *pZMiss, *pError;
    T local_error = 0.0;

    pError = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pZPredict = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
    pZMiss = (T *) STARPU_MATRIX_GET_PTR(apBuffers[2]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num);
    for (int i = 0; i < rows_num; i++) {
        local_error += pow((pZPredict[i] - pZMiss[i]), 2);
    }
    *pError += local_error;
}
