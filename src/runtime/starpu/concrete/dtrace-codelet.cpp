
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dtrace-codelet.cpp
 * @brief A class for starpu codelet dtrace.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/dtrace-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet DTRACECodelet<T>::cl_dtrace = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dtrace_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dtrace_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 3,
        .modes        = {STARPU_W, STARPU_W, STARPU_W},
        .name         = "dtrace"
};

template<typename T>
void DTRACECodelet<T>::InsertTask(void *apDescA, void *apDescNum, void *apDescTrace) {
    int row, rows_num;
    auto pDescriptor_A = (CHAM_desc_t *) apDescA;

    for (row = 0; row < pDescriptor_A->mt; row++) {
        rows_num = row == pDescriptor_A->mt - 1 ? pDescriptor_A->m - row * pDescriptor_A->mb : pDescriptor_A->mb;
        starpu_insert_task(&this->cl_dtrace,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescA, row, row),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescNum, 0, 0),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescTrace, row, 0),
                           0);
    }
}

template<typename T>
void DTRACECodelet<T>::cl_dtrace_function(void *apBuffers[], void *apCodeletArguments) {
    int rows_num;
    T *pDescriptor_A, *pSum, *pTrace;

    pDescriptor_A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pSum = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
    pTrace = (T *) STARPU_MATRIX_GET_PTR(apBuffers[2]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num);
    T local_sum = core_dtrace(pDescriptor_A, rows_num, pTrace);
    *pSum += local_sum;
}

template<typename T>
double DTRACECodelet<T>::core_dtrace(const T *pDescriptor, const int &aSize, T *pTrace) {
    T result = 0.0;
    for (int i = 0; i < aSize; i++) {
        result += pDescriptor[i + i * aSize];
        pTrace[i] = pDescriptor[i + i * aSize];
    }
    return result;
}
