
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ddotp-codelet.cpp
 * @brief A class for starpu codelet ddotp.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/ddotp-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet DDOTPCodelet<T>::cl_ddotp = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_ddotp_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_ddotp_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 2,
        .modes        = {STARPU_RW, STARPU_R},
        .name         = "ddotp"
};

template<typename T>
void DDOTPCodelet<T>::InsertTask(void *apDescA, void *apDescProduct) {

    int row, rows_num;
    auto pDesc_A = (CHAM_desc_t *) apDescA;
    auto pDesc_product = (CHAM_desc_t *) apDescProduct;

    auto desc_mt = pDesc_A->mt;
    auto desc_m = pDesc_A->m;
    auto desc_mb = pDesc_A->mb;

    for (row = 0; row < desc_mt; row++) {
        rows_num = row == desc_mt - 1 ? desc_m - row * desc_mb : desc_mb;

        starpu_insert_task(&this->cl_ddotp,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr(pDesc_product, 0, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr(pDesc_A, row, 0),
                           0);
    }
}

template<typename T>
void DDOTPCodelet<T>::cl_ddotp_function(void *apBuffers[], void *apCodeletArguments) {
    int rows_num;
    T *pDescriptor_A, *pDot_product;

    pDot_product = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pDescriptor_A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
    starpu_codelet_unpack_args(apCodeletArguments, &rows_num);
    T local_dot = cblas_ddot(rows_num, (double *) pDescriptor_A, 1, (double *) pDescriptor_A, 1);
    *pDot_product += local_dot;
}
