
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmloe-mmom-codelet.cpp
 * @brief A class for starpu codelet dmloe-mmom.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/dmloe-mmom-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet DmloeMmomCodelet<T>::cl_dmloe_mmom = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dmloe_mmom_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dmloe_mmom_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers    = 5,
        .modes        = {STARPU_R, STARPU_R, STARPU_R, STARPU_RW, STARPU_RW},
        .name        = "dmloe_mmom"
};

template<typename T>
void DmloeMmomCodelet<T>::InsertTask(void *apDescExpr1, void *apDescExpr2, void *apDescExpr3, void *apDescMLOE,
                                     void *apDescMMOM) {
    int row, col, rows_num, cols_num;

    for (col = 0; col < ((CHAM_desc_t *) apDescExpr1)->nt; col++) {
        cols_num = col == ((CHAM_desc_t *) apDescExpr1)->nt - 1 ? ((CHAM_desc_t *) apDescExpr1)->n -
                                                                  col * ((CHAM_desc_t *) apDescExpr1)->nb
                                                                : ((CHAM_desc_t *) apDescExpr1)->nb;
        for (row = 0; row < ((CHAM_desc_t *) apDescExpr1)->mt; row++) {

            rows_num = row == ((CHAM_desc_t *) apDescExpr1)->mt - 1 ? ((CHAM_desc_t *) apDescExpr1)->m -
                                                                      row * ((CHAM_desc_t *) apDescExpr1)->mb
                                                                    : ((CHAM_desc_t *) apDescExpr1)->mb;
            starpu_insert_task(&this->cl_dmloe_mmom,
                               STARPU_VALUE, &rows_num, sizeof(int),
                               STARPU_VALUE, &cols_num, sizeof(int),
                               STARPU_R,
                               (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescExpr1, row, col),
                               STARPU_R,
                               (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescExpr2, row, col),
                               STARPU_R,
                               (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescExpr3, row, col),
                               STARPU_RW,
                               (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescMLOE, row, col),
                               STARPU_RW,
                               (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescMMOM, row, col),
                               0);
        }
    }
}

template<typename T>
void DmloeMmomCodelet<T>::cl_dmloe_mmom_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num, cols_num;
    T *pExpr1, *pExpr2, *pExpr3, *pMloe, *pMmom;

    pExpr1 = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pExpr2 = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
    pExpr3 = (T *) STARPU_MATRIX_GET_PTR(apBuffers[2]);
    pMloe = (T *) STARPU_MATRIX_GET_PTR(apBuffers[3]);
    pMmom = (T *) STARPU_MATRIX_GET_PTR(apBuffers[4]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num, &cols_num);
    T expr1_ = 0, expr2_ = 0, expr3_ = 0;

    for (int i = 0; i < rows_num * cols_num; i += 2) {
        expr1_ += pExpr1[i];
        expr2_ += pExpr2[i];
        expr3_ += pExpr3[i];
    }

    if (expr2_ == 0.0) {
        *pMloe -= 1.0;
    } else {
        *pMloe += (expr1_ / expr2_) - 1.0;
    }

    if (expr2_ == 0.0) {
        *pMmom -= 1.0;
    } else {
        *pMmom += (expr3_ / expr1_) - 1.0;
    }
}
