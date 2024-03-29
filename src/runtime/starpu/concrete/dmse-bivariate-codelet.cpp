
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmse-bivariate-codelet.cpp
 * @brief A class for starpu codelet dmse-bivariate.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/dmse-bivariate-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet DMSEBivariateCodelet<T>::cl_dmse_bivariate = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dmse_bivariate_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dmse_bivariate_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 5,
        .modes        = {STARPU_RW, STARPU_RW, STARPU_RW, STARPU_R, STARPU_R},
        .name         = "dmse-bivariate"
};

template<typename T>
void DMSEBivariateCodelet<T>::InsertTask(void *apDescZMiss, void *apDescZPre, void *apDescsError, void *apDescsError1,
                                         void *apDescsError2) {
    int row, rows_num;

    auto pDesc_ZPre = (CHAM_desc_t *) apDescZPre;
    auto desc_mt = pDesc_ZPre->mt;
    auto desc_m = pDesc_ZPre->m;
    auto desc_mb = pDesc_ZPre->mb;

    for (row = 0; row < desc_mt; row++) {
        rows_num = row == desc_mt - 1 ? desc_m - row * desc_mb : desc_mb;
        starpu_insert_task(&this->cl_dmse_bivariate,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescsError1, 0, 0),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescsError2, 0, 0),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescsError, 0, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZPre, row, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZMiss, row, 0),
                           0);
    }
}

template<typename T>
void DMSEBivariateCodelet<T>::cl_dmse_bivariate_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num;
    T *pZpre, *pZmiss, *pSerror1, *pSerror2, *pSerror;
    T local_serror1 = 0.0, local_serror2 = 0.0, local_serror = 0.0;

    pSerror1 = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pSerror2 = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
    pSerror = (T *) STARPU_MATRIX_GET_PTR(apBuffers[2]);
    pZpre = (T *) STARPU_MATRIX_GET_PTR(apBuffers[3]);
    pZmiss = (T *) STARPU_MATRIX_GET_PTR(apBuffers[4]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num);

    for (int i = 0; i < rows_num; i++) {
        if (i % 2 == 0) {
            local_serror1 += pow((pZpre[i] - pZmiss[i]), 2);
        } else
            local_serror2 += pow((pZpre[i] - pZmiss[i]), 2);
        local_serror += pow((pZpre[i] - pZmiss[i]), 2);
    }
    *pSerror1 += local_serror1;
    *pSerror2 += local_serror2;
    *pSerror += local_serror;
}
