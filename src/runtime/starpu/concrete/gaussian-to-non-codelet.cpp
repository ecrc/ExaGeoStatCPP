
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file gaussian-to-non-codelet.cpp
 * @brief A class for starpu codelet gaussian-to-non.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#include <starpu.h>

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/concrete/gaussian-to-non-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet GaussianCodelet<T>::cl_gaussian_to_non = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_gaussian_to_non_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where        = STARPU_CPU,
        .cpu_funcs    = {cl_gaussian_to_non_function},
        .cuda_funcs={},
        .cuda_flags={0},
#endif
       .nbuffers    = 1,
        .modes        = {STARPU_RW},
        .name        = "gaussian_to_non"
};

template<typename T>
void GaussianCodelet<T>::InsertTask(void *apDesc, T *apTheta) {
    int row, rows_num;
    auto pDescriptor_Z = (CHAM_desc_t *) apDesc;

    for (row = 0; row < pDescriptor_Z->mt; row++) {
        rows_num = row == pDescriptor_Z->mt - 1 ? pDescriptor_Z->m - row * pDescriptor_Z->mb : pDescriptor_Z->mb;
        starpu_insert_task(&this->cl_gaussian_to_non,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr(pDescriptor_Z, row, 0),
                           STARPU_VALUE, &apTheta[0], sizeof(T),
                           STARPU_VALUE, &apTheta[1], sizeof(T),
                           STARPU_VALUE, &apTheta[2], sizeof(T),
                           STARPU_VALUE, &apTheta[3], sizeof(T),
                           STARPU_VALUE, &apTheta[4], sizeof(T),
                           STARPU_VALUE, &apTheta[5], sizeof(T),
                           0);
    }
}

template<typename T>
void GaussianCodelet<T>::cl_gaussian_to_non_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num;
    T *pDescriptorZ, *pTheta;

    pTheta = new T[6];
    pDescriptorZ = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num, &pTheta[0], &pTheta[1], &pTheta[2], &pTheta[3],
                               &pTheta[4],
                               &pTheta[5]);
    //core function to convert Z tile from Gaussian to non-Gaussian.
    core_gaussian_to_non(pDescriptorZ, pTheta, rows_num);
    delete[] pTheta;
}

template<typename T>
void GaussianCodelet<T>::core_gaussian_to_non(T *apDescriptorZ, const T *apLocalTheta, const int &aSize) {

    T xi = apLocalTheta[2];
    T omega = apLocalTheta[3];
    T g = apLocalTheta[4];
    T h = apLocalTheta[5];

    int i;
    if (h < 0) {
        throw std::runtime_error("The kurtosis parameter cannot be negative");
    }
    if (g == 0) {
        for (i = 0; i < aSize; i++)
            apDescriptorZ[i] = xi + omega * apDescriptorZ[i] * (exp(0.5 * h * pow(apDescriptorZ[i], 2)));
    } else {
        for (i = 0; i < aSize; i++)
            apDescriptorZ[i] =
                    xi + omega * (exp(g * apDescriptorZ[i]) - 1) * (exp(0.5 * h * pow(apDescriptorZ[i], 2))) / g;
    }
}


