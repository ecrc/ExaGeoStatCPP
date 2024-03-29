
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file non-gaussian-codelet-loglike-codelet.cpp
 * @brief A class for starpu codelet non-gaussian-loglike.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-26
**/

#include <starpu.h>

#include <runtime/starpu/concrete/non-gaussian-loglike-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet NonGaussianLoglike<T>::cl_non_gaussian_loglike = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_non_gaussian_loglike_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_non_gaussian_loglike_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 2,
        .modes        = {STARPU_R, STARPU_RW},
        .name         = "non_gaussian_loglike"
};

template<typename T>
void NonGaussianLoglike<T>::InsertTask(void *apDescZ, void *apDescSum, const T *apTheta,
                                       std::unique_ptr<StarPuHelpers> &aStarPuHelpers) {
    auto desc_mt = aStarPuHelpers->GetMT(apDescZ);
    auto desc_m = aStarPuHelpers->GetM(apDescZ);
    auto desc_mb = aStarPuHelpers->GetMB(apDescZ);

    int row, rows_num;
    for (row = 0; row < desc_mt; row++) {
        rows_num = row == desc_mt - 1 ? desc_m - row * desc_mb : desc_mb;
        starpu_insert_task(&this->cl_non_gaussian_loglike,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_R, aStarPuHelpers->ExaGeoStatDataGetAddr(apDescZ, row, 0),
                           STARPU_RW, aStarPuHelpers->ExaGeoStatDataGetAddr(apDescSum, 0, 0),
                           STARPU_VALUE, &apTheta[0], sizeof(double),
                           STARPU_VALUE, &apTheta[1], sizeof(double),
                           STARPU_VALUE, &apTheta[2], sizeof(double),
                           STARPU_VALUE, &apTheta[3], sizeof(double),
                           STARPU_VALUE, &apTheta[4], sizeof(double),
                           STARPU_VALUE, &apTheta[5], sizeof(double),
                           0);
    }
}

template<typename T>
void NonGaussianLoglike<T>::cl_non_gaussian_loglike_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num;
    T *pDescriptor_Z, *pDescriptor_sum;

    auto *pTheta = new T[6];
    pDescriptor_Z = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pDescriptor_sum = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num, &pTheta[0], &pTheta[1], &pTheta[2], &pTheta[3],
                               &pTheta[4],
                               &pTheta[5]);
    T local_sum = core_non_gaussian_loglike_helper(pDescriptor_Z, pTheta, rows_num);
    *pDescriptor_sum += local_sum;
    delete[] pTheta;
}

template<typename T>
double NonGaussianLoglike<T>::core_non_gaussian_loglike_helper(const T *apDescriptorZ, const T *apLocalTheta, const int &aSize) {
    T g = apLocalTheta[4];
    T h = apLocalTheta[5];

    int i;
    T sum = 0.0;
    if (h < 0) {
        throw std::runtime_error("The kurtosis parameter cannot be negative");

    }
    for (i = 0; i < aSize; i++) {
        if (g == 0)
            sum += log(1 + h * pow(apDescriptorZ[i], 2)) + 0.5 * h * pow(apDescriptorZ[i], 2);
        else {
            sum += log(exp(g * apDescriptorZ[i]) + (exp(g * apDescriptorZ[i]) - 1) * h * apDescriptorZ[i] / g) + 0.5 * h * pow(apDescriptorZ[i], 2);
        }
    }
    return sum;
}
