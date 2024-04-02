
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file non-gaussian-codelet-transform-codelet.cpp
 * @brief A class for starpu codelet non-gaussian-transform.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-26
**/

#include <complex>
#include <starpu.h>

#include <runtime/starpu/concrete/non-gaussian-transform-codelet.hpp>

using namespace exageostat::runtime;

template<typename T>
struct starpu_codelet NonGaussianTransform<T>::cl_non_gaussian_transform = {
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_non_gaussian_transform_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_non_gaussian_transform_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers     = 1,
        .modes        = {STARPU_RW},
        .name         = "non_gaussian_transform"
};

template<typename T>
void
NonGaussianTransform<T>::InsertTask(void *apDescZ, const T *apTheta, std::unique_ptr<StarPuHelpers> &apStarPuHelpers) {
    int row, rows_num;

    auto desc_mt = apStarPuHelpers->GetMT(apDescZ);
    auto desc_m = apStarPuHelpers->GetM(apDescZ);
    auto desc_mb = apStarPuHelpers->GetMB(apDescZ);

    for (row = 0; row < desc_mt; row++) {
        rows_num = row == desc_mt - 1 ? desc_m - row * desc_mb : desc_mb;
        starpu_insert_task(&this->cl_non_gaussian_transform,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_RW, apStarPuHelpers->ExaGeoStatDataGetAddr(apDescZ, row, 0),
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
void NonGaussianTransform<T>::cl_non_gaussian_transform_function(void **apBuffers, void *apCodeletArguments) {
    int rows_num;
    T *pDescriptorZ, *pTheta;

    pTheta = new T[6];
    pDescriptorZ = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num, &pTheta[0], &pTheta[1], &pTheta[2], &pTheta[3],
                               &pTheta[4],
                               &pTheta[5]);
    core_non_gaussian_transform_helper(pDescriptorZ, pTheta, rows_num);
    delete[] pTheta;
}

template<typename T>
void
NonGaussianTransform<T>::core_non_gaussian_transform_helper(T *apDescripZ, const T *apLocalTheta, const int &aSize) {

    T xi = apLocalTheta[2];
    T omega = apLocalTheta[3];
    T g = apLocalTheta[4];
    T h = apLocalTheta[5];
    T eps = 1.0e-5;

    for (int i = 0; i < aSize; i++)
        apDescripZ[i] = newton_raphson(apDescripZ[i], xi, omega, g, h, eps);
}

template<typename T>
double
NonGaussianTransform<T>::newton_raphson(const T apDescriptorZ, const T aTransLocation, const T aTransScale,
                                        const T aTransShape, const T aTransKurtosis, const T aEpsilon) {
    int itr, max_itr;
    T x0 = 0, x1, all_err, diff;
    all_err = aEpsilon;
    max_itr = 1000;
    for (itr = 1; itr <= max_itr; itr++) {
        diff = tukeyGHTransfor(apDescriptorZ, x0, aTransLocation, aTransScale, aTransShape, aTransKurtosis) /
               tukeyGHDiferencial(x0, aTransScale, aTransShape, aTransKurtosis);
        x1 = x0 - diff;
        if (fabs(diff) < all_err)
            return x1;
        x0 = x1;
    }
    return x1;
}

template<typename T>
double NonGaussianTransform<T>::tukeyGHTransfor(const T aOriginalValue, const T aCurrentValue, const T aTransLocation,
                                                const T aTransScale, const T aTransShape, const T aTransKurtosis) {
    if (aTransShape == 0)
        return aOriginalValue - aTransLocation -
               aTransScale * aCurrentValue * exp(0.5 * aTransKurtosis * aCurrentValue * aCurrentValue);
    else
        return aOriginalValue - aTransLocation - (aTransScale * (exp(aTransShape * aCurrentValue) - 1) *
                                                  (exp(0.5 * aTransKurtosis * aCurrentValue * aCurrentValue)) /
                                                  aTransShape);
}

template<typename T>
double NonGaussianTransform<T>::tukeyGHDiferencial(const T aCurrentValue, const T aTransScale, const T aTransShape,
                                                   const T aTransKurtosis) {
    if (aTransShape == 0)
        return -aTransScale * exp((aTransKurtosis * aCurrentValue * aCurrentValue) / 2.0) -
               aTransScale * aTransKurtosis * aCurrentValue * aCurrentValue *
               exp((aTransKurtosis * aCurrentValue * aCurrentValue) / 2.0);
    else
        return -aTransScale * exp(aTransShape * aCurrentValue) *
               exp((aTransKurtosis * aCurrentValue * aCurrentValue) / 2.0) -
               (aTransKurtosis * aCurrentValue * exp((aTransKurtosis * aCurrentValue * aCurrentValue) / 2.0) *
                (aTransScale * exp(aTransShape * aCurrentValue) - aTransScale)) / aTransShape;
}
