
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmdet-codelet.cpp
 * @brief A class for starpu codelet dmdet.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-21
**/

#include <starpu.h>

#include <runtime/starpu/concrete/dmdet-codelet.hpp>

using namespace exageostat::runtime;
using namespace exageostat::common;

template<typename T>
struct starpu_codelet DMDETCodelet<T>::cl_dmdet{
#ifdef USE_CUDA
        .where= STARPU_CPU | STARPU_CUDA,
        .cpu_funcs={cl_dmdet_function},
        .cuda_funcs={},
        .cuda_flags={0},
#else
        .where=STARPU_CPU,
        .cpu_funcs={cl_dmdet_function},
        .cuda_funcs={},
        .cuda_flags={(0)},
#endif
        .nbuffers    = 2,
        .modes        = {STARPU_R, STARPU_RW},
        .name        = "dmdet"
};

template<typename T>
void DMDETCodelet<T>::InsertTask(const Computation &aComputation, void *apDescA, void *apDescDet,
                                 std::unique_ptr<StarPuHelpers> &aStarPuHelpers) {
    int row, rows_num;
    auto desc_mt = aStarPuHelpers->GetMT(apDescA);
    auto desc_m = aStarPuHelpers->GetM(apDescA);
    auto desc_mb = aStarPuHelpers->GetMB(apDescA);

    for (row = 0; row < desc_mt; row++) {
        rows_num = row == desc_mt - 1 ? desc_m - row * desc_mb : desc_mb;
        starpu_insert_task(&this->cl_dmdet,
                           STARPU_VALUE, &rows_num, sizeof(int),
                           STARPU_R,
                           aStarPuHelpers->ExaGeoStatDataGetAddr(apDescA, row, aComputation != TILE_LOW_RANK ? row : 0),
                           STARPU_RW, aStarPuHelpers->ExaGeoStatDataGetAddr(apDescDet, 0, 0),
                           0);
    }
}

template<typename T>
void DMDETCodelet<T>::cl_dmdet_function(void *apBuffers[], void *apCodeletArguments) {
    int rows_num;
    T *pDescriptor_A, *pDeterminant;

    pDescriptor_A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
    pDeterminant = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);

    starpu_codelet_unpack_args(apCodeletArguments, &rows_num);
    T local_det = core_dmdet(pDescriptor_A, rows_num);
    *pDeterminant += local_det;
}

template<typename T>
T DMDETCodelet<T>::core_dmdet(const T *apDescriptor, const int &aSize) {
    T result = 0.0;
    for (int i = 0; i < aSize; i++) {
        if (apDescriptor[i + i * aSize] > 0)
            result += log(apDescriptor[i + i * aSize]);
    }
    return result;
}
