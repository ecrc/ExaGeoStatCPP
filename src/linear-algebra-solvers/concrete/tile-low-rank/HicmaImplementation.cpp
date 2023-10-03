
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.cpp
 * @brief Sets up the HiCMA descriptors needed for the tile low rank computations in ExaGeoStat.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#include <linear-algebra-solvers/concrete/hicma/tile-low-rank/HicmaImplementation.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::hardware;
using namespace exageostat::configurations;


template<typename T>
T HicmaImplementation<T>::ExaGeoStatMLETile(const ExaGeoStatHardware &apHardware, ExaGeoStatData<T> &aData,
                                            Configurations &aConfigurations, const double *theta,
                                            T *apMeasurementsMatrix) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatLapackCopyTile(const UpperLower &aUpperLower, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");

}

template<typename T>
void
HicmaImplementation<T>::ExaGeoStatOptionsInit(void *apOptoins, void *apContext, void *apSequence,
                                              void *apRequest) {

    HICMA_RUNTIME_options_init((HICMA_option_t *) apOptoins, (HICMA_context_t *) apContext,
                               (HICMA_sequence_t *) apSequence,
                               (HICMA_request_t *) apRequest);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatOptionsFree(void *apOptions) {
    HICMA_RUNTIME_options_ws_free((HICMA_option_t *) apOptions);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatSequenceWait(void *apSequence) {
    HICMA_Sequence_Wait((HICMA_sequence_t *) apSequence);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatCreateSequence(void *apSequence) {
    int status = HICMA_Sequence_Create((HICMA_sequence_t **) apSequence);
    if (status != HICMA_SUCCESS) {
        throw std::runtime_error("HICMA_Sequence_Create Failed!");
    }
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) {
    HICMA_RUNTIME_options_finalize((HICMA_option_t *) apOptions, (HICMA_context_t *) apContext);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatPotrfTile(const UpperLower &aUpperLower, void *apA, int aDiagThick) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatTrsmTile(const Side &aSide, const UpperLower &aUpperLower, const Trans &aTrans,
                                                const Diag &aDiag, const T &aAlpha,
                                                void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
HicmaImplementation<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apSequence,
                                                       void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
HicmaImplementation<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                       void *apDescD, void *apSequence,
                                                       void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                          void *apDescDet) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const UpperLower &aUpperLower) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA, const UpperLower &aUpperLower) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatLaSetTile(const UpperLower &aUpperLower, T alpha, T beta,
                                                 void *apDescriptor) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatGeaddTile(const Trans &aTrans, const T &aAlpha, void *apDescA, const T &aBeta,
                                                 void *apDescB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatTrmmTile(const Side &aSide, const UpperLower &aUpperLower, const Trans &aTrans,
                                                const Diag &aDiag, const T &alpha,
                                                void *apDescA, void *apDescB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void
HicmaImplementation<T>::ExaGeoStatPosvTile(const UpperLower &aUpperLower, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void *HicmaImplementation<T>::ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) {
    return HICMA_RUNTIME_data_getaddr((HICMA_desc_t *) apA, aAm, aAn);
}