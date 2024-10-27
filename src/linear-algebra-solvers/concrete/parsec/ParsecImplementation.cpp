
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecImplementation.cpp
 * @brief This file contains the implementation of ParsecImplementation class.
 * @details ParsecImplementation is a concrete implementation of the LinearAlgebraMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-10-15
**/

#include <linear-algebra-solvers/concrete/parsec/ParsecImplementation.hpp>

using namespace exageostat::linearAlgebra;
using namespace exageostat::configurations;

template<typename T>
void ParsecImplementation<T>::ExaGeoStatSYRK(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

}

template<typename T>
void ParsecImplementation<T>::ExaGeoStatTLRCholesky(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

}

template<typename T>
void
ParsecImplementation<T>::ExaGeoStatNorm(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) {

}

template<typename T>
double
ParsecImplementation<T>::CalculateMSE(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) {

    return 0;
}

template<typename T>
T ParsecImplementation<T>::ModelingOperations(std::unique_ptr<ExaGeoStatData<T>> &aData,
                                              configurations::Configurations &aConfigurations, const double *apTheta,
                                              T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

}

template<typename T>
void ParsecImplementation<T>::ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) {

}


template<typename T>
void ParsecImplementation<T>::ExaGeoStatSequenceWait(void *apSequence) {

}


template<typename T>
void ParsecImplementation<T>::ExaGeoStatCreateSequence(void *apSequence) {

}


template<typename T>
void
ParsecImplementation<T>::ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apCD,
                                             void *apCrk, const int &aMaxRank, const int &aAcc) {

}


template<typename T>
void ParsecImplementation<T>::ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                                 const common::Trans &aTrans, const common::Diag &aDiag,
                                                 const T &aAlpha, void *apA, void *apCD, void *apCrk, void *apZ,
                                                 const int &aMaxRank) {

}


template<typename T>
void ParsecImplementation<T>::CopyDescriptors(void *apSourceDesc, void *apDestinationDesc, const int &aSize,
                                              const common::CopyDirection &aDirection) {

}