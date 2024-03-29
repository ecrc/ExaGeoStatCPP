
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecFunctions.cpp
 * @brief A class for static functions used in Parsec runtime
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-03-10
**/

#include <runtime/RuntimeFunctions.hpp>

using namespace std;
using namespace exageostat::common;
using namespace exageostat::runtime;
using namespace exageostat::dataunits;

//TODO: implement parsec functions
template<typename T>
void RuntimeFunctions<T>::CovarianceMatrix(DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                           const int &aTriangularPart, Locations<T> *apLocation1,
                                           Locations<T> *apLocation2, Locations<T> *apLocation3,
                                           T *apLocalTheta, const int &aDistanceMetric,
                                           const kernels::Kernel<T> *apKernel, void *apContext) {}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatMLETileAsyncMLOEMMOM(void *apDescExpr1, void *apDescExpr2, void *apDescExpr3,
                                                         void *apDescMLOE, void *apDescMMOM, void *apSequence,
                                                         void *apRequest, void *apContext) {}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatMLEMSPETileAsync(void *apDescZPredict, void *apDescZMiss, void *apDescError,
                                                     void *apSequence, void *apRequest, void *apContext) {}

template<typename T>
void RuntimeFunctions<T>::CopyDescriptorZ(DescriptorData<T> &aDescriptorData, void *apDescriptor, T *apDoubleVector,
                                          void *apContext) {}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatGaussianToNonTileAsync(DescriptorData<T> &aDescriptorData, void *apDesc, T *apTheta,
                                                      void *apContext) {}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatMeasureDetTileAsync(const Computation &aComputation, void *apDescA, void *apSequence,
                                                   void *apRequest, void *apDescDet, void *apContext) {}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB,
                                                         void *apDescC, void *apSequence, void *apRequest,
                                                         void *apContext) {}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apDescD,
                                                         void *apSequence, void *apRequest, void *apContext) {}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatMLETraceTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescNum,
                                                 void *apDescTrace, void *apContext) {}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatDoubleDotProduct(void *apDescA, void *apDescProduct, void *apSequence,
                                                void *apRequest,
                                                void *apContext) {}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatMLEMSPEBivariateTileAsync(void *apDescZPre, void *apDescZMiss, void *apDescsError1,
                                                         void *apDescsError2, void *apDescsError, void *apSequence,
                                                         void *apRequest, void *apContext) {}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatNonGaussianLogLikeTileAsync(const Computation &aComputation, void *apDescZ,
                                                                void *apDescSum, const T *apTheta, void *apSequence,
                                                                void *apRequest, void *apContext) {}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatNonGaussianTransformTileAsync(const Computation &aComputation, void *apDescZ,
                                                             const T *apTheta,
                                                             void *apSequence, void *apRequest, void *apContext) {}
