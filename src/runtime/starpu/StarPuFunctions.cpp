
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarPuFunctions.cpp
 * @brief A class for static functions that make use of starpu_codelets.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-28
**/

#include <data-units/ExaGeoStatData.hpp>
#include <runtime/RuntimeFunctions.hpp>
#include <runtime/starpu/StarPuCodeletsHeaders.hpp>
#include <runtime/starpu/helpers/StarPuHelpersFactory.hpp>

using namespace exageostat::common;
using namespace exageostat::runtime;
using namespace exageostat::dataunits;

template<typename T>
void RuntimeFunctions<T>::CovarianceMatrix(dataunits::DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                           const int &aTriangularPart,
                                           dataunits::Locations<T> *apLocation1, dataunits::Locations<T> *apLocation2,
                                           dataunits::Locations<T> *apLocation3, T *apLocalTheta,
                                           const int &aDistanceMetric,
                                           const kernels::Kernel<T> *apKernel) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, aDescriptorData.GetSequence(),
                                         aDescriptorData.GetRequest());

    DCMGCodelet<T> cl;
    cl.InsertTask(apDescriptor, aTriangularPart, apLocation1, apLocation2, apLocation3, apLocalTheta, aDistanceMetric,
                  apKernel);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatMLETileAsyncMLOEMMOM(void *apDescExpr1, void *apDescExpr2, void *apDescExpr3,
                                                         void *apDescMLOE,
                                                         void *apDescMMOM, void *apSequence, void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    DmloeMmomCodelet<T> cl;
    cl.InsertTask(apDescExpr1, apDescExpr2, apDescExpr3, apDescMLOE, apDescMMOM);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatMLEMSPETileAsync(void *apDescZPredict, void *apDescZMiss, void *apDescError,
                                                     void *apSequence,
                                                     void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    DMSECodelet<T> cl;
    cl.InsertTask(apDescError, apDescZPredict, apDescZMiss);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::CopyDescriptorZ(dataunits::DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                          T *apDoubleVector) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, aDescriptorData.GetSequence(),
                                         aDescriptorData.GetRequest());

    DZCPYCodelet<T> cl;
    cl.InsertTask(apDescriptor, apDoubleVector);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData<T> &aDescriptorData, void *apDesc,
                                                           T *apTheta) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, aDescriptorData.GetSequence(),
                                         aDescriptorData.GetRequest());

    GaussianCodelet<T> cl;
    cl.InsertTask(apDesc, apTheta);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatMeasureDetTileAsync(const common::Computation &aComputation, void *apDescA,
                                                   void *apSequence, void *apRequest,
                                                   void *apDescDet) {
    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(aComputation);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    DMDETCodelet<T> cl;
    cl.InsertTask(aComputation, apDescA, apDescDet, starpu_helper);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apSequence,
                                                         void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    STRIDEVECCodelet<T> cl;
    cl.InsertTask(apDescA, apDescB, apDescC);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apDescD,
                                                         void *apSequence,
                                                         void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    TriStrideVecCodelet<T> cl;
    cl.InsertTask(apDescA, apDescB, apDescC, apDescD);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatMLETraceTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescNum,
                                                      void *apDescTrace) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    DTRACECodelet<T> cl;
    cl.InsertTask(apDescA, apDescNum, apDescTrace);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatDoubleDotProduct(void *apDescA, void *apDescProduct, void *apSequence, void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    DDOTPCodelet<T> cl;
    cl.InsertTask(apDescA, apDescProduct);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatMLEMSPEBivariateTileAsync(void *apDescZPre, void *apDescZMiss, void *apDescError1,
                                                              void *apDescError2,
                                                              void *apDescError, void *apSequence, void *apRequest) {
    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(EXACT_DENSE);
    auto *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    DMSEBivariateCodelet<T> cl;
    cl.InsertTask(apDescZPre, apDescError, apDescError1, apDescError2, apDescZMiss);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void RuntimeFunctions<T>::ExaGeoStatNonGaussianLogLikeTileAsync(const common::Computation &aComputation, void *apDescZ,
                                                                void *apDescSum,
                                                                const T *apTheta, void *apSequence, void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(aComputation);
    void *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    NonGaussianLoglike<T> cl;
    cl.InsertTask(apDescZ, apDescSum, apTheta, starpu_helper);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}

template<typename T>
void
RuntimeFunctions<T>::ExaGeoStatNonGaussianTransformTileAsync(const common::Computation &aComputation, void *apDescZ,
                                                             const T *apTheta,
                                                             void *apSequence, void *apRequest) {

    auto starpu_helper = StarPuHelpersFactory::CreateStarPuHelper(aComputation);
    void *pOptions = starpu_helper->GetOptions();
    starpu_helper->ExaGeoStatOptionsInit(pOptions, apSequence, apRequest);

    NonGaussianTransform<T> cl;
    cl.InsertTask(apDescZ, apTheta, starpu_helper);

    starpu_helper->ExaGeoStatOptionsFree(pOptions);
    starpu_helper->ExaGeoStatOptionsFinalize(pOptions);
    starpu_helper->DeleteOptions(pOptions);

}
