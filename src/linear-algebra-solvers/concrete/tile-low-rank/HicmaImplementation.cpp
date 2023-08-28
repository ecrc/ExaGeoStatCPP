
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.cpp
 * @brief Sets up the HiCMA descriptors needed for the tile low rank computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-03-26
**/

#include <lapacke.h>

extern "C" {
#include <hicma.h>
#include <control/hicma_context.h>
}

#include <linear-algebra-solvers/concrete/tile-low-rank/HicmaImplementation.hpp>
#include <data-units/DescriptorData.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

template<typename T>
void
HicmaImplementation<T>::InitiateDescriptors(Configurations &apConfigurations, DescriptorData <T> &aDescriptorData,
                                            T *apMeasurementsMatrix) {

    // Check for Initialise the Hicma context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    // Create a Hicma sequence
    HICMA_sequence_t *pSequence;
    HICMA_request_t request[2] = {HICMA_SUCCESS, HICMA_SUCCESS};
    HICMA_Sequence_Create(&pSequence);

    int N = apConfigurations.GetProblemSize();
    int lts = apConfigurations.GetLowTileSize();
    int p_grid = apConfigurations.GetPGrid();
    int q_grid = apConfigurations.GetQGrid();
    bool is_OOC = apConfigurations.GetIsOOC();
    int max_rank = apConfigurations.GetMaxRank();
    int nZmiss = apConfigurations.GetUnknownObservationsNb();
    T mean_square_error = apConfigurations.GetMeanSquareError();
    int approximation_mode = apConfigurations.GetApproximationMode();
    string actual_observations_path = apConfigurations.GetActualObservationsFilePath();

    int nZobs_value;
    if (actual_observations_path.empty()) {
        nZobs_value = N - nZmiss;
    } else {
        nZobs_value = N;
    }

    apConfigurations.SetKnownObservationsValues(nZobs_value);
    int nZobs = apConfigurations.GetKnownObservationsValues();

    // For distributed system and should be removed
    T *Zcpy = new T[N];

    int MBC, NBC, MC, NC;
    int MBD, NBD, MD, ND;
    int MBUV, NBUV, MUV, NUV;
    int MBrk, NBrk, Mrk, Nrk;

    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    }

    //CDense Descriptor
    if (approximation_mode == 1) {
        MBC = lts;
        NBC = lts;
        MC = N;
        NC = N;
    } else {
        MBC = 1;
        NBC = 1;
        MC = lts;
        NC = lts;
    }

    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C, is_OOC, nullptr, float_point, MBC, NBC,
                                  MBC * NBC, MC, NC, 0, 0, MC, NC, p_grid, q_grid);
    //CAD Descriptor
    MBD = lts;
    NBD = lts;
    MD = N;
    ND = MBD;
    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CD, is_OOC, nullptr, float_point, MBD, NBD,
                                  MBD * NBD,
                                  MD, ND, 0, 0, MD, ND, p_grid, q_grid);

    //CUV Descriptor
    MBUV = lts;
    NBUV = 2 * max_rank;
    int N_over_lts_times_lts = N / lts * lts;
    if (N_over_lts_times_lts < N) {
        MUV = N_over_lts_times_lts + lts;
    } else if (N_over_lts_times_lts == N) {
        MUV = N_over_lts_times_lts;
    } else {
        throw range_error("Invalid value. This case should not happen, Please make sure of N and lts values.");
    }

    T expr = (T) MUV / (T) lts;
    NUV = 2 * expr * max_rank;
    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CUV, is_OOC, nullptr, float_point, MBUV,
                                  NBUV,
                                  MBUV * NBUV, MUV, NUV, 0, 0, MUV, NUV, p_grid, q_grid);

    //Crk Descriptor
    MBrk = 1;
    NBrk = 1;
    auto *desc_cuv = aDescriptorData.GetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CUV).hicma_desc;
    Mrk = desc_cuv->mt;
    Nrk = desc_cuv->mt;
    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CRK, is_OOC, nullptr, float_point, MBrk,
                                  NBrk,
                                  MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, nullptr, float_point, lts, lts,
                                  lts * lts, N, 1, 0, 0, N, 1, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, nullptr, float_point, lts,
                                  lts,
                                  lts * lts,
                                  N, 1, 0, 0, N, 1, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_DETERMINANT, is_OOC, nullptr,
                                  float_point, lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    if (nZmiss != 0) {
        if (actual_observations_path.empty()) {
            aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z_OBSERVATIONS, is_OOC,
                                          &Zcpy[nZmiss],
                                          float_point, lts, lts, lts * lts, nZobs, 1, 0, 0, nZobs,
                                          1, p_grid, q_grid);
        } else {

            aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z_OBSERVATIONS, is_OOC, nullptr,
                                          float_point,
                                          lts, lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);

        }
        //C12AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZmiss;
        ND = MBD;
        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z_Actual, is_OOC, nullptr, float_point,
                                      lts,
                                      lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);


        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C12D, is_OOC, nullptr, float_point, MBD,
                                      NBD,
                                      MBD * NBD, MD, ND, 0, 0, MD, ND, p_grid, q_grid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * max_rank;

        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C12UV, is_OOC, nullptr, float_point,
                                      MBUV,
                                      NBUV,
                                      MBUV * NBUV, MBUV, NBUV, 0, 0, MBUV, NBUV, p_grid, q_grid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        auto *desc_c12uv = aDescriptorData.GetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C12UV).hicma_desc;
        Mrk = desc_c12uv->mt;
        Nrk = desc_c12uv->mt;
        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C12RK, is_OOC, nullptr, float_point,
                                      MBrk,
                                      NBrk,
                                      MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, p_grid, q_grid);

        //C22D Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZobs;
        ND = MBD;

        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C22D, is_OOC, nullptr, float_point, MBD,
                                      NBD,
                                      MBD * NBD, MD, ND, 0, 0, MD, ND, p_grid, q_grid);

        //C22UV Descriptor
        MBUV = lts;
        NBUV = 2 * max_rank;
        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C22UV, is_OOC, nullptr, float_point,
                                      MBUV,
                                      NBUV,
                                      MBUV * NBUV, MBUV, NBUV, 0, 0, MBUV, NBUV, p_grid, q_grid);

        //C22Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        auto *desc_c22uv = aDescriptorData.GetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C22UV).hicma_desc;
        Mrk = desc_c22uv->mt;
        Nrk = desc_c22uv->mt;

        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_C22RK, is_OOC, nullptr, float_point,
                                      MBrk,
                                      NBrk,
                                      MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, p_grid, q_grid);

        //Other descriptors
        aDescriptorData.SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_MSE, is_OOC, nullptr, float_point, lts,
                                      lts, lts * lts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    }

    aDescriptorData.SetSequence(pSequence);
    aDescriptorData.SetRequest(request);

    //stop gsl error handler
    gsl_set_error_handler_off();
    delete[] Zcpy;
    aDescriptorData.SetIsDescriptorInitiated(true);
}

template<typename T>
void
HicmaImplementation<T>::ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData <T> *apDescriptorData,
                                                         void *apDesc, T *apTheta) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void
HicmaImplementation<T>::CovarianceMatrixCodelet(DescriptorData <T> *apDescriptorData, void *apDescriptor,
                                                int &aTriangularPart,
                                                dataunits::Locations <T> *apLocation1,
                                                dataunits::Locations <T> *apLocation2,
                                                dataunits::Locations <T> *apLocation3,
                                                T *aLocalTheta, int aDistanceMetric,
                                                kernels::Kernel <T> *apKernel) {

    // Check for Initialise the Hicma context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    HICMA_option_t options;
    HICMA_RUNTIME_options_init(&options, (HICMA_context_t * )
    this->mpContext,
            (HICMA_sequence_t *) apDescriptorData->GetSequence(),
            (HICMA_request_t *) apDescriptorData->GetRequest());
    int tempmm, tempnn;

    auto *HICMA_apDescriptor = (HICMA_desc_t *) apDescriptor;
    HICMA_desc_t A = *HICMA_apDescriptor;
    struct starpu_codelet *cl = &this->cl_dcmg;
    int m, n, m0 = 0, n0 = 0;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (aTriangularPart == HicmaUpperLower) {
            m = 0;
        } else {
            m = A.m == A.n ? n : 0;
        }
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;

            // Register the data with StarPU
            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) HICMA_RUNTIME_data_getaddr(HICMA_apDescriptor, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations < T > *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations < T > *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations < T > *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(exageostat::kernels::Kernel < T > *),
                               0);
        }
    }
    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, (HICMA_context_t * )
    this->mpContext);

    HICMA_Sequence_Wait((HICMA_sequence_t *) apDescriptorData->GetSequence());

}

template<typename T>
void HicmaImplementation<T>::GenerateObservationsVector(Configurations &apConfigurations,
                                                        DescriptorData <T> *apDescriptorData,
                                                        BaseDescriptor aDescriptor,
                                                        Locations <T> *apLocation1, Locations <T> *apLocation2,
                                                        Locations <T> *apLocation3, int aDistanceMetric) {

    // Check for Initialise the Hicma context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }
    int N = apConfigurations.GetProblemSize();
    int seed = apConfigurations.GetSeed();
    int iseed[4] = {seed, seed, seed, 1};
    auto *pDescriptor = aDescriptor.hicma_desc;

    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *Nrand = new T[N];
    LAPACKE_dlarnv(3, iseed, N, (double *) Nrand);


    //Generate the co-variance matrix C
    auto *theta = new T[apConfigurations.GetInitialTheta().size()];
    for (int i = 0; i < apConfigurations.GetInitialTheta().size(); i++) {
        theta[i] = apConfigurations.GetInitialTheta()[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    // Register and create a kernel object
    Kernel <T> *kernel = exageostat::plugins::PluginRegistry < Kernel < T >> ::Create(apConfigurations.GetKernelName());

    this->CovarianceMatrixCodelet(apDescriptorData, pDescriptor, upper_lower, apLocation1, apLocation2, apLocation3,
                                  theta,
                                  aDistanceMetric, kernel);
    delete[] theta;
    VERBOSE("Done.")

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto *HICMA_descriptorZ = apDescriptorData->GetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z).hicma_desc;
    CopyDescriptorZ(apDescriptorData, HICMA_descriptorZ, Nrand);
    VERBOSE("Done.")
    delete[] Nrand;
    //// RESET OF THE IMPLEMENTATION WILL BE ADDED AFTER FINALIZING ALL MODULES WITH EXACT.
}

template<typename T>
void
HicmaImplementation<T>::CopyDescriptorZ(DescriptorData <T> *apDescriptorData, void *apDescriptor, T *apDoubleVector) {
    throw std::runtime_error("unimplemented for now");
}


template<typename T>
T HicmaImplementation<T>::ExaGeoStatMleTile(const hardware::ExaGeoStatHardware &apHardware,
                                            dataunits::ExaGeoStatData <T> &apData,
                                            configurations::Configurations &apConfigurations, const double *theta,
                                            T *apMeasurementsMatrix) {

    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatLapackCopyTile(exageostat::common::UpperLower aUpperLower, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");

}

template<typename T>
int
HicmaImplementation<T>::ExaGeoStatLapackToDescriptor(exageostat::common::UpperLower aUpperLower, void *apAf77, int aLda,
                                                     void *apA) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatSequenceWait(void *apSequence) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatPotrfTile(exageostat::common::UpperLower aUpperLower, void *apA) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans,
                                               common::Diag aDiag, T aAlpha, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha,
                                               void *apA, void *apB, T aBeta, void *apC) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                           void *apSequence, void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                          void *apDescDet) {
    throw std::runtime_error("unimplemented for now");
}
