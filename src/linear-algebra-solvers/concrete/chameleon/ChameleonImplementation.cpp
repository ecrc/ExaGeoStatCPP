
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.cpp
 * @brief This file contains the declaration of ChameleonImplementation class.
 * @details ChameleonImplementation is a concrete implementation of LinearAlgebraMethods class for dense or diagonal-super tile matrices..
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifdef USE_MPI

#include <mpi.h>

#endif

#include <cblas.h>

#include <linear-algebra-solvers/concrete/chameleon/ChameleonImplementation.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::runtime;
using namespace exageostat::results;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;

template<typename T>
T ChameleonImplementation<T>::ExaGeoStatMLETile(std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                configurations::Configurations &aConfigurations, const double *theta,
                                                T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

    double value;
    if (!aData->GetDescriptorData()->GetIsDescriptorInitiated()) {
        this->InitiateDescriptors(aConfigurations, *aData->GetDescriptorData(), aKernel.GetVariablesNumber(),
                                  apMeasurementsMatrix);
    }
    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    if (!aData->GetDescriptorData()->GetSequence()) {
        RUNTIME_sequence_t *sequence;
        this->ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    }

    auto pSequence = (RUNTIME_sequence_t *) aData->GetDescriptorData()->GetSequence();
    //Initialization
    T loglik = 0.0, logdet, variance, variance1 = 1, variance2 = 1, variance3, dot_product = 0, dot_product1 = 0, dot_product2 = 0, dot_product3, n, dzcpy_time, time_facto, time_solve, logdet_calculate, matrix_gen_time;
    double accumulated_executed_time, accumulated_flops;

    int nhrs, i;
    T flops = 0.0;
    T *univariate_theta, *univariate2_theta, *univariate3_theta, nu12, rho, sigma_square12;

    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = aKernel.GetParametersNumbers();
    auto median_locations = Locations<T>(1, aData->GetLocations()->GetDimension());

    auto *CHAM_desc_Zobs = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_desc_X = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_X).chameleon_desc;
    auto *CHAM_desc_part1 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_PART1).chameleon_desc;
    auto *CHAM_desc_part2 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_PART2).chameleon_desc;
    auto *CHAM_desc_part2_vector = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_PART2_VECTOR).chameleon_desc;
    auto *CHAM_desc_XtX = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_XtX).chameleon_desc;


    T *part1 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PART1);
    T *part2 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PART2);
    *part2 = 0;
    n = CHAM_desc_X->m;
    int iter_count = aData->GetMleIterations();
    int no_years=751;
    double* localtheta =  (double *) malloc((3+no_years) * sizeof(double));
    int M = 10;
    int t = 8760;
    localtheta[0]= theta[0];
    localtheta[1]= t;
    localtheta[2]= M;

    int upper_lower = EXAGEOSTAT_UPPER_LOWER;
    int distance_metric = aConfigurations.GetDistanceMetric();
    for(int ii=0;ii<no_years;ii++) {
        localtheta[3+ii] = apMeasurementsMatrix[ii];
    }
//    fprintf(stderr, "%f %f %f %f\n", localtheta[0], localtheta[1], localtheta[2], localtheta[3]);

    RuntimeFunctions<T>::CovarianceMatrix(*aData->GetDescriptorData(), CHAM_desc_X, upper_lower,
                                              aData->GetLocations(), aData->GetLocations(), &median_locations,
                                          (T*)localtheta, distance_metric, &aKernel);
    this->ExaGeoStatSequenceWait(pSequence);
    //calculate part1
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, CHAM_desc_Zobs, CHAM_desc_Zobs, 0.0, CHAM_desc_part1);

    //Calculate part2
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, CHAM_desc_X, CHAM_desc_Zobs, 0, CHAM_desc_part2_vector);

    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, CHAM_desc_X, CHAM_desc_X, 0, CHAM_desc_XtX);

    this->ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_desc_XtX, aConfigurations.GetBand(), nullptr, nullptr, 0, 0);

    this->ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                             CHAM_desc_XtX, nullptr, nullptr, CHAM_desc_part2_vector, 0);

    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, CHAM_desc_part2_vector, CHAM_desc_part2_vector, 0, CHAM_desc_part2);

    value =(-1.0) * log(*part1-*part2);

//    LOGGER("\t" << iter_count + 1 << " - Model Parameters (", true)
//    for (i=0; i < num_params; i++) {
//        LOGGER_PRECISION(theta[i])
//        if (i < num_params - 1) {
//            LOGGER_PRECISION(", ")
//        }
//        if (aConfigurations.GetLogger()) {
//            fprintf(aConfigurations.GetFileLogPath(), "%.8f, ", theta[i]);
//        }
//    }
//
//    LOGGER_PRECISION(")----> LogLi: " << loglik << "\n", 18)
//
//    if (aConfigurations.GetLogger()) {
//        fprintf(aConfigurations.GetFileLogPath(), ")----> LogLi: %.18f\n", loglik);
//    }
//    VERBOSE("---- Facto Time: " << time_facto)
//    VERBOSE("---- Log Determent Time: " << logdet_calculate)
//    VERBOSE("---- dtrsm Time: " << time_solve)
//    VERBOSE("---- Matrix Generation Time: " << matrix_gen_time)
//    VERBOSE("---- Total Time: " << time_facto + logdet_calculate + time_solve)
//    VERBOSE("---- Gflop/s: " << flops / 1e9 / (time_facto + time_solve))

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&value,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast(theta,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	if ( MORSE_My_Mpi_Rank() == 0 )
	{
#endif

    fprintf(stderr, "%d- data->part1:%f, data->part2: %f,  theta[0]:%0.14f ------- lh: %0.18f\n",   iter_count + 1, *part1, *part2, theta[0], value);
#if defined(CHAMELEON_USE_MPI)
    }
#endif
    aData->SetMleIterations(aData->GetMleIterations() + 1);

    // for experiments and benchmarking
    accumulated_executed_time =
            Results::GetInstance()->GetTotalModelingExecutionTime() + time_facto + logdet_calculate +
            time_solve;
    Results::GetInstance()->SetTotalModelingExecutionTime(accumulated_executed_time);
    accumulated_flops =
            Results::GetInstance()->GetTotalModelingFlops() + (flops / 1e9 / (time_facto + time_solve));
    Results::GetInstance()->SetTotalModelingFlops(accumulated_flops);

    Results::GetInstance()->SetMLEIterations(iter_count + 1);
    Results::GetInstance()->SetMaximumTheta(vector<double>(theta, theta + num_params));
    Results::GetInstance()->SetLogLikValue(loglik);

    return value;
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatLapackCopyTile(const UpperLower &aUpperLower, void *apA, void *apB) {
    int status = CHAMELEON_dlacpy_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dlacpy_Tile Failed!");
    }
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatTrsmTile(const Side &aSide, const UpperLower &aUpperLower,
                                                    const Trans &aTrans, const Diag &aDiag,
                                                    const T &aAlpha, void *apA, void *apCD, void *apCrk, void *apZ,
                                                    const int &aMaxRank) {
    int status = CHAMELEON_dtrsm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans,
                                      (cham_diag_t) aDiag, aAlpha, (CHAM_desc_t *) apA, (CHAM_desc_t *) apZ);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dtrsm_Tile Failed!");
    }
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatSequenceWait(void *apSequence) {
    int status = CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apSequence);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Sequence_Wait Failed!");
    }
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatCreateSequence(void *apSequence) {
    int status = CHAMELEON_Sequence_Create(((RUNTIME_sequence_t **) apSequence));
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Sequence_Create Failed!");
    }
}
