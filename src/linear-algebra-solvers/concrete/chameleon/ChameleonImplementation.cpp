
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
    aData->CalculateMedianLocations(kernel_name, median_locations);

    auto *CHAM_desc_C = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_C).chameleon_desc;
    auto *CHAM_desc_sub_C11 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_C11).chameleon_desc;
    auto *CHAM_desc_sub_C12 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_sub_C22 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_C22).chameleon_desc;
    auto *CHAM_desc_Z = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_desc_Z1 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_Z_1).chameleon_desc;
    auto *CHAM_desc_Z2 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_Z_2).chameleon_desc;
    auto *CHAM_desc_Z3 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_Z_3).chameleon_desc;
    auto *CHAM_desc_Zcpy = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_Z_COPY).chameleon_desc;
    auto *CHAM_desc_det = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_DETERMINANT).chameleon_desc;
    auto *CHAM_desc_product = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *CHAM_desc_sum = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_SUM).chameleon_desc;
    auto *CHAM_desc_product1 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_PRODUCT_1).chameleon_desc;
    auto *CHAM_desc_product2 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_PRODUCT_2).chameleon_desc;
    auto *CHAM_desc_product3 = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_PRODUCT_3).chameleon_desc;

    T *determinant = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT);
    *determinant = 0;
    T *product = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT);
    *product = 0;
    T *sum;
    if (aConfigurations.GetIsNonGaussian()) {
        sum = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, DESCRIPTOR_SUM);
        *sum = 0;
    }
    n = CHAM_desc_C->m;
    nhrs = CHAM_desc_Z->n;

    string recovery_file = aConfigurations.GetRecoveryFile();
    int iter_count = aData->GetMleIterations();

    if (recovery_file.empty() ||
        !(this->Recover((char *) (recovery_file.c_str()), iter_count, (T *) theta, &loglik, num_params))) {
        START_TIMING(dzcpy_time);
        if (iter_count == 0) {
            // Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            this->ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
        } else {
            VERBOSE("\tRe-store the original Z vector...")
            this->ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Zcpy, CHAM_desc_Z);
            VERBOSE("\tDone.")
        }
        STOP_TIMING(dzcpy_time);
    }

    //Generate new co-variance matrix C based on new theta
    VERBOSE("\tGenerate New Covariance Matrix...")
    START_TIMING(matrix_gen_time);

    if (kernel_name == "bivariate_matern_parsimonious2" ||
        kernel_name == "bivariate_matern_parsimonious2_profile") {

        int upper_lower = EXAGEOSTAT_UPPER_LOWER;
        univariate_theta = new T[3];
        univariate2_theta = new T[3];
        univariate3_theta = new T[3];
        univariate_theta[0] = theta[0];
        univariate_theta[1] = theta[2];
        univariate_theta[2] = theta[3];

        int distance_metric = aConfigurations.GetDistanceMetric();

        RuntimeFunctions<T>::CovarianceMatrix(*aData->GetDescriptorData(), CHAM_desc_sub_C11, upper_lower,
                                              aData->GetLocations(), aData->GetLocations(), &median_locations,
                                              univariate_theta, distance_metric, &aKernel);

        nu12 = 0.5 * (theta[3] + theta[4]);
        rho = theta[5] * sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) / (tgamma(theta[3]) * tgamma(theta[4]))) *
              tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(theta[0] * theta[1]);
        univariate2_theta[0] = sigma_square12;
        univariate2_theta[1] = theta[2];
        univariate2_theta[2] = nu12;

        RuntimeFunctions<T>::CovarianceMatrix(*aData->GetDescriptorData(), CHAM_desc_sub_C12, upper_lower,
                                              &median_locations, aData->GetLocations(), &median_locations,
                                              univariate2_theta, 0, &aKernel);
        STOP_TIMING(matrix_gen_time);

        univariate3_theta[0] = theta[1];
        univariate3_theta[1] = theta[2];
        univariate3_theta[2] = theta[4];

        RuntimeFunctions<T>::CovarianceMatrix(*aData->GetDescriptorData(), CHAM_desc_sub_C22, upper_lower,
                                              &median_locations, aData->GetLocations(), &median_locations,
                                              univariate2_theta, 0, &aKernel);
    } else {
        int upper_lower = EXAGEOSTAT_LOWER;
        RuntimeFunctions<T>::CovarianceMatrix(*aData->GetDescriptorData(), CHAM_desc_C, upper_lower,
                                              aData->GetLocations(), aData->GetLocations(), &median_locations,
                                              (T *) theta, 0, &aKernel);
    }
    this->ExaGeoStatSequenceWait(pSequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE("\tDone.")

    VERBOSE("\tCholesky factorization of Sigma...")
    START_TIMING(time_facto);
    this->ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_desc_C, aConfigurations.GetBand(), nullptr, nullptr, 0, 0);
    STOP_TIMING(time_facto);
    flops += flops_dpotrf(n);
    VERBOSE("\tDone.")

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("\tCalculating the log determinant ...")
    START_TIMING(logdet_calculate);
    RuntimeFunctions<T>::ExaGeoStatMeasureDetTileAsync(aConfigurations.GetComputation(), CHAM_desc_C, pSequence,
                                                       &request_array, CHAM_desc_det);
    this->ExaGeoStatSequenceWait(pSequence);

    logdet = 2 * (*determinant);
    STOP_TIMING(logdet_calculate);
    VERBOSE("\tDone.")

    if (aConfigurations.GetIsNonGaussian()) {
        VERBOSE("Transform Z vector to Gaussian field ...")
        RuntimeFunctions<T>::ExaGeoStatNonGaussianTransformTileAsync(aConfigurations.GetComputation(), CHAM_desc_Z,
                                                                     (T *) theta, pSequence, &request_array[0]);
        this->ExaGeoStatSequenceWait(pSequence);
        VERBOSE("\tDone.")

        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z, CHAM_desc_Z, 0, CHAM_desc_product);

        VERBOSE("Calculate non-Gaussian loglik ...")
        RuntimeFunctions<T>::ExaGeoStatNonGaussianLogLikeTileAsync(aConfigurations.GetComputation(), CHAM_desc_Z,
                                                                   CHAM_desc_sum, (T *) theta, pSequence,
                                                                   &request_array[0]);
        this->ExaGeoStatSequenceWait(pSequence);
        VERBOSE("\tDone.")
    }

    // Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("\tSolving the linear system ...")
    START_TIMING(time_solve);
    this->ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                             CHAM_desc_C, nullptr, nullptr, CHAM_desc_Z, 0);
    STOP_TIMING(time_solve);
    flops += flops_dtrsm(ChamLeft, n, nhrs);
    VERBOSE("\tDone.")

    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...")
    RuntimeFunctions<T>::ExaGeoStatDoubleDotProduct(CHAM_desc_Z, CHAM_desc_product,
                                                    pSequence, request_array);
    ExaGeoStatSequenceWait(pSequence);

    if (kernel_name == "BivariateMaternParsimonious2Profile") {
        loglik =
                -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet - (T) (n / 2.0) * log(2.0 * PI);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0, CHAM_desc_product1);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0, CHAM_desc_product2);
        variance1 = (1.0 / (n / 2)) * dot_product1;
        variance2 = (1.0 / (n / 2)) * dot_product2;

    } else if (kernel_name == "BivariateMaternParsimoniousProfile") {
        loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet -
                 (double) (n / 2.0) * log(2.0 * PI);
        RuntimeFunctions<T>::ExaGeoStaStrideVectorTileAsync(CHAM_desc_Z, CHAM_desc_Z1, CHAM_desc_Z2, pSequence,
                                                            &request_array[0]);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0, CHAM_desc_product1);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0, CHAM_desc_product2);
        variance1 = (1.0 / (n / 2)) * dot_product1;
        variance2 = (1.0 / (n / 2)) * dot_product2;
    } else if (kernel_name == "TrivariateMaternParsimoniousProfile") {

        loglik = -(n / 3.0) + (n / 3.0) * log(n / 3.0) - (n / 3.0) * log(dot_product) - 0.5 * logdet -
                 (double) (n / 3.0) * log(2.0 * PI);
        //to be optimized
        RuntimeFunctions<T>::ExaGeoStaStrideVectorTileAsync(CHAM_desc_Z, CHAM_desc_Z1, CHAM_desc_Z2, CHAM_desc_Z3,
                                                            pSequence, &request_array[0]);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0, CHAM_desc_product1);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0, CHAM_desc_product2);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z3, CHAM_desc_Z3, 0, CHAM_desc_product3);
        variance1 = (1.0 / (n / 3.0)) * dot_product1;
        variance2 = (1.0 / (n / 3.0)) * dot_product2;
    } else {
        dot_product = *product;
        loglik = -0.5 * dot_product - 0.5 * logdet;
        if (aConfigurations.GetIsNonGaussian()) {
            loglik = loglik - *sum - n * log(theta[3]) - (double) (n / 2.0) * log(2.0 * PI);
        } else {
            loglik = loglik - (double) (n / 2.0) * log(2.0 * PI);
        }
    }
    VERBOSE("\tDone.")

    //Distribute the values in the case of MPI
#ifdef USE_MPI
    MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    LOGGER("\t" << iter_count + 1 << " - Model Parameters (", true)

    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), "\t %d- Model Parameters (", iter_count + 1);
    }

    if ((aConfigurations.GetKernelName() == "BivariateMaternParsimoniousProfile") ||
        (aConfigurations.GetKernelName() == "BivariateMaternParsimonious2Profile")) {
        LOGGER(setprecision(8) << variance1 << setprecision(8) << variance2)
        if (aConfigurations.GetLogger()) {
            fprintf(aConfigurations.GetFileLogPath(), "%.8f, %.8f,", variance1, variance2);
        }
        i = 2;
    } else {
        i = 0;
    }
    for (; i < num_params; i++) {
        LOGGER_PRECISION(theta[i])
        if (i < num_params - 1) {
            LOGGER_PRECISION(", ")
        }
        if (aConfigurations.GetLogger()) {
            fprintf(aConfigurations.GetFileLogPath(), "%.8f, ", theta[i]);
        }
    }

    LOGGER_PRECISION(")----> LogLi: " << loglik << "\n", 18)

    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), ")----> LogLi: %.18f\n", loglik);
    }
    VERBOSE("---- Facto Time: " << time_facto)
    VERBOSE("---- Log Determent Time: " << logdet_calculate)
    VERBOSE("---- dtrsm Time: " << time_solve)
    VERBOSE("---- Matrix Generation Time: " << matrix_gen_time)
    VERBOSE("---- Total Time: " << time_facto + logdet_calculate + time_solve)
    VERBOSE("---- Gflop/s: " << flops / 1e9 / (time_facto + time_solve))

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

    return loglik;
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
