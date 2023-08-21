
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDense.cpp
 * @brief Dense Tile implementation of linear algebra methods.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <lapacke.h>

// Include Chameleon libraries
extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/descriptor.h>
#include <control/context.h>
}

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>
#include <data-units/DescriptorData.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors(Configurations &aConfigurations,
                                                          DescriptorData<T> &aDescriptorData) {

    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    // Get the problem size and other configuration parameters
    int N = aConfigurations.GetProblemSize();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));

    // Create a Chameleon sequence
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request_array[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&pSequence);


    // Set the floating point precision based on the template type
    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else if (sizeof(T) == SIZE_OF_DOUBLE) {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    } else {
        throw runtime_error("Unsupported for now!");
    }

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C, is_OOC, nullptr, float_point, dts, dts,
                                  dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, nullptr, float_point, dts, dts,
                                  dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, nullptr, float_point, dts,
                                  dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT, is_OOC, nullptr, float_point,
                                  dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT, is_OOC, nullptr,
                                  float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    if (float_point == EXAGEOSTAT_REAL_DOUBLE) {
        auto *CHAM_descC = aDescriptorData.GetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C11, is_OOC, nullptr, float_point, dts,
                                      dts,
                                      dts * dts, N, N, 0, 0, CHAM_descC->m / 2, CHAM_descC->n / 2, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr, float_point, dts,
                                      dts,
                                      dts * dts, N, N, CHAM_descC->m / 2, 0, CHAM_descC->m / 2, CHAM_descC->n / 2,
                                      p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, is_OOC, nullptr, float_point, dts,
                                      dts,
                                      dts * dts, N, N, CHAM_descC->m / 2, CHAM_descC->n / 2, CHAM_descC->m / 2,
                                      CHAM_descC->n / 2, p_grid, q_grid);

        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_1, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_2, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1, is_OOC, nullptr,
                                      float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2, is_OOC, nullptr,
                                      float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }

    aDescriptorData.SetSequence(pSequence);
    aDescriptorData.SetRequest(request_array);

    //stop gsl error handler
    gsl_set_error_handler_off();
    free(Zcpy);
}

template<typename T>
void ChameleonImplementationDense<T>::GenerateObservationsVector(Configurations &aConfigurations,
                                                                 DescriptorData<T> *apDescriptorData,
                                                                 BaseDescriptor aDescriptor, Locations<T> *apLocation1,
                                                                 Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                                 int aDistanceMetric) {

    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    const int N = aConfigurations.GetProblemSize();
    int seed = aConfigurations.GetSeed();
    int iseed[4] = {seed, seed, seed, 1};
    auto *p_descriptor = aDescriptor.chameleon_desc;

    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *Nrand = (T *) malloc(N * sizeof(T));
    LAPACKE_dlarnv(3, iseed, N, (double *) Nrand);

    //Generate the co-variance matrix C
    auto *theta = (T *) malloc(aConfigurations.GetInitialTheta().size() * sizeof(T));
    for (int i = 0; i < aConfigurations.GetInitialTheta().size(); i++) {
        theta[i] = aConfigurations.GetInitialTheta()[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    // Register and create a kernel object
    Kernel<T> *kernel = exageostat::plugins::PluginRegistry<Kernel<T>>::Create(aConfigurations.GetKernelName());

    this->CovarianceMatrixCodelet(apDescriptorData, p_descriptor, upper_lower, apLocation1, apLocation2, apLocation3,
                                  theta,
                                  aDistanceMetric, kernel);

    delete kernel;
    VERBOSE("Done.\n")

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto *CHAM_descriptorZ = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
    CopyDescriptorZ(apDescriptorData, CHAM_descriptorZ, Nrand);
    VERBOSE("Done.\n")

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")
    int potential_failure = CHAMELEON_dpotrf_Tile(ChamLower, p_descriptor);
    FAILURE_LOGGER(potential_failure, "Factorization cannot be performed..\nThe matrix is not positive definite")
    VERBOSE("Done.\n")

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, p_descriptor, CHAM_descriptorZ);
    VERBOSE("Done.\n")


    if (aConfigurations.GetKernelName() == "UnivariateMaternNonGaussian") {
        //Gaussian to non-gaussian transformation
        VERBOSE("Convert Z Gaussian to non-Gaussian (Synthetic Dataset Generation Phase) .....")
        ExaGeoStatGaussianToNonTileAsync(apDescriptorData, CHAM_descriptorZ, theta);
        VERBOSE(" Done.\n")
    }
    free(theta);
    const int P = aConfigurations.GetP();
    if (aConfigurations.GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = (T*) malloc(N * sizeof(T));
        CHAMELEON_Tile_to_Lapack( CHAM_descriptorZ, pMatrix, N);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, aConfigurations->GetLoggerPath(), apLocation1);
        }
        free(pMatrix);
#else
        pMatrix = (T *) CHAM_descriptorZ->mat;
        string path = aConfigurations.GetLoggerPath();
        DiskWriter<T>::WriteVectorsToDisk(*pMatrix, N, P, path, *apLocation1);
#endif
        VERBOSE(" Done.\n")
    }

    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, p_descriptor);
    free(Nrand);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
T ChameleonImplementationDense<T>::ExaGeoStatMleTile(hardware::ExaGeoStatHardware &aHardware, ExaGeoStatData<T> *apData,
                                                     configurations::Configurations *apConfigurations,
                                                     const double *theta) {

    this->SetContext(aHardware.GetContext());

    //Creating sequence and request.
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request_array[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&pSequence);
    apData->GetDescriptorData()->SetSequence(pSequence);
    apData->GetDescriptorData()->SetRequest(request_array);

    //Initialization
    T loglik = 0.0, logdet;
    T variance, variance1 = 1, variance2 = 1, variance3, dot_product = 0, dot_product1 = 0, dot_product2 = 0, dot_product3, n;
    T dzcpy_time, time_facto, time_solve, logdet_calculate, matrix_gen_time;

    double avg_executed_time_per_iteration = 0, avg_flops_per_iter = 0.0;

    int nhrs, success, i, num_params = 0;
    T flops = 0.0;
    T *univariate_theta;
    T *univariate2_theta;
    T *univariate3_theta;
    T nu12;
    T rho;
    T sigma_square12;

    auto kernel_name = apConfigurations->GetKernelName();
    auto median_locations = Locations<T>(1, apData->GetLocations()->GetDimension());
    apData->CalculateMedianLocations(kernel_name, median_locations);

    auto *CHAM_desc_C = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_C).chameleon_desc;
    auto *CHAM_desc_sub_C11 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_C11).chameleon_desc;
    auto *CHAM_desc_sub_C12 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_sub_C22 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_C22).chameleon_desc;

    auto *CHAM_desc_Z = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_desc_Z1 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_1).chameleon_desc;
    auto *CHAM_desc_Z2 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_2).chameleon_desc;
    auto *CHAM_desc_Zcpy = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_Z_COPY).chameleon_desc;

    auto *CHAM_desc_det = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_DETERMINANT).chameleon_desc;

    auto *CHAM_desc_product = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *CHAM_desc_product1 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                          DescriptorName::DESCRIPTOR_PRODUCT_1).chameleon_desc;
    auto *CHAM_desc_product2 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                          DescriptorName::DESCRIPTOR_PRODUCT_2).chameleon_desc;

    auto sequence = (RUNTIME_sequence_t *) apData->GetDescriptorData()->GetSequence();
    auto request = (RUNTIME_request_t *) apData->GetDescriptorData()->GetRequest();

    T *determinant = apData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_det);
    *determinant = 0;

    T *product = apData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_product);
    *product = 0;

    n = CHAM_desc_C->m;
    nhrs = CHAM_desc_Z->n;

    START_TIMING(dzcpy_time);
    int iter_count = apData->GetMleIterations();
    if (iter_count == 0) {
        // Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        ExaGeoStatLapackCopyTile(UpperLower::EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
    }
    string recovery_file = apConfigurations->GetRecoveryFile();

    if (recovery_file.empty() ||
        !(this->recover((char *) (recovery_file.c_str()), iter_count, (T *) theta, &loglik, num_params))) {
        if (iter_count == 0) {
            // Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
        } else {
            VERBOSE("Re-store the original Z vector...")
            ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Zcpy, CHAM_desc_Z);
            VERBOSE(" Done.\n")
        }
        STOP_TIMING(dzcpy_time);

        //Generate new co-variance matrix C based on new theta
        VERBOSE("Generate New Covariance Matrix...")
        START_TIMING(matrix_gen_time);

        if (kernel_name == "bivariate_matern_parsimonious2" ||
            kernel_name == "bivariate_matern_parsimonious2_profile") {

            int upper_lower = EXAGEOSTAT_UPPER_LOWER;
            kernels::Kernel<T> *kernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T>>::Create(
                    kernel_name);
            num_params = kernel->GetParametersNumbers();
            univariate_theta = (T *) malloc(3 * sizeof(T));
            univariate2_theta = (T *) malloc(3 * sizeof(T));
            univariate3_theta = (T *) malloc(3 * sizeof(T));
            univariate_theta[0] = theta[0];
            univariate_theta[1] = theta[2];
            univariate_theta[2] = theta[3];

            int distance_metric = apConfigurations->GetDistanceMetric();

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_sub_C11, upper_lower,
                                          apData->GetLocations(), apData->GetLocations(),
                                          &median_locations, univariate_theta,
                                          distance_metric, kernel);

            nu12 = 0.5 * (theta[3] + theta[4]);

            rho = theta[5] * sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) /
                                  (tgamma(theta[3]) * tgamma(theta[4]))) *
                  tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(theta[0] * theta[1]);

            univariate2_theta[0] = sigma_square12;
            univariate2_theta[1] = theta[2];
            univariate2_theta[2] = nu12;

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_sub_C12, upper_lower,
                                          &median_locations,
                                          apData->GetLocations(),
                                          &median_locations, univariate2_theta,
                                          0, kernel);

            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n")

            univariate3_theta[0] = theta[1];
            univariate3_theta[1] = theta[2];
            univariate3_theta[2] = theta[4];

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_sub_C22, upper_lower,
                                          &median_locations,
                                          apData->GetLocations(),
                                          &median_locations, univariate2_theta,
                                          0, kernel);
            delete kernel;
        } else {
            int upper_lower = EXAGEOSTAT_LOWER;

            kernels::Kernel<T> *kernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T>>::Create(
                    kernel_name);
            num_params = kernel->GetParametersNumbers();

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_C, upper_lower, apData->GetLocations(),
                                          apData->GetLocations(),
                                          &median_locations, (T *) theta,
                                          0, kernel);
            delete kernel;
        }

        ExaGeoStatSequenceWait(sequence);
        STOP_TIMING(matrix_gen_time);

        VERBOSE("Done.\n")
        VERBOSE("Cholesky factorization of Sigma...")
        START_TIMING(time_facto);
        success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_desc_C);
        //// TODO: Contact chameleon team
        FAILURE_LOGGER(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n")
        STOP_TIMING(time_facto);
        flops = flops + FLOPS_DPOTRF(n);
        VERBOSE(" Done.\n")

        //Calculate log(|C|) --> log(square(|L|))
        VERBOSE("Calculating the log determinant ...")
        START_TIMING(logdet_calculate);
        ExaGeoStatMeasureDetTileAsync(CHAM_desc_C, sequence, &request, CHAM_desc_det);
        ExaGeoStatSequenceWait(sequence);

        logdet = 2 * (*determinant);
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n")

        // Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("Solving the linear system ...\n")
        START_TIMING(time_solve);
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_C,
                           CHAM_desc_Z);
        STOP_TIMING(time_solve);
        flops = flops + FLOPS_DTRSM(ChamLeft, n, nhrs);
        VERBOSE(" Done.\n")

        //Calculate MLE likelihood
        VERBOSE("Calculating the MLE likelihood function ...")
        ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z, CHAM_desc_Z, 0, CHAM_desc_product);

        if (kernel_name == "bivariate_matern_parsimonious2_profile") {

            loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet -
                     (T) (n / 2.0) * log(2.0 * PI);
            ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0,
                               CHAM_desc_product1);
            ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0,
                               CHAM_desc_product2);

            variance1 = (1.0 / (n / 2)) * dot_product1;
            variance2 = (1.0 / (n / 2)) * dot_product2;

        } else if (kernel_name == "bivariate_matern_parsimonious_profile") {
            loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet -
                     (double) (n / 2.0) * log(2.0 * PI);

            ExaGeoStaStrideVectorTileAsync(CHAM_desc_Z, CHAM_desc_Z1, CHAM_desc_Z2, sequence, &request[0]);

            ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0,
                               CHAM_desc_product1);
            ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0,
                               CHAM_desc_product2);
            variance1 = (1.0 / (n / 2)) * dot_product1;
            variance2 = (1.0 / (n / 2)) * dot_product2;

        } else {
            dot_product = *product;
            loglik = -0.5 * dot_product - 0.5 * logdet - (double) (n / 2.0) * log(2.0 * PI);
        }
        VERBOSE(" Done.\n")
    }
    fprintf(stderr, " %3d- Model Parameters (", iter_count + 1);

    if (apConfigurations->GetLogger()) {
        ///should we be leaving the file opening and closing to the user?
        fprintf(apConfigurations->GetFileLogPath(), " %3d- Model Parameters (", iter_count + 1);
    }

    if ((apConfigurations->GetKernelName() == "bivariate_matern_parsimonious_profile") ||
        (apConfigurations->GetKernelName() == "bivariate_matern_parsimonious2_profile")) {
        fprintf(stderr, "%.8f, %.8f,", variance1, variance2);
        if (apConfigurations->GetLogger()) {
            fprintf(apConfigurations->GetFileLogPath(), "%.8f, %.8f,", variance1,
                    variance2);
        }
        i = 2;
    } else {
        i = 0;
    }

    for (; i < num_params; i++) {
        fprintf(stderr, "%.8f", theta[i]);
        if (i < num_params - 1) {
            fprintf(stderr, ",");
        }

        if (apConfigurations->GetLogger()) {
            fprintf(apConfigurations->GetFileLogPath(), "%.8f, ", theta[i]);
        }
    }
    fprintf(stderr, ")----> LogLi: %.18f\n", loglik);

    if (apConfigurations->GetLogger()) {
        fprintf(apConfigurations->GetFileLogPath(), ")----> LogLi: %.18f\n", loglik);
    }


    fprintf(stderr, " ---- Facto Time: %6.2f \n", time_facto);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);

    apData->SetMleIterations(apData->GetMleIterations() + 1);
    // for experiments
    if (Configurations::GetRunMode() == VERBOSE_MODE) {
        avg_executed_time_per_iteration += /*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
        avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    }

    return loglik;
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatLapackCopyTile(exageostat::common::UpperLower aUpperLower, void *apA,
                                                              void *apB) {
    return CHAMELEON_dlacpy_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int
ChameleonImplementationDense<T>::ExaGeoStatLapackToDescriptor(exageostat::common::UpperLower aUpperLower, void *apAf77,
                                                              int aLda, void *apA) {
    return CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apAf77, aLda, (CHAM_desc_t *) apA);
}


template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatSequenceWait(void *apSequence) {
    return CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apSequence);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatPotrfTile(exageostat::common::UpperLower aUpperLower, void *apA) {
    return CHAMELEON_dpotrf_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower,
                                                        common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA,
                                                        void *apB) {
    return CHAMELEON_dtrsm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans,
                                (cham_diag_t) aDiag, aAlpha, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha,
                                                        void *apA, void *apB, T aBeta, void *apC) {
    return CHAMELEON_dgemm_Tile((cham_trans_t) aTransA, (cham_trans_t) aTransB, aAlpha, (CHAM_desc_t *) apA,
                                (CHAM_desc_t *) apB, aBeta, (CHAM_desc_t *) apC);
}

template<typename T>
void
ChameleonImplementationDense<T>::CopyDescriptorZ(DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                 T *apDoubleVector) {

    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->mpContext,
                         (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    int m, m0;
    int tempmm;
    auto A = (CHAM_desc_t *) apDescriptor;
    struct starpu_codelet *cl = &this->cl_dzcpy;

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &apDoubleVector, sizeof(double),
                           STARPU_W, RUNTIME_data_getaddr(A, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dzcpy",
#endif
                           0);
    }
    RUNTIME_options_ws_free(&options);

}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                                   void *apDescDet) {

    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->mpContext,
                         (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);

    int m;
    int tempmm;
    auto Z = (CHAM_desc_t *) apDescA;
    auto det = (CHAM_desc_t *) apDescDet;
    struct starpu_codelet *cl = &this->cl_dmdet;

    for (m = 0; m < Z->mt; m++) {
        tempmm = m == Z->mt - 1 ? Z->m - m * Z->mb : Z->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, RUNTIME_data_getaddr(Z, m, m),
                           STARPU_RW, RUNTIME_data_getaddr(det, 0, 0),
                           0);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                                    void *apSequence, void *apRequest) {
    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->mpContext,
                         (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);

    int m, m0;
    int tempmm;
    auto A = (CHAM_desc_t *) apDescA;
    auto B = (CHAM_desc_t *) apDescB;
    auto C = (CHAM_desc_t *) apDescC;
    struct starpu_codelet *cl = &this->cl_stride_vec;
    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &m, sizeof(int),
                           STARPU_R, RUNTIME_data_getaddr(A, m, 0),
                           STARPU_W, RUNTIME_data_getaddr(B, (int) floor(m / 2.0), 0),
                           STARPU_W, RUNTIME_data_getaddr(C, (int) floor(m / 2.0), 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "stride_vec",
#endif
                           0);

    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
void ChameleonImplementationDense<T>::CovarianceMatrixCodelet(DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                              int &aTriangularPart,
                                                              dataunits::Locations<T> *apLocation1,
                                                              dataunits::Locations<T> *apLocation2,
                                                              dataunits::Locations<T> *apLocation3, T *aLocalTheta,
                                                              int aDistanceMetric, kernels::Kernel<T> *apKernel) {

    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->mpContext,
                         (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    int tempmm, tempnn;
    auto *CHAM_apDescriptor = (CHAM_desc_t *) apDescriptor;

    struct starpu_codelet *cl = &this->cl_dcmg;
    int m, n, m0 = 0, n0 = 0;

    for (n = 0; n < CHAM_apDescriptor->nt; n++) {
        tempnn = n == CHAM_apDescriptor->nt - 1 ? CHAM_apDescriptor->n - n * CHAM_apDescriptor->nb
                                                : CHAM_apDescriptor->nb;
        if (aTriangularPart == ChamUpperLower) {
            m = 0;
        } else {
            m = CHAM_apDescriptor->m == CHAM_apDescriptor->n ? n : 0;
        }
        for (; m < CHAM_apDescriptor->mt; m++) {

            tempmm = m == CHAM_apDescriptor->mt - 1 ? CHAM_apDescriptor->m - m * CHAM_apDescriptor->mb
                                                    : CHAM_apDescriptor->mb;
            m0 = m * CHAM_apDescriptor->mb;
            n0 = n * CHAM_apDescriptor->nb;
            // Register the data with StarPU
            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(kernels::Kernel<T> *),
                               0);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);

    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apDescriptorData->GetSequence());
}

template<typename T>
void
ChameleonImplementationDense<T>::ExaGeoStatGaussianToNonTileAsync(DescriptorData<T> *apDescriptorData, void *apDesc,
                                                                  T *apTheta) {

    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->mpContext,
                         (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    int m, m0;
    int tempmm;
    auto desc_z = (CHAM_desc_t *) apDesc;
    struct starpu_codelet *cl = &this->cl_gaussian_to_non;


    for (m = 0; m < desc_z->mt; m++) {
        tempmm = m == desc_z->mt - 1 ? desc_z->m - m * desc_z->mb : desc_z->mb;
        m0 = m * desc_z->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr(desc_z, m, 0),
                           STARPU_VALUE, &apTheta[0], sizeof(T),
                           STARPU_VALUE, &apTheta[1], sizeof(T),
                           STARPU_VALUE, &apTheta[2], sizeof(T),
                           STARPU_VALUE, &apTheta[3], sizeof(T),
                           STARPU_VALUE, &apTheta[4], sizeof(T),
                           STARPU_VALUE, &apTheta[5], sizeof(T),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "gaussian_to_non",
#endif
                           0);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apDescriptorData->GetSequence());
}
