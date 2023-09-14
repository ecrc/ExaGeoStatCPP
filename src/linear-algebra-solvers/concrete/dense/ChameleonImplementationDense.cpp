
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDense.cpp
 * @brief Dense Tile implementation of linear algebra methods.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <lapacke.h>

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>
#include <data-units/DescriptorData.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors(Configurations &aConfigurations,
                                                          DescriptorData<T> &aDescriptorData, T *apMeasurementsMatrix) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    // Get the problem size and other configuration parameters
    int n = aConfigurations.GetProblemSize();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = new T[n];

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    if (!aDescriptorData.GetSequence()) {
        RUNTIME_sequence_t *pSequence;
        CHAMELEON_Sequence_Create(&pSequence);
        aDescriptorData.SetSequence(pSequence);
        RUNTIME_request_t request_array[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
        aDescriptorData.SetRequest(request_array);
    }

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
                                  dts * dts, n, n, 0, 0, n, n, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, apMeasurementsMatrix, float_point,
                                  dts, dts, dts * dts, n, 1, 0, 0, n, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, nullptr, float_point, dts,
                                  dts, dts * dts, n, 1, 0, 0, n, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT, is_OOC, nullptr, float_point,
                                  dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT, is_OOC, nullptr, float_point, dts,
                                  dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    if (float_point == EXAGEOSTAT_REAL_DOUBLE) {
        auto *CHAM_descC = aDescriptorData.GetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C11, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, n, n, 0, 0, CHAM_descC->m / 2, CHAM_descC->n / 2, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, n, n, CHAM_descC->m / 2, 0, CHAM_descC->m / 2, CHAM_descC->n / 2,
                                      p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, n, n, CHAM_descC->m / 2, CHAM_descC->n / 2, CHAM_descC->m / 2,
                                      CHAM_descC->n / 2, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_1, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, n / 2, 1, 0, 0, n / 2, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_2, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, n / 2, 1, 0, 0, n / 2, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1, is_OOC, nullptr, float_point,
                                      dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2, is_OOC, nullptr, float_point,
                                      dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }

    //stop gsl error handler
    gsl_set_error_handler_off();
    aDescriptorData.SetIsDescriptorInitiated(true);
    delete[] Zcpy;
}

template<typename T>
void ChameleonImplementationDense<T>::InitiatePredictionDescriptors(
        Configurations &aConfigurations, ExaGeoStatData<T> &aData) {

    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    // Get the problem size and other configuration parameters
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();
    int z_miss_number = aConfigurations.GetUnknownObservationsNb();
    int n_z_obs = aConfigurations.CalculateZObsNumber();

    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else if (sizeof(T) == SIZE_OF_DOUBLE) {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    } else {
        throw runtime_error("Unsupported for now!");
    }
    aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_OBSERVATIONS, is_OOC, nullptr,
                                             float_point, dts, dts, dts * dts, n_z_obs, 1, 0, 0, n_z_obs, 1, p_grid,
                                             q_grid);

    //TODO: probably unnecessary check since it's only called when mspe is on.
    if (aConfigurations.GetIsMSPE()) {
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_Actual, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, z_miss_number, 1, 0, 0,
                                                 z_miss_number, 1,
                                                 p_grid, q_grid);
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MSE, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MSE_1, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MSE_2, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }

    aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_MISS, is_OOC, nullptr,
                                             float_point, dts, dts, dts * dts, z_miss_number, 1, 0, 0, z_miss_number, 1,
                                             p_grid,
                                             q_grid);
    try {
        descriptor::ExaGeoStatDescriptor<T> exaGeoStatDescriptor;
        exaGeoStatDescriptor.DestroyDescriptor(CHAMELEON_DESCRIPTOR,
                                               aData.GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                        DESCRIPTOR_C12).chameleon_desc);
        exaGeoStatDescriptor.DestroyDescriptor(CHAMELEON_DESCRIPTOR,
                                               aData.GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                        DESCRIPTOR_C22).chameleon_desc);
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, z_miss_number, n_z_obs, 0, 0,
                                                 z_miss_number,
                                                 n_z_obs, p_grid, q_grid);
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, n_z_obs, n_z_obs, 0, 0, n_z_obs,
                                                 n_z_obs,
                                                 p_grid, q_grid);

    } catch (...) {
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, z_miss_number, n_z_obs, 0, 0,
                                                 z_miss_number,
                                                 n_z_obs, p_grid, q_grid);
        aData.GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, is_OOC, nullptr,
                                                 float_point, dts, dts, dts * dts, n_z_obs, n_z_obs, 0, 0, n_z_obs,
                                                 n_z_obs,
                                                 p_grid, q_grid);
    }

}

template<typename T>
void ChameleonImplementationDense<T>::GenerateObservationsVector(Configurations &aConfigurations,
                                                                 DescriptorData<T> *apDescriptorData,
                                                                 BaseDescriptor aDescriptor, Locations<T> *apLocation1,
                                                                 Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                                 int aDistanceMetric) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    const int n = aConfigurations.GetProblemSize();
    int seed = aConfigurations.GetSeed();
    int iseed[4] = {seed, seed, seed, 1};
    auto *p_descriptor = aDescriptor.chameleon_desc;

    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *Nrand = new T[n];
    LAPACKE_dlarnv(3, iseed, n, (double *) Nrand);

    //Generate the co-variance matrix C
    auto *theta = new T[aConfigurations.GetInitialTheta().size()];
    for (int i = 0; i < aConfigurations.GetInitialTheta().size(); i++) {
        theta[i] = aConfigurations.GetInitialTheta()[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(apDescriptorData, p_descriptor, upper_lower, apLocation1, apLocation2, apLocation3,
                                  theta, aDistanceMetric, aConfigurations.GetKernelName());
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
    delete[] theta;
    const int P = aConfigurations.GetP();
    if (aConfigurations.GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = new T[n];
        CHAMELEON_Tile_to_Lapack( CHAM_descriptorZ, pMatrix, N);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &n, &P, aConfigurations->GetLoggerPath(), apLocation1);
        }
        delete[] pMatrix;
#else
        pMatrix = (T *) CHAM_descriptorZ->mat;
        string path = aConfigurations.GetLoggerPath();
        DiskWriter<T>::WriteVectorsToDisk(*pMatrix, n, P, path, *apLocation1);
#endif
        VERBOSE(" Done.\n")
    }

    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, p_descriptor);
    delete[] Nrand;
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
T ChameleonImplementationDense<T>::ExaGeoStatMleTile(const hardware::ExaGeoStatHardware &aHardware,
                                                     ExaGeoStatData<T> &aData, Configurations &aConfigurations,
                                                     const double *theta, T *apMeasurementsMatrix) {

    this->SetContext(aHardware.GetContext());
    if (!aData.GetDescriptorData()->GetIsDescriptorInitiated()) {
        this->InitiateDescriptors(aConfigurations, *aData.GetDescriptorData(), apMeasurementsMatrix);
    }

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    if (!aData.GetDescriptorData()->GetSequence()) {
        RUNTIME_sequence_t *sequence;
        CHAMELEON_Sequence_Create(&sequence);
        aData.GetDescriptorData()->SetSequence(sequence);
        aData.GetDescriptorData()->SetRequest(request_array);
    }
    auto pSequence = (RUNTIME_sequence_t *) aData.GetDescriptorData()->GetSequence();

    //Initialization
    T loglik = 0.0, logdet, variance, variance1 = 1, variance2 = 1, variance3, dot_product = 0, dot_product1 = 0, dot_product2 = 0, dot_product3, n, dzcpy_time, time_facto, time_solve, logdet_calculate, matrix_gen_time;
    double avg_executed_time_per_iteration = 0, avg_flops_per_iter = 0.0;

    int nhrs, success, i;
    T flops = 0.0;
    T *univariate_theta, *univariate2_theta, *univariate3_theta, nu12, rho, sigma_square12;

    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = kernels::KernelsConfigurations::GetParametersNumberKernelMap()[kernel_name];
    auto median_locations = Locations<T>(1, aData.GetLocations()->GetDimension());
    aData.CalculateMedianLocations(kernel_name, median_locations);

    auto *CHAM_desc_C = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_C).chameleon_desc;
    auto *CHAM_desc_sub_C11 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C11).chameleon_desc;
    auto *CHAM_desc_sub_C12 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_sub_C22 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C22).chameleon_desc;
    auto *CHAM_desc_Z = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_desc_Z1 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z_1).chameleon_desc;
    auto *CHAM_desc_Z2 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z_2).chameleon_desc;
    auto *CHAM_desc_Zcpy = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_COPY).chameleon_desc;
    auto *CHAM_desc_det = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_DETERMINANT).chameleon_desc;
    auto *CHAM_desc_product = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *CHAM_desc_product1 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT_1).chameleon_desc;
    auto *CHAM_desc_product2 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT_2).chameleon_desc;

    T *determinant = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_det);
    *determinant = 0;
    T *product = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_product);
    *product = 0;
    n = CHAM_desc_C->m;
    nhrs = CHAM_desc_Z->n;

    START_TIMING(dzcpy_time);
    int iter_count = aData.GetMleIterations();
    if (iter_count == 0) {
        // Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        ExaGeoStatLapackCopyTile(UpperLower::EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
    }
    string recovery_file = aConfigurations.GetRecoveryFile();

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
            univariate_theta = new T[3];
            univariate2_theta = new T[3];
            univariate3_theta = new T[3];
            univariate_theta[0] = theta[0];
            univariate_theta[1] = theta[2];
            univariate_theta[2] = theta[3];

            int distance_metric = aConfigurations.GetDistanceMetric();

            this->CovarianceMatrixCodelet(aData.GetDescriptorData(), CHAM_desc_sub_C11, upper_lower,
                                          aData.GetLocations(), aData.GetLocations(), &median_locations,
                                          univariate_theta, distance_metric, kernel_name);
            nu12 = 0.5 * (theta[3] + theta[4]);
            rho = theta[5] *
                  sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) / (tgamma(theta[3]) * tgamma(theta[4]))) *
                  tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(theta[0] * theta[1]);
            univariate2_theta[0] = sigma_square12;
            univariate2_theta[1] = theta[2];
            univariate2_theta[2] = nu12;

            this->CovarianceMatrixCodelet(aData.GetDescriptorData(), CHAM_desc_sub_C12, upper_lower, &median_locations,
                                          aData.GetLocations(), &median_locations, univariate2_theta, 0, kernel_name);
            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n")

            univariate3_theta[0] = theta[1];
            univariate3_theta[1] = theta[2];
            univariate3_theta[2] = theta[4];

            this->CovarianceMatrixCodelet(aData.GetDescriptorData(), CHAM_desc_sub_C22, upper_lower, &median_locations,
                                          aData.GetLocations(), &median_locations, univariate2_theta, 0, kernel_name);
        } else {
            int upper_lower = EXAGEOSTAT_LOWER;
            this->CovarianceMatrixCodelet(aData.GetDescriptorData(), CHAM_desc_C, upper_lower, aData.GetLocations(),
                                          aData.GetLocations(), &median_locations, (T *) theta, 0, kernel_name);
        }

        ExaGeoStatSequenceWait(pSequence);
        STOP_TIMING(matrix_gen_time);

        VERBOSE("Done.\n")
        VERBOSE("Cholesky factorization of Sigma...")
        START_TIMING(time_facto);
        success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_desc_C);
        //// TODO: Contact chameleon team
        FAILURE_LOGGER(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n")
        STOP_TIMING(time_facto);
        flops = flops + flops_dpotrf(n);
        VERBOSE(" Done.\n")

        //Calculate log(|C|) --> log(square(|L|))
        VERBOSE("Calculating the log determinant ...")
        START_TIMING(logdet_calculate);
        ExaGeoStatMeasureDetTileAsync(CHAM_desc_C, pSequence, &request_array, CHAM_desc_det);
        ExaGeoStatSequenceWait(pSequence);

        logdet = 2 * (*determinant);
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n")

        // Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("Solving the linear system ...\n")
        START_TIMING(time_solve);
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_C,
                           CHAM_desc_Z);
        STOP_TIMING(time_solve);
        flops = flops + flops_dtrsm(ChamLeft, n, nhrs);
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
            ExaGeoStaStrideVectorTileAsync(CHAM_desc_Z, CHAM_desc_Z1, CHAM_desc_Z2, pSequence, &request_array[0]);
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

    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), " %3d- Model Parameters (", iter_count + 1);
    }
    if ((aConfigurations.GetKernelName() == "bivariate_matern_parsimonious_profile") ||
        (aConfigurations.GetKernelName() == "bivariate_matern_parsimonious2_profile")) {
        fprintf(stderr, "%.8f, %.8f,", variance1, variance2);
        if (aConfigurations.GetLogger()) {
            fprintf(aConfigurations.GetFileLogPath(), "%.8f, %.8f,", variance1, variance2);
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

        if (aConfigurations.GetLogger()) {
            fprintf(aConfigurations.GetFileLogPath(), "%.8f, ", theta[i]);
        }
    }
    fprintf(stderr, ")----> LogLi: %.18f\n", loglik);

    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), ")----> LogLi: %.18f\n", loglik);
    }

    fprintf(stderr, " ---- Facto Time: %6.2f \n", time_facto);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);

    aData.SetMleIterations(aData.GetMleIterations() + 1);
    // for experiments
    if (Configurations::GetRunMode() == VERBOSE_MODE) {
        avg_executed_time_per_iteration += /*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
        avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    }

    return loglik;
}

template<typename T>
T *ChameleonImplementationDense<T>::ExaGeoStatMLEPredictTILE(exageostat::dataunits::ExaGeoStatData<T> &aData, T *apTheta, int aZMissNumber,
                                                             int aZObsNumber, T *apZObs, T *apZActual, T *apZMiss,
                                                             const hardware::ExaGeoStatHardware &aHardware,
                                                             configurations::Configurations &aConfiguration,
                                                             exageostat::dataunits::Locations<T> &aMissLocations,
                                                             exageostat::dataunits::Locations<T> &aObsLocations) {

    int i;
    this->SetContext(aHardware.GetContext());
    this->InitiatePredictionDescriptors(aConfiguration, aData);

    double time_solve, mat_gen_time, time_gemm, time_mse = 0.0, flops = 0.0;
    int num_params;

    auto *CHAM_desc_Zmiss = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_Z_MISS).chameleon_desc;
    auto *CHAM_desc_C12 = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_C22 = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_C22).chameleon_desc;
    auto *CHAM_desc_mse = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_MSE).chameleon_desc;
    auto *CHAM_desc_mse1 = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_MSE_1).chameleon_desc;
    auto *CHAM_desc_mse2 = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_MSE_2).chameleon_desc;
    auto *CHAM_desc_Zactual = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_Z_Actual).chameleon_desc;
    auto *CHAM_desc_Zobs = aData.GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_OBSERVATIONS).chameleon_desc;

    T *mspe = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mse);
    *mspe = 0;
    T *mspe_1 = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mse1);
    *mspe_1 = 0;
    T *mspe_2 = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mse2);
    *mspe_2 = 0;

    auto kernel_name = aConfiguration.GetKernelName();
    num_params = kernels::KernelsConfigurations::GetParametersNumberKernelMap()[kernel_name];
    auto median_locations = Locations<T>(1, aData.GetLocations()->GetDimension());
    aData.CalculateMedianLocations(kernel_name, median_locations);

    //Creating sequence and request.
    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    RUNTIME_sequence_t *sequence;
    if (!aData.GetDescriptorData()->GetSequence()) {
        CHAMELEON_Sequence_Create(&sequence);
        aData.GetDescriptorData()->SetSequence(sequence);
        aData.GetDescriptorData()->SetRequest(request_array);
    } else {
        sequence = (RUNTIME_sequence_t *) aData.GetDescriptorData()->GetSequence();
    }

    void *request = aData.GetDescriptorData()->GetRequest();

    //Copy data to vectors
    VERBOSE("Copy measurments vector to descZobs descriptor...")
    ExaGeoStatLap2Desc(apZObs, aZObsNumber, CHAM_desc_Zobs, UpperLower::EXAGEOSTAT_UPPER_LOWER);
    VERBOSE(" Done.\n")

    if (apZActual) {
        //Copy data to vectors
        VERBOSE("Copy actual measurments vector to descZactual descriptor...")
        ExaGeoStatLap2Desc(apZActual, aZMissNumber, CHAM_desc_Zactual, UpperLower::EXAGEOSTAT_UPPER_LOWER);
        VERBOSE(" Done.\n")
    }

    printf("estimated parameters:");
    for (i = 0; i < num_params; i++) {
        printf("%.8f,", apTheta[i]);
    }
    printf(")\n");

    ExaGeoStatSequenceWait(sequence);

    START_TIMING(mat_gen_time);
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(aData.GetDescriptorData(), CHAM_desc_C22, upper_lower, &aObsLocations, &aObsLocations,
                                  &median_locations, apTheta, 0, kernel_name);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE(" Done.\n")
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)")
    this->CovarianceMatrixCodelet(aData.GetDescriptorData(), CHAM_desc_C12, upper_lower, &aMissLocations,
                                  &aObsLocations, &median_locations, apTheta, 0, kernel_name);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE(" Done.\n")
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)")
    ExaGeoStatPosvTile(EXAGEOSTAT_LOWER, CHAM_desc_C22, CHAM_desc_Zobs);
    flops = flops + flops_dpotrf(aZObsNumber);
    flops = flops + flops_dtrsm(ChamLeft, aZObsNumber, aZObsNumber);
    VERBOSE(" Done.\n")
    STOP_TIMING(time_solve);


    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)")
    ExaGeoStatGemmTile(EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_C12, CHAM_desc_Zobs, 0, CHAM_desc_Zmiss);
    flops = flops + flops_dgemm(aZMissNumber, aZObsNumber, aZObsNumber);
    VERBOSE(" Done.\n")
    STOP_TIMING(time_gemm);

    ExaGeoStatDesc2Lap(apZMiss, aZMissNumber, CHAM_desc_Zmiss, EXAGEOSTAT_UPPER_LOWER);

    if (apZActual) {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n")
        if (kernel_name == "BivariateMaternParsimonious" || kernel_name == "BivariateMaternParsimonious2" ||
            kernel_name == "BivariateMaternParsimoniousProfile") {
            //TODO:  EXAGEOSTAT_MLE_dmse_bivariate_Tile_Async
            throw runtime_error("Bivariate Kernels are not supported yet.");
        } else {
            ExaGeoStatMleMseTileAsync(CHAM_desc_Zactual, CHAM_desc_Zmiss, CHAM_desc_mse, sequence, request);
        }
        ExaGeoStatSequenceWait(sequence);
        VERBOSE(" Done.\n")
        STOP_TIMING(time_mse);

        *mspe /= aZMissNumber;
        *mspe_1 /= aZMissNumber / 2;
        *mspe_2 /= aZMissNumber / 2;

    } else {
        *mspe = -1;
    }

    if (aConfiguration.GetLogger()) {
        fprintf(aConfiguration.GetFileLogPath(),
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, ""Flops: %.8f, Mean Square Error (MSE): %.8f\n\n",
                aZMissNumber, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), *mspe);
    }
    for (i = 0; i < aZMissNumber; i++)
        printf("(%3.6f, %3.6f)\n ", apZActual[i], apZMiss[i]);

    fprintf(stderr,
            "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, ""Flops: %.8f, Mean Square Error (MSE): %.8f\n\n",
            aZMissNumber, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), mspe[0]);

    T *all_mspe = new T[3];
    all_mspe[0] = *mspe;
    all_mspe[1] = *mspe_1;
    all_mspe[2] = *mspe_2;

    return all_mspe;
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
int
ChameleonImplementationDense<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha, void *apA,
                                                    void *apB, T aBeta, void *apC) {
    return CHAMELEON_dgemm_Tile((cham_trans_t) aTransA, (cham_trans_t) aTransB, aAlpha, (CHAM_desc_t *) apA,
                                (CHAM_desc_t *) apB, aBeta, (CHAM_desc_t *) apC);
}

template<typename T>
void
ChameleonImplementationDense<T>::CopyDescriptorZ(DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                 T *apDoubleVector) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
                                 this->mpContext,
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

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
            this->mpContext, (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);

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
    RUNTIME_options_finalize(&options, (CHAM_context_t *)
            this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                                    void *apSequence, void *apRequest) {
    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
            this->mpContext, (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);

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
    RUNTIME_options_finalize(&options, (CHAM_context_t *)
            this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
void ChameleonImplementationDense<T>::CovarianceMatrixCodelet(DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                              int &aTriangularPart, Locations<T> *apLocation1,
                                                              Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                              T *aLocalTheta, int aDistanceMetric,
                                                              const string &aKernelName) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
                                 this->mpContext, (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    kernels::Kernel<T> *pKernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T >>::Create(aKernelName);

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
                               STARPU_VALUE, &apLocation1, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocation2, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocation3, sizeof(Locations<T> *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &pKernel, sizeof(kernels::Kernel<T> *),
                               0);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apDescriptorData->GetSequence());
    delete pKernel;
}

template<typename T>
void
ChameleonImplementationDense<T>::ExaGeoStatGaussianToNonTileAsync(DescriptorData<T> *apDescriptorData, void *apDesc,
                                                                  T *apTheta) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
                                 this->mpContext, (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
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

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatMleMseTileAsync(void *apDescZPredict, void *apDescZMiss,
                                                               void *apDescError, void *apSequence,
                                                               void *apRequest) {
    // Check for Initialise the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
                                 this->mpContext,
                         (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);
    int m, m0;
    int tempmm;
    auto Zpre = (CHAM_desc_t *) apDescZPredict;
    struct starpu_codelet *cl = &this->cl_dmse;

    for (m = 0; m < Zpre->mt; m++) {
        tempmm = m == Zpre->mt - 1 ? Zpre->m - m * Zpre->mb : Zpre->mb;

        m0 = m * Zpre->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescError, 0, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZPredict, m, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZMiss, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dmse",
#endif
                           0);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
    ExaGeoStatSequenceWait(apSequence);
    return CHAMELEON_SUCCESS;
}

template<typename T>
int
ChameleonImplementationDense<T>::ExaGeoStatPosvTile(common::UpperLower aUpperLower, void *apA, void *apB) {
    return CHAMELEON_dposv_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}


template<typename T>
void
ChameleonImplementationDense<T>::ExaGeoStatLap2Desc(T *apA, int aLDA, void *apDescA, common::UpperLower aUpperLower) {
    CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, apA, aLDA, (CHAM_desc_t *) apDescA);
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatDesc2Lap(T *apA, int aLDA, void *apDescA,
                                                         common::UpperLower aUpperLower) {
    CHAMELEON_Desc2Lap((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apDescA, apA, aLDA);
}

template<typename T>
void ChameleonImplementationDense<T>::GetZObs(T *apZ, int aSize, DescriptorData<T> &aDescData) {
    auto z_desc = (CHAM_desc_t *) aDescData.GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
    ExaGeoStatDesc2Lap(apZ, aSize, z_desc, UpperLower::EXAGEOSTAT_UPPER_LOWER);
}