
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraMethods.cpp
 * @brief Implementation of linear algebra methods.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <lapacke.h>

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

using namespace std;

using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::dataunits;

// Define a method to set up the Chameleon descriptors
template<typename T>
void LinearAlgebraMethods<T>::InitiateDescriptors(Configurations &aConfigurations, DescriptorData<T> &aDescriptorData,
                                                  const int &aP, T *apMeasurementsMatrix) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    // Get the problem size and other configuration parameters
    int full_problem_size = aConfigurations.GetProblemSize() * aP;
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    if (!aDescriptorData.GetSequence()) {
        RUNTIME_sequence_t *pSequence;
        ExaGeoStatCreateSequence(&pSequence);
        aDescriptorData.SetSequence(pSequence);
        RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
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
                                  dts * dts, full_problem_size, full_problem_size, 0, 0, full_problem_size,
                                  full_problem_size, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, apMeasurementsMatrix, float_point,
                                  dts, dts, dts * dts, full_problem_size, 1, 0, 0, full_problem_size, 1, p_grid,
                                  q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, nullptr, float_point, dts,
                                  dts, dts * dts, full_problem_size, 1, 0, 0, full_problem_size, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT, is_OOC, nullptr, float_point,
                                  dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid, false);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT, is_OOC, nullptr, float_point, dts,
                                  dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid, false);

    if (float_point == EXAGEOSTAT_REAL_DOUBLE) {
        auto *CHAM_descC = aDescriptorData.GetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C11, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, full_problem_size, full_problem_size, 0, 0, CHAM_descC->m / 2,
                                      CHAM_descC->n / 2, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, full_problem_size, full_problem_size, CHAM_descC->m / 2, 0,
                                      CHAM_descC->m / 2, CHAM_descC->n / 2,
                                      p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, full_problem_size, full_problem_size, CHAM_descC->m / 2,
                                      CHAM_descC->n / 2, CHAM_descC->m / 2,
                                      CHAM_descC->n / 2, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_1, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, full_problem_size / 2, 1, 0, 0, full_problem_size / 2, 1, p_grid,
                                      q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_2, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, full_problem_size / 2, 1, 0, 0, full_problem_size / 2, 1, p_grid,
                                      q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1, is_OOC, nullptr, float_point,
                                      dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2, is_OOC, nullptr, float_point,
                                      dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }

    if (aConfigurations.GetIsNonGaussian()) {
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_SUM, is_OOC, nullptr, float_point,
                                      dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }

    //stop gsl error handler
    gsl_set_error_handler_off();
    aDescriptorData.SetIsDescriptorInitiated(true);
}

template<typename T>
void LinearAlgebraMethods<T>::InitiateFisherDescriptors(Configurations &aConfigurations,
                                                        dataunits::DescriptorData<T> &aDescriptorData) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    // Get the problem size and other configuration parameters
    int full_problem_size = aConfigurations.GetProblemSize();
    int num_params = kernels::KernelsConfigurations::GetParametersNumberKernelMap()[aConfigurations.GetKernelName()];
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    if (!aDescriptorData.GetSequence()) {
        RUNTIME_sequence_t *pSequence;
        ExaGeoStatCreateSequence(&pSequence);
        aDescriptorData.SetSequence(pSequence);
        RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
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

    if (!aDescriptorData.GetIsDescriptorInitiated()) {
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C, is_OOC, nullptr, float_point, dts,
                                      dts, dts * dts, full_problem_size, full_problem_size, 0, 0, full_problem_size,
                                      full_problem_size, p_grid, q_grid);
    }

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_A, is_OOC, nullptr, float_point, dts,
                                  dts, dts * dts, num_params, num_params, 0, 0, num_params, num_params, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_CJ, is_OOC, nullptr, float_point,
                                  dts, dts, dts * dts, full_problem_size, full_problem_size, 0, 0, full_problem_size,
                                  full_problem_size, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_CK, is_OOC, nullptr, float_point,
                                  dts, dts, dts * dts, full_problem_size, full_problem_size, 0, 0, full_problem_size,
                                  full_problem_size, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_RESULTS, is_OOC, nullptr,
                                  float_point, dts, dts, dts * dts, full_problem_size, full_problem_size, 0, 0,
                                  full_problem_size, full_problem_size, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_C_TRACE, is_OOC, nullptr,
                                  float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_C_DIAG, is_OOC, nullptr, float_point,
                                  dts, dts, dts * dts, full_problem_size, 1, 0, 0, full_problem_size, 1, p_grid,
                                  q_grid);
}

template<typename T>
void LinearAlgebraMethods<T>::InitiatePredictionDescriptors(
        Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData, const int &aP) {

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
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_OBSERVATIONS, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, n_z_obs, 1, 0, 0, n_z_obs, 1, p_grid,
                                              q_grid);
    if (aConfigurations.GetIsMSPE()) {
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_Actual, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, z_miss_number, 1, 0, 0,
                                                  z_miss_number, 1, p_grid, q_grid);
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MSPE, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MSPE_1, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MSPE_2, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_MISS, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, z_miss_number, 1, 0, 0, z_miss_number,
                                              1, p_grid, q_grid);

    descriptor::ExaGeoStatDescriptor<T> exaGeoStatDescriptor;
    auto *CHAM_descC12 = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12).chameleon_desc;
    if (CHAM_descC12) {
        exaGeoStatDescriptor.DestroyDescriptor(CHAMELEON_DESCRIPTOR, CHAM_descC12);
    }
    auto *CHAM_descC22 = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22).chameleon_desc;
    if (CHAM_descC22) {
        exaGeoStatDescriptor.DestroyDescriptor(CHAMELEON_DESCRIPTOR, CHAM_descC22);
    }
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, n_z_obs, n_z_obs, 0, 0, n_z_obs,
                                              n_z_obs, p_grid, q_grid);
    if (aConfigurations.GetIsNonGaussian()) {
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_R, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, n_z_obs, z_miss_number, 0, 0,
                                                  n_z_obs, z_miss_number, p_grid, q_grid);
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_R_COPY, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, n_z_obs, z_miss_number, 0, 0,
                                                  n_z_obs, z_miss_number, p_grid, q_grid);
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, z_miss_number, z_miss_number, 0, 0,
                                                  z_miss_number, z_miss_number, p_grid, q_grid);
    } else {
        aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, is_OOC, nullptr,
                                                  float_point, dts, dts, dts * dts, z_miss_number, n_z_obs, 0, 0,
                                                  z_miss_number, n_z_obs, p_grid, q_grid);
    }
}

template<typename T>
void LinearAlgebraMethods<T>::InitiateMLOEMMOMDescriptors(Configurations &aConfigurations,
                                                          std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                          const int &aP) {

    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    int n_z_obs = aConfigurations.CalculateZObsNumber();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();

    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else if (sizeof(T) == SIZE_OF_DOUBLE) {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    } else {
        throw runtime_error("Unsupported for now!");
    }
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_k_T, is_OOC, nullptr,
                                              float_point,
                                              dts, dts, dts * dts, n_z_obs, aP, 0, 0, n_z_obs, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_k_A, is_OOC, nullptr,
                                              float_point,
                                              dts, dts, dts * dts, n_z_obs, aP, 0, 0, n_z_obs, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_k_A_TMP, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, n_z_obs, aP, 0, 0, n_z_obs, aP,
                                              p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_k_T_TMP, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, n_z_obs, aP, 0, 0, n_z_obs, aP,
                                              p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_EXPR_1, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, aP, aP, 0, 0, aP, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_EXPR_2, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, aP, aP, 0, 0, aP, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_EXPR_3, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, aP, aP, 0, 0, aP, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_EXPR_4, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, aP, aP, 0, 0, aP, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MLOE, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_MMOM, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_TRUTH_ALPHA, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, aP, aP, 0, 0, aP, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_TIMATED_ALPHA, is_OOC, nullptr,
                                              float_point, dts, dts, dts * dts, aP, aP, 0, 0, aP, aP, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_K_T, is_OOC, nullptr,
                                              float_point,
                                              dts, dts, dts * dts, n_z_obs, n_z_obs, 0, 0, n_z_obs,
                                              n_z_obs, p_grid, q_grid);
    aData->GetDescriptorData()->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_K_A, is_OOC, nullptr,
                                              float_point,
                                              dts, dts, dts * dts, n_z_obs, n_z_obs, 0, 0, n_z_obs,
                                              n_z_obs, p_grid, q_grid);
    //stop gsl error handler
    gsl_set_error_handler_off();
}

template<typename T>
void LinearAlgebraMethods<T>::GenerateSyntheticData(Configurations &aConfigurations,
                                                    const ExaGeoStatHardware &aHardware,
                                                    std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                    const kernels::Kernel<T> &aKernel) {

    this->mpContext = aHardware.GetChameleonContext();
    this->InitiateDescriptors(aConfigurations, *aData->GetDescriptorData(), aKernel.GetVariablesNumber());
    auto median_locations = Locations<T>(1, aData->GetLocations()->GetDimension());
    this->GenerateObservationsVector(aConfigurations, aData, aData->GetLocations(), aData->GetLocations(),
                                     &median_locations, aConfigurations.GetDistanceMetric(), aKernel);
}

template<typename T>
void LinearAlgebraMethods<T>::GenerateObservationsVector(Configurations &aConfigurations,
                                                         std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                         Locations<T> *apLocation1, Locations<T> *apLocation2,
                                                         Locations<T> *apLocation3, const int &aDistanceMetric,
                                                         const kernels::Kernel<T> &aKernel) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    const int P = aKernel.GetVariablesNumber();
    const int full_problem_size = aConfigurations.GetProblemSize() * P;
    int seed = aConfigurations.GetSeed();
    int initial_seed[4] = {seed, seed, seed, 1};
    T time_facto = 0.0, time_trmm = 0.0, matrix_gen_time = 0.0, flops = 0;

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    RUNTIME_sequence_t *sequence;
    if (!aData->GetDescriptorData()->GetSequence()) {
        ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    } else {
        sequence = (RUNTIME_sequence_t *) aData->GetDescriptorData()->GetSequence();
    }
    auto *CHAM_descC = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
    //normal random generation of e -- ei~N(0, 1) to generate Z
    auto *randomN = new T[full_problem_size];
    LAPACKE_dlarnv(3, initial_seed, full_problem_size, (double *) randomN);

    //Generate the co-variance matrix C
    auto *theta = new T[aConfigurations.GetInitialTheta().size()];
    for (int i = 0; i < aConfigurations.GetInitialTheta().size(); i++) {
        theta[i] = aConfigurations.GetInitialTheta()[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    START_TIMING(matrix_gen_time);

    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_descC, upper_lower, apLocation1, apLocation2,
                                  apLocation3, theta, aDistanceMetric, &aKernel);
    ExaGeoStatSequenceWait(sequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE("Done.")

    //Copy randomN to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto *CHAM_descZ = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
    CopyDescriptorZ(*aData->GetDescriptorData(), CHAM_descZ, randomN);
    VERBOSE("Done.")

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")
    START_TIMING(time_facto);
    ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_descC, 0, nullptr, nullptr, 0, 0);
    STOP_TIMING(time_facto);
    flops = flops + flops_dpotrf(full_problem_size);
    VERBOSE("Done.")

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
    START_TIMING(time_trmm);
    ExaGeoStatTrmmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_descC,
                       CHAM_descZ);
    STOP_TIMING(time_trmm);
    flops = flops + flops_dtrmm(ChamLeft, full_problem_size, CHAM_descZ->n);
    VERBOSE("Done.")

    if (aConfigurations.GetIsNonGaussian()) {
        //Gaussian to non-gaussian transformation
        VERBOSE("Convert Z Gaussian to non-Gaussian (Synthetic Dataset Generation Phase) .....")
        ExaGeoStatGaussianToNonTileAsync(*aData->GetDescriptorData(), CHAM_descZ, theta);
        VERBOSE("Done.")
    }
    delete[] theta;
    if (aConfigurations.GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = new T[full_problem_size];
        string path = aConfigurations.GetLoggerPath();
        ExaGeoStatDesc2Lap(pMatrix, full_problem_size, CHAM_descZ, EXAGEOSTAT_UPPER_LOWER);
        if ( CHAMELEON_Comm_rank() == 0 ){
            helpers::DiskWriter<T>::WriteVectorsToDisk(*pMatrix, full_problem_size, P, path, *apLocation1);
        }
        delete[] pMatrix;
#else
        pMatrix = (T *) CHAM_descZ->mat;
        string path = aConfigurations.GetLoggerPath();
        helpers::DiskWriter<T>::WriteVectorsToDisk(*pMatrix, full_problem_size, P, path, *apLocation1);
#endif
        VERBOSE("Done.")
    }

    ExaGeoStatLaSetTile(EXAGEOSTAT_UPPER_LOWER, 0, 0, CHAM_descC);
    delete[] randomN;
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")

    int total_flops = flops / 1e9 / (time_facto + time_trmm);
    LOGGER(" ---- Facto Time: " << time_facto)
    LOGGER(" ---- dtrmm Time: " << time_trmm)
    LOGGER(" ---- Matrix Generation Time: " << matrix_gen_time)
    LOGGER(" ---- Total Time: " << time_facto + time_trmm)
    LOGGER(" ---- Gflop/s: " << total_flops)

    results::Results::GetInstance()->SetTotalDataGenerationExecutionTime(time_facto + time_trmm);
    results::Results::GetInstance()->SetTotalDataGenerationFlops(total_flops);
    results::Results::GetInstance()->SetGeneratedLocationsNumber(full_problem_size / aConfigurations.GetTimeSlot());
    results::Results::GetInstance()->SetIsLogger(aConfigurations.GetLogger());
    results::Results::GetInstance()->SetLoggerPath(aConfigurations.GetLoggerPath());

}

template<typename T>
T *LinearAlgebraMethods<T>::ExaGeoStatMLEPredictTile(std::unique_ptr<ExaGeoStatData<T>> &aData, T *apTheta,
                                                     const int &aZMissNumber, const int &aZObsNumber, T *apZObs,
                                                     T *apZActual, T *apZMiss, const ExaGeoStatHardware &aHardware,
                                                     Configurations &aConfiguration, Locations<T> &aMissLocations,
                                                     Locations<T> &aObsLocations, const kernels::Kernel<T> &aKernel) {

    int i;
    this->SetContext(aHardware.GetChameleonContext());
    this->InitiatePredictionDescriptors(aConfiguration, aData, aKernel.GetVariablesNumber());

    double time_solve, mat_gen_time, time_gemm, time_mspe = 0.0, flops = 0.0;
    int num_params;

    auto *CHAM_desc_Zmiss = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_Z_MISS).chameleon_desc;
    auto *CHAM_desc_C12 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_C22 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_C22).chameleon_desc;
    auto *CHAM_desc_mspe = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_MSPE).chameleon_desc;
    auto *CHAM_desc_mspe1 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_MSPE_1).chameleon_desc;
    auto *CHAM_desc_mspe2 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_MSPE_2).chameleon_desc;
    auto *CHAM_desc_Zactual = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_Z_Actual).chameleon_desc;
    auto *CHAM_desc_Zobs = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_Z_OBSERVATIONS).chameleon_desc;

    T *mspe = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mspe);
    *mspe = 0;
    T *mspe_1 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mspe1);
    *mspe_1 = 0;
    T *mspe_2 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mspe2);
    *mspe_2 = 0;

    auto kernel_name = aConfiguration.GetKernelName();
    num_params = aKernel.GetParametersNumbers();
    auto median_locations = Locations<T>(1, aData->GetLocations()->GetDimension());
    aData->CalculateMedianLocations(kernel_name, median_locations);

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    RUNTIME_sequence_t *sequence;
    if (!aData->GetDescriptorData()->GetSequence()) {
        ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    } else {
        sequence = (RUNTIME_sequence_t *) aData->GetDescriptorData()->GetSequence();
    }
    void *request = aData->GetDescriptorData()->GetRequest();

    //Copy data to vectors
    VERBOSE("Copy measurements vector to descZobs descriptor...")
    ExaGeoStatLap2Desc(apZObs, aZObsNumber, CHAM_desc_Zobs, UpperLower::EXAGEOSTAT_UPPER_LOWER);
    VERBOSE("Done.")

    if (apZActual) {
        //Copy data to vectors
        VERBOSE("Copy actual measurements vector to descZactual descriptor...")
        ExaGeoStatLap2Desc(apZActual, aZMissNumber, CHAM_desc_Zactual, UpperLower::EXAGEOSTAT_UPPER_LOWER);
        VERBOSE("Done.")
    }

    LOGGER("- Estimated Parameters (", true)
    for (i = 0; i < num_params; i++) {
        LOGGER_PRECISION(apTheta[i])
        if (i != num_params - 1) {
            LOGGER_PRECISION(", ")
        }
    }
    LOGGER_PRECISION(")")
    LOGGER("")

    START_TIMING(mat_gen_time);
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_C22, upper_lower, &aObsLocations,
                                  &aObsLocations, &median_locations, apTheta, 0, &aKernel);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE("Done.")
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)")
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_C12, upper_lower, &aMissLocations,
                                  &aObsLocations, &median_locations, apTheta, 0, &aKernel);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE("Done.")
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)")
    ExaGeoStatPosvTile(EXAGEOSTAT_LOWER, CHAM_desc_C22, CHAM_desc_Zobs);
    flops = flops + flops_dpotrf(aZObsNumber);
    flops = flops + flops_dtrsm(ChamLeft, aZObsNumber, aZObsNumber);
    VERBOSE("Done.")
    STOP_TIMING(time_solve);

    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)")
    CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_desc_C12, CHAM_desc_Zobs, 0, CHAM_desc_Zmiss);
    flops = flops + flops_dgemm(aZMissNumber, aZObsNumber, aZObsNumber);
    VERBOSE("Done.")
    STOP_TIMING(time_gemm);
    ExaGeoStatDesc2Lap(apZMiss, aZMissNumber, CHAM_desc_Zmiss, EXAGEOSTAT_UPPER_LOWER);

    if (apZActual) {
        START_TIMING(time_mspe);
        VERBOSE("Calculate Mean Square Prediction Error (MSPE) ... (Prediction Stage)")
        if (kernel_name == "BivariateMaternParsimonious") {
            ExaGeoStatMLEMSPEBivariateTileAsync(CHAM_desc_Zactual, CHAM_desc_Zmiss, CHAM_desc_mspe1, CHAM_desc_mspe2,
                                                CHAM_desc_mspe, sequence, &request_array[0]);
        } else {
            this->ExaGeoStatMLEMSPETileAsync(CHAM_desc_Zactual, CHAM_desc_Zmiss, CHAM_desc_mspe, sequence, request);
        }
        ExaGeoStatSequenceWait(sequence);
        VERBOSE("Done.")
        STOP_TIMING(time_mspe);

        *mspe /= aZMissNumber;
        *mspe_1 /= aZMissNumber / 2;
        *mspe_2 /= aZMissNumber / 2;

    } else {
        *mspe = -1;
    }

    if (aConfiguration.GetLogger()) {
        fprintf(aConfiguration.GetFileLogPath(),
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, ""Flops: %.8f, Mean Square Prediction Error (MSPE): %.8f\n\n",
                aZMissNumber, (mat_gen_time + time_solve + time_mspe), (flops / 1e9 / (time_solve)), *mspe);
    }
    VERBOSE("- Z Actual .. Z Miss")
    for (i = 0; i < aZMissNumber; i++) {
        VERBOSE(" (" << apZActual[i] << ", " << apZMiss[i] << ")")
    }

    results::Results::GetInstance()->SetMSPEExecutionTime(time_solve + time_gemm);
    results::Results::GetInstance()->SetMSPEFlops((flops / 1e9 / (time_solve + time_gemm)));
    results::Results::GetInstance()->SetMSPEError(*mspe);

    T *all_mspe = new T[3];
    all_mspe[0] = *mspe;
    all_mspe[1] = *mspe_1;
    all_mspe[2] = *mspe_2;

    return all_mspe;
}

template<typename T>
T *LinearAlgebraMethods<T>::ExaGeoStatMLENonGaussianPredictTile(std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                                T *apTheta, const int &aZMissNumber,
                                                                const int &aZObsNumber, T *apZObs, T *apZActual,
                                                                T *apZMiss,
                                                                const ExaGeoStatHardware &aHardware,
                                                                Configurations &aConfiguration,
                                                                dataunits::Locations<T> &aMissLocations,
                                                                dataunits::Locations<T> &aObsLocations,
                                                                const kernels::Kernel<T> &aKernel) {

    int i;
    this->SetContext(aHardware.GetChameleonContext());
    this->InitiatePredictionDescriptors(aConfiguration, aData, aKernel.GetVariablesNumber());

    double time_solve, mat_gen_time, time_mse, mat_gen_time_2, dposv_time, gemms_time, time_trsm, time_gemm, time_mspe = 0.0, flops = 0.0;
    int num_params;

    auto *CHAM_desc_Zmiss = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_Z_MISS).chameleon_desc;
    auto *CHAM_desc_C12 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_C22 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_C22).chameleon_desc;
    auto *CHAM_desc_mspe = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_MSPE).chameleon_desc;
    auto *CHAM_desc_mspe1 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_MSPE_1).chameleon_desc;
    auto *CHAM_desc_mspe2 = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_MSPE_2).chameleon_desc;
    auto *CHAM_desc_Zactual = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_Z_Actual).chameleon_desc;
    auto *CHAM_desc_Zobs = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_Z_OBSERVATIONS).chameleon_desc;
    auto *CHAM_desc_R = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_R).chameleon_desc;
    auto *CHAM_desc_Rcopy = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_R_COPY).chameleon_desc;
    auto *CHAM_desc_product = aData->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;

    T *mspe = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mspe);
    *mspe = 0;
    T *mspe_1 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mspe1);
    *mspe_1 = 0;
    T *mspe_2 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mspe2);
    *mspe_2 = 0;

    auto kernel_name = aConfiguration.GetKernelName();
    num_params = aKernel.GetParametersNumbers();
    auto median_locations = Locations<T>(1, aData->GetLocations()->GetDimension());
    aData->CalculateMedianLocations(kernel_name, median_locations);

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    RUNTIME_sequence_t *sequence;
    if (!aData->GetDescriptorData()->GetSequence()) {
        ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    } else {
        sequence = (RUNTIME_sequence_t *) aData->GetDescriptorData()->GetSequence();
    }
    void *request = aData->GetDescriptorData()->GetRequest();

    //Copy data to vectors
    VERBOSE("Copy measurements vector to descZobs descriptor...")
    ExaGeoStatLap2Desc(apZObs, aZObsNumber, CHAM_desc_Zobs, UpperLower::EXAGEOSTAT_UPPER_LOWER);
    VERBOSE("Done.")

    if (apZActual) {
        //Copy data to vectors
        VERBOSE("Copy actual measurements vector to descZactual descriptor...")
        ExaGeoStatLap2Desc(apZActual, aZMissNumber, CHAM_desc_Zactual, UpperLower::EXAGEOSTAT_UPPER_LOWER);
        VERBOSE("Done.")
    }

    LOGGER("- Estimated Parameters (", true)
    for (i = 0; i < num_params; i++) {
        LOGGER_PRECISION(apTheta[i])
        if (i != num_params - 1) {
            LOGGER_PRECISION(", ")
        }
    }
    LOGGER_PRECISION(")")
    LOGGER("")

    //Convert the non-Gaussian observation to Gaussian
    this->ExaGeoStatNonGaussianTransformTileAsync(CHAM_desc_Zobs, CHAM_desc_product, apTheta, sequence,
                                                  &request_array[0]);

    START_TIMING(mat_gen_time);
    VERBOSE("Generate R_theta Covariance Matrix... (Prediction Stage)")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_C22, upper_lower, &aObsLocations,
                                  &aObsLocations, &median_locations, apTheta, 0, &aKernel);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE("Done.")
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    VERBOSE("Calculate dposv R_theta Covariance Matrix... (Prediction Stage)");
    ExaGeoStatPosvTile(EXAGEOSTAT_LOWER, CHAM_desc_C22, CHAM_desc_Zobs);
    flops = flops + flops_dpotrf(aZObsNumber);
    flops = flops + 2 * flops_dtrsm(ChamLeft, CHAM_desc_C22->m, CHAM_desc_Zobs->n);
    VERBOSE("Done.")
    STOP_TIMING(time_solve);

    START_TIMING(mat_gen_time_2);
    VERBOSE("Generate R_theta Covariance Matrix... (Prediction Stage)")
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_R, upper_lower, &aObsLocations,
                                  &aMissLocations, &median_locations, apTheta, 0, &aKernel);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE("Done.")
    STOP_TIMING(mat_gen_time_2);

    ExaGeoStatLapackCopyTile(EXAGEOSTAT_LOWER, CHAM_desc_R, CHAM_desc_Rcopy);

    START_TIMING(time_trsm);
    VERBOSE("Calculate dposv r_theta Covariance Matrix... (Prediction Stage)");
    ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_C22,
                       nullptr, nullptr, CHAM_desc_R, 0);
    ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_C22,
                       nullptr, nullptr, CHAM_desc_R, 0);
    flops = flops + 2 * flops_dtrsm(ChamLeft, CHAM_desc_C22->m, CHAM_desc_Zobs->n);
    VERBOSE("Done.")
    STOP_TIMING(time_trsm);

    START_TIMING(time_gemm);
    VERBOSE("For each missing location, Generate correlation vector CHAMELEON_descr (Prediction Stage) .....");
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Rcopy, CHAM_desc_Zobs, 0, CHAM_desc_Zmiss);
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Rcopy, CHAM_desc_R, 0, CHAM_desc_C12);
    STOP_TIMING(time_gemm);

    auto r = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_C12);
    auto Zmiss2 = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_Zmiss);;

    auto ng_nu = new T[aZMissNumber];
    auto ng_sigma_sq = new T[aZMissNumber];

    double mu = -1;
    double sigma_sq = -1;

    for (i = 0; i < aZMissNumber; i++) {
        mu = Zmiss2[i];
        sigma_sq = 1 - r[i * aZMissNumber + i];
        ng_nu[i] = mu;
        ng_sigma_sq[i] = sigma_sq;
        //r^t X Z
        apZMiss[i] = apTheta[2] + (apTheta[3] / (apTheta[4] * sqrt(1 - apTheta[5] * sigma_sq))) *
                                  exp(apTheta[5] * pow(mu, 2) / (2 * (1 - apTheta[5] * sigma_sq))) *
                                  (exp((pow(apTheta[4], 2) * sigma_sq + 2 * apTheta[4] * mu) /
                                       (2 * (1 - apTheta[5] * sigma_sq))) - 1);
    }
    VERBOSE(" Done.")

    ExaGeoStatLap2Desc(apZMiss, aZMissNumber, CHAM_desc_Zmiss, EXAGEOSTAT_UPPER_LOWER);

    if (apZActual != nullptr) {
        START_TIMING(time_mspe);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
        this->ExaGeoStatMLEMSPETileAsync(CHAM_desc_Zactual, CHAM_desc_Zmiss, CHAM_desc_mspe, sequence, request);
        ExaGeoStatSequenceWait(sequence);
        VERBOSE(" Done.")
        STOP_TIMING(time_mspe);
        *mspe /= aZMissNumber;
    } else {
        *mspe = -1;
    }

    if(helpers::CommunicatorMPI::GetInstance()->GetRank() == 0) {
        if (aConfiguration.GetLogger()) {
            fprintf(aConfiguration.GetFileLogPath(),
                    "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, ""Flops: %.8f, Mean Square Prediction Error (MSPE): %.8f\n\n",
                    aZMissNumber, (mat_gen_time + mat_gen_time_2 + time_solve + time_mspe),
                    (flops / 1e9 / (time_solve)),
                    *mspe);
        }
        VERBOSE("- Z Actual .. Z Miss")
        for (i = 0; i < aZMissNumber; i++) {
            VERBOSE(" (" << apZActual[i] << ", " << apZMiss[i] << ")")
        }

        results::Results::GetInstance()->SetMSPEExecutionTime(time_solve + time_gemm + time_trsm);
        results::Results::GetInstance()->SetMSPEFlops((flops / 1e9 / (time_solve + time_gemm)));
        results::Results::GetInstance()->SetMSPEError(*mspe);
    }

    T *all_mspe = new T[3];
    all_mspe[0] = *mspe;
    all_mspe[1] = *mspe_1;
    all_mspe[2] = *mspe_2;
    delete[] ng_sigma_sq;
    delete[] ng_nu;
    return all_mspe;
}

template<typename T>
void LinearAlgebraMethods<T>::ExaGeoStatMLETileMLOEMMOM(Configurations &aConfigurations,
                                                        std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                        const ExaGeoStatHardware &aHardware, T *apTruthTheta,
                                                        T *apEstimatedTheta, Locations<T> &aMissLocations,
                                                        Locations<T> &aObsLocations,
                                                        const kernels::Kernel<T> &aKernel) {

    this->SetContext(aHardware.GetChameleonContext());
    this->InitiateMLOEMMOMDescriptors(aConfigurations, aData, aKernel.GetVariablesNumber());
    auto kernel_name = aConfigurations.GetKernelName();
    auto median_locations = Locations<T>(1, aData->GetLocations()->GetDimension());
    aData->CalculateMedianLocations(kernel_name, median_locations);

    int n_z_miss = aConfigurations.GetUnknownObservationsNb();
    int num_par = aKernel.GetParametersNumbers();
    LOGGER("- Truth Theta: ", true)
    for (int num = 0; num < num_par; num++) {
        LOGGER_PRECISION(apTruthTheta[num] << " ")
    }
    LOGGER("")
    LOGGER("- Estimated Theta: ", true)
    for (int num = 0; num < num_par; num++) {
        LOGGER_PRECISION(apEstimatedTheta[num] << " ")
    }
    LOGGER("")
    int p;
    double all_time, cholesky1, cholesky2, matrix_gen, vecs_gen, copy_vecs, trsm1, trsm2, trsm3, trsm4, gevv1, gevv2, gevv3, gevv4, gevv5;

    auto loe = new T[n_z_miss];
    auto mom = new T[n_z_miss];

    auto *CHAM_desc_k_t = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                    DESCRIPTOR_k_T).chameleon_desc;
    auto *CHAM_desc_k_a = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                    DESCRIPTOR_k_A).chameleon_desc;
    auto *CHAM_desc_K_t = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                    DESCRIPTOR_K_T).chameleon_desc;
    auto *CHAM_desc_K_a = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                    DESCRIPTOR_K_A).chameleon_desc;
    auto *CHAM_desc_k_a_tmp = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                        DESCRIPTOR_k_A_TMP).chameleon_desc;
    auto *CHAM_desc_k_t_tmp = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                        DESCRIPTOR_k_T_TMP).chameleon_desc;
    auto *CHAM_desc_expr1 = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                      DESCRIPTOR_EXPR_1).chameleon_desc;
    auto *CHAM_desc_expr2 = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                      DESCRIPTOR_EXPR_2).chameleon_desc;
    auto *CHAM_desc_expr3 = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                      DESCRIPTOR_EXPR_3).chameleon_desc;
    auto *CHAM_desc_expr4 = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                      DESCRIPTOR_EXPR_4).chameleon_desc;
    auto *CHAM_desc_mloe = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                     DESCRIPTOR_MLOE).chameleon_desc;
    auto *CHAM_desc_mmom = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                     DESCRIPTOR_MMOM).chameleon_desc;
    auto *CHAM_desc_estimated_alpha = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                DESCRIPTOR_TIMATED_ALPHA).chameleon_desc;
    auto *CHAM_desc_truth_alpha = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                            DESCRIPTOR_TRUTH_ALPHA).chameleon_desc;

    T *mloe = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mloe);
    *mloe = 0;
    T *mmom = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_mmom);
    *mmom = 0;

    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    RUNTIME_sequence_t *sequence;
    if (!aData->GetDescriptorData()->GetSequence()) {
        ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    } else {
        sequence = (RUNTIME_sequence_t *) aData->GetDescriptorData()->GetSequence();
    }
    void *request = aData->GetDescriptorData()->GetRequest();

    auto lmiss = new Locations<T>(n_z_miss, aData->GetLocations()->GetDimension());
    T nu12;
    T rho;
    T sigma_square12;

    T flops = 0.0;
    START_TIMING(all_time);

    int m = CHAM_desc_estimated_alpha->m;

    auto truth_alpha = new T[m * m];
    auto estimated_alpha = new T[m * m];
    auto temp1 = new T[m * m];
    auto temp2 = new T[m * m];
    auto temp3 = new T[m * m];

    if (m == 1) {
        truth_alpha[0] = apTruthTheta[0];
        estimated_alpha[0] = apEstimatedTheta[0];
    }

    if (m == 2) {
        double truth_nu12 = 0.5 * (apTruthTheta[3] + apTruthTheta[4]);
        double truth_rho = apTruthTheta[5] * sqrt((tgamma(apTruthTheta[3] + 1) * tgamma(apTruthTheta[4] + 1)) /
                                                  (tgamma(apTruthTheta[3]) * tgamma(apTruthTheta[4]))) *
                           tgamma(truth_nu12) / tgamma(truth_nu12 + 1);
        double estimated_nu12 = 0.5 * (apEstimatedTheta[3] + apEstimatedTheta[4]);
        double estimated_rho = apEstimatedTheta[5] *
                               sqrt((tgamma(apEstimatedTheta[3] + 1) * tgamma(apEstimatedTheta[4] + 1)) /
                                    (tgamma(apEstimatedTheta[3]) * tgamma(apEstimatedTheta[4]))) *
                               tgamma(estimated_nu12) / tgamma(estimated_nu12 + 1);

        truth_alpha[0] = apTruthTheta[0];
        estimated_alpha[0] = apEstimatedTheta[0];
        truth_alpha[1] = truth_alpha[3] = truth_rho * sqrt(apTruthTheta[0] * apTruthTheta[1]);
        estimated_alpha[1] = estimated_alpha[3] = estimated_rho * sqrt(apEstimatedTheta[0] * apEstimatedTheta[1]);
        truth_alpha[2] = apTruthTheta[1];
        estimated_alpha[2] = apEstimatedTheta[1];
    }

    this->ExaGeoStatLap2Desc(truth_alpha, m, CHAM_desc_truth_alpha, EXAGEOSTAT_UPPER_LOWER);
    this->ExaGeoStatLap2Desc(estimated_alpha, m, CHAM_desc_estimated_alpha, EXAGEOSTAT_UPPER_LOWER);

    START_TIMING(matrix_gen);
    VERBOSE("Create K_a and K_t Covariance Matrices (MLOE-MMOM).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_K_a, upper_lower, &aObsLocations,
                                  &aObsLocations, &median_locations, apEstimatedTheta, 0, &aKernel);
    this->ExaGeoStatSequenceWait(sequence);
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_K_t, upper_lower, &aObsLocations,
                                  &aObsLocations, &median_locations, apTruthTheta, 0, &aKernel);
    this->ExaGeoStatSequenceWait(sequence);
    VERBOSE("Done.")
    STOP_TIMING(matrix_gen);

    //Cholesky factorization for the Co-variance matrix CHAM_desc_K_a
    START_TIMING(cholesky1);
    VERBOSE("Cholesky factorization of CHAM_desc_K_a (MLOE-MMOM) .....")
    ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_desc_K_a, aConfigurations.GetBand(), nullptr, nullptr, 0, 0);
    VERBOSE("Done.")
    STOP_TIMING(cholesky1);
    flops = flops + flops_dpotrf(CHAM_desc_K_a->m);

    START_TIMING(cholesky2);
    //Cholesky factorization for the Co-variance matrix CHAM_desc_K_t
    VERBOSE("Cholesky factorization of CHAM_desc_K_t (MLOE-MMOM) .....")
    ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_desc_K_t, aConfigurations.GetBand(), nullptr, nullptr, 0, 0);
    VERBOSE("Done.")
    STOP_TIMING(cholesky2);
    flops = flops + flops_dpotrf(CHAM_desc_K_t->m);

    T total_loop_time = 0.0;
    T loop_time;
    bool verbose;
    VERBOSE("* Verbose messages are printed every 10 iteration *")
    for (p = 0; p < n_z_miss; p++) {
        verbose = p % 10 == 0;
        lmiss->GetLocationX()[0] = aMissLocations.GetLocationX()[p];
        lmiss->GetLocationY()[0] = aMissLocations.GetLocationY()[p];

        if (verbose) {
            VERBOSE("Generate two vectors k_a and k_t (MLOE-MMOM).....")
        }
        START_TIMING(vecs_gen);
        upper_lower = EXAGEOSTAT_UPPER_LOWER;
        this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_k_t, upper_lower, &aObsLocations, lmiss,
                                      &median_locations, apTruthTheta, 0, &aKernel);
        this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_k_a, upper_lower, &aObsLocations, lmiss,
                                      &median_locations, apEstimatedTheta, 0, &aKernel);
        this->ExaGeoStatSequenceWait(sequence);

        STOP_TIMING(vecs_gen);
        //Copy CHAM_desc_k_a to CHAM_descK_atmp  (MLOE-MMOM)
        if (verbose) {
            VERBOSE("Copy CHAM_desc_k_a to CHAM_descK_atmp  (MLOE-MMOM).....")
        }
        START_TIMING(copy_vecs);
        ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_k_t, CHAM_desc_k_t_tmp);
        ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_k_a, CHAM_desc_k_a_tmp);
        STOP_TIMING(copy_vecs);
        if (verbose) {
            VERBOSE("Done.")
        }

        START_TIMING(loop_time);
        START_TIMING(trsm1);
        // Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
        if (verbose) {
            VERBOSE("Solving the linear system k_a = TRSM(l_a^-1, k_a) ...(MLOE-MMOM)")
        }
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                           CHAM_desc_K_a, nullptr, nullptr, CHAM_desc_k_a, 0);
        if (verbose) {
            VERBOSE("Done.")
        }
        flops = flops + flops_dtrsm(ChamLeft, CHAM_desc_K_a->m, CHAM_desc_k_a->n);
        STOP_TIMING(trsm1);

        START_TIMING(trsm2);
        // Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
        if (verbose) {
            VERBOSE("Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)")
        }
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                           CHAM_desc_K_t, nullptr, nullptr, CHAM_desc_k_t, 0);
        flops = flops + flops_dtrsm(ChamLeft, CHAM_desc_K_t->m, CHAM_desc_k_t->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        STOP_TIMING(trsm2);

        START_TIMING(trsm3);
        // Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
        if (verbose) {
            VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)")
        }
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_K_a,
                           nullptr, nullptr, CHAM_desc_k_a, 0);
        flops = flops + flops_dtrsm(ChamLeft, CHAM_desc_K_a->m, CHAM_desc_k_a->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        STOP_TIMING(trsm3);

        START_TIMING(trsm4);
        // Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
        if (verbose) {
            VERBOSE("Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)")
        }
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_K_t,
                           nullptr, nullptr, CHAM_desc_k_t, 0);
        flops = flops + flops_dtrsm(ChamLeft, CHAM_desc_K_t->m, CHAM_desc_k_t->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        STOP_TIMING(trsm4);

        START_TIMING(gevv2);
        // Calculate dgemm value= CHAM_desc_k_t^T * CHAM_desc_k_a
        if (verbose) {
            VERBOSE("Calculate dgemm CHAM_desc_expr1 = CHAM_desc_k_t^T * CHAM_desc_k_a... (MLOE-MMOM)")
        }
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_k_t_tmp, CHAM_desc_k_a, 0, CHAM_desc_expr1);
        flops = flops + flops_dgemm(CHAM_desc_k_t_tmp->m, CHAM_desc_k_a->n, CHAM_desc_expr1->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        STOP_TIMING(gevv2);
        START_TIMING(gevv3);
        // Calculate dgemm value= CHAM_desc_k_a^T * CHAM_desc_k_a_tmp
        if (verbose) {
            VERBOSE("Calculate dgemm CHAM_desc_expr1 = CHAM_desc_k_a^T * CHAM_desc_k_a... (MLOE-MMOM)")
        }
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_k_a_tmp, CHAM_desc_k_a, 0, CHAM_desc_expr4);
        flops = flops + flops_dgemm(CHAM_desc_k_a_tmp->m, CHAM_desc_k_a->n, CHAM_desc_expr4->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        STOP_TIMING(gevv3);

        START_TIMING(gevv1);
        // Calculate dgemm value= CHAM_desc_k_a^T * CHAM_desc_k_t
        if (verbose) {
            VERBOSE("Calculate dgemm CHAM_desc_expr4 = CHAM_desc_k_a^T * CHAM_desc_k_t... (Prediction Stage)")
        }
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_k_t_tmp, CHAM_desc_k_t, 0, CHAM_desc_expr3);
        flops = flops + flops_dgemm(CHAM_desc_k_t_tmp->m, CHAM_desc_k_t->n, CHAM_desc_expr3->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        STOP_TIMING(gevv1);

        // Calculate dgemm CHAM_desc_k_a= CHAM_desc_K_t * CHAM_desc_k_a (use k_t as k_a)
        START_TIMING(gevv4);
        ExaGeoStatTrmmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_K_t,
                           CHAM_desc_k_a);
        STOP_TIMING(gevv4);

        // Calculate dgemm value= CHAM_desc_k_a^T * CHAM_desc_k_t
        if (verbose) {
            VERBOSE("Calculate dgemm CHAM_desc_expr1 = CHAM_desc_k_a^T * CHAM_desc_k_a... (Prediction Stage)")
        }
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_k_a, CHAM_desc_k_a, 0, CHAM_desc_expr2);
        flops = flops + flops_dgemm(CHAM_desc_k_a_tmp->m, CHAM_desc_k_t->n, CHAM_desc_expr2->n);
        if (verbose) {
            VERBOSE("Done.")
        }
        START_TIMING(gevv5);
        STOP_TIMING(gevv5);

        STOP_TIMING(loop_time);
        total_loop_time += loop_time;

        ExaGeoStatGeaddTile(EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_truth_alpha, -2, CHAM_desc_expr1);
        ExaGeoStatGeaddTile(EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_expr1, 1, CHAM_desc_expr2);
        ExaGeoStatGeaddTile(EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_truth_alpha, -1, CHAM_desc_expr3);
        ExaGeoStatGeaddTile(EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_estimated_alpha, -1, CHAM_desc_expr4);

        VERBOSE("- Matrix Generation Time: " << matrix_gen << " Vectors Generation Time: " << vecs_gen
                                             << " First Cholesky factorization Time: " << cholesky1
                                             << " First Cholesky factorization Time: " << cholesky2)
        VERBOSE("- First Trsm time: " << trsm1 << " Second Trsm time: " << trsm2 << " Third Trsm time: " << trsm3
                                      << " Fourth Trsm time: " << trsm4)
        VERBOSE("- First gemm time: " << gevv1 << " Second gemm time: " << gevv2 << " Third gemm time: " << gevv3
                                      << " Fourth gemm time: " << gevv4 << " Fifth gemm time: " << gevv5)
        ExaGeoStatMLETileAsyncMLOEMMOM(CHAM_desc_expr2, CHAM_desc_expr3, CHAM_desc_expr4, CHAM_desc_mloe,
                                       CHAM_desc_mmom, sequence, request);
        this->ExaGeoStatSequenceWait(sequence);
    }
    LOGGER(" ---- MLOE-MMOM Gflop/s: " << flops / 1e9 / (total_loop_time + cholesky1 + cholesky2))

    *mloe /= n_z_miss;
    *mmom /= n_z_miss;
    STOP_TIMING(all_time);
    LOGGER(" ---- MLOE = " << *mloe)
    LOGGER(" ---- MMOM = " << *mmom)
    LOGGER(" ---- MLOE MMOM Time: " << all_time << " seconds.")

    results::Results::GetInstance()->SetMLOE(*mloe);
    results::Results::GetInstance()->SetMMOM(*mmom);
    results::Results::GetInstance()->SetExecutionTimeMLOEMMOM(all_time);
    results::Results::GetInstance()->SetMatrixGenerationTimeMLOEMMOM(matrix_gen);
    results::Results::GetInstance()->SetFactoTimeMLOEMMOM(cholesky1 + cholesky2);
    results::Results::GetInstance()->SetLoopTimeMLOEMMOM(total_loop_time);
    results::Results::GetInstance()->SetFlopsMLOEMMOM((flops / 1e9 / (total_loop_time + cholesky1 + cholesky2)));

    delete[] loe;
    delete[] mom;
    delete[] temp1;
    delete[] temp2;
    delete[] temp3;
    delete[] estimated_alpha;
    delete[] truth_alpha;
    delete lmiss;
}

template<typename T>
T *LinearAlgebraMethods<T>::ExaGeoStatFisherTile(Configurations &aConfigurations,
                                                 std::unique_ptr<ExaGeoStatData<T>> &aData,
                                                 const ExaGeoStatHardware &aHardware, T *apTheta,
                                                 const kernels::Kernel<T> &aKernel) {

    this->SetContext(aHardware.GetChameleonContext());
    this->InitiateFisherDescriptors(aConfigurations, *aData->GetDescriptorData());

    auto *CHAM_desc_A = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_A).chameleon_desc;
    auto *CHAM_desc_C = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_C).chameleon_desc;
    auto *CHAM_desc_CJ = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_CJ).chameleon_desc;
    auto *CHAM_desc_results = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_RESULTS).chameleon_desc;
    auto *CHAM_desc_CK = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_CK).chameleon_desc;
    auto *CHAM_desc_C_diag = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C_DIAG).chameleon_desc;
    auto *CHAM_desc_C_trace = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_C_TRACE).chameleon_desc;

    auto trace = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_C_trace);
    *trace = 0.0;

    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    RUNTIME_sequence_t *sequence;
    if (!aData->GetDescriptorData()->GetSequence()) {
        ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    } else {
        sequence = (RUNTIME_sequence_t *) aData->GetDescriptorData()->GetSequence();
    }

    void *request = aData->GetDescriptorData()->GetRequest();

    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = aKernel.GetParametersNumbers();
    auto median_locations = Locations<T>(1, aData->GetLocations()->GetDimension());
    aData->CalculateMedianLocations(kernel_name, median_locations);
    double time = 0.0;

    START_TIMING(time);
    VERBOSE("Generate covariance matrix  CHAM_desc_C  (Fisher Matrix Generation).....")
    this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_C, EXAGEOSTAT_LOWER, aData->GetLocations(),
                                  aData->GetLocations(), &median_locations, apTheta,
                                  aConfigurations.GetDistanceMetric(),
                                  &aKernel);
    ExaGeoStatSequenceWait(sequence);
    VERBOSE("Done.")

    VERBOSE("Calculate Cholesky decomposition  (Fisher Matrix Generation).....")
    ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_desc_C, 0, nullptr, nullptr, 0, 0);
    VERBOSE("Done.")

    //Allocate memory for A, and initialize it with 0s.
    auto A = new T[num_params * num_params]();

    for (int j = 0; j < num_params; j++) {

        if (j == 0) {
            kernel_name = "UnivariateMaternDdsigmaSquare";
        } else if (j == 1) {
            kernel_name = "UnivariateMaternDbeta";
        } else if (j == 2) {
            kernel_name = "UnivariateMaternDnu";
        } else if (j == 3) {
            kernel_name = "UnivariateMaternNuggetsStationary";
        }

        auto pKernel_cj = plugins::PluginRegistry<kernels::Kernel<T>>::Create(kernel_name,
                                                                              aConfigurations.GetTimeSlot());

        VERBOSE("Generate covariance matrix  CHAM_desc_CJ  (Fisher Matrix Generation).....")
        this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_CJ, EXAGEOSTAT_UPPER_LOWER,
                                      aData->GetLocations(), aData->GetLocations(), &median_locations, apTheta,
                                      aConfigurations.GetDistanceMetric(), pKernel_cj);
        ExaGeoStatSequenceWait(sequence);
        VERBOSE("Done.")

        VERBOSE("Compute triangular solve  CHAM_desc_CJ  (Fisher Matrix Generation).....")
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                           CHAM_desc_C,
                           nullptr, nullptr, CHAM_desc_CJ, 0);
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                           CHAM_desc_C,
                           nullptr, nullptr, CHAM_desc_CJ, 0);
        VERBOSE("Done.")

        delete pKernel_cj;
        for (int k = j; k < num_params; k++) {

            if (k == 0) {
                kernel_name = "UnivariateMaternDdsigmaSquare";
            } else if (k == 1) {
                kernel_name = "UnivariateMaternDbeta";
            } else if (k == 2) {
                kernel_name = "UnivariateMaternDnu";
            } else if (k == 3) {
                kernel_name = "UnivariateMaternNuggetsStationary";
            }

            auto pKernel_ck = plugins::PluginRegistry<kernels::Kernel<T>>::Create(kernel_name,
                                                                                  aConfigurations.GetTimeSlot());

            VERBOSE("Generate covariance matrix  CHAM_desc_CK  (Fisher Matrix Generation).....")
            this->CovarianceMatrixCodelet(*aData->GetDescriptorData(), CHAM_desc_CK, EXAGEOSTAT_UPPER_LOWER,
                                          aData->GetLocations(), aData->GetLocations(), &median_locations, apTheta,
                                          aConfigurations.GetDistanceMetric(), pKernel_ck);
            ExaGeoStatSequenceWait(sequence);
            VERBOSE("Done.")

            VERBOSE("Compute tringular solve  CHAM_desc_CK  (Fisher Matrix Generation).....")
            ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                               CHAM_desc_C,
                               nullptr, nullptr, CHAM_desc_CK, 0);
            ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                               CHAM_desc_C,
                               nullptr, nullptr, CHAM_desc_CK, 0);
            VERBOSE("Done.")

            VERBOSE("Compute matrix-matrix multiplication  CHAM_desc_CK  (Fisher Matrix Generation).....")
            CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_desc_CJ, CHAM_desc_CK, 0, CHAM_desc_results);
            VERBOSE("Done.")

            VERBOSE("Compute the trace/diagonal of CHAM_desc_CK (Fisher Matrix Generation).....")
            ExaGeoStatMLETraceTileAsync(CHAM_desc_results, sequence, &request_array[0], CHAM_desc_C_trace,
                                        CHAM_desc_C_diag);
            ExaGeoStatSequenceWait(sequence);
            VERBOSE("Done.")

            A[k + num_params * j] = 0.5 * *trace;
            *trace = 0;
            delete pKernel_ck;
        }
    }

    STOP_TIMING(time);

    VERBOSE("Copy A array to descriptor CHAM_desc_A (Fisher Matrix Generation).....")
    ExaGeoStatLap2Desc(A, num_params, CHAM_desc_A, EXAGEOSTAT_UPPER_LOWER);
    VERBOSE("Done.")

    VERBOSE("Calculate Cholesky decomposition  (Fisher Matrix Generation).....")
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', num_params, (double *) A, num_params);
    VERBOSE("Done.")
    VERBOSE("Generate Identity Matrix (I) (Fisher Matrix Generation).....")

    //Allocate memory for A, and initialize it with 0s.
    auto I_matrix = new T[num_params * num_params + 1]();
    LAPACKE_dlaset(LAPACK_COL_MAJOR, 'L', num_params, num_params, 0, 1, (double *) I_matrix, num_params);
    VERBOSE("Done.")

    cblas_dtrsm(
            CblasColMajor,
            CblasLeft,
            CblasLower,
            CblasNoTrans,
            CblasNonUnit,
            num_params, num_params, 1.0, (double *) A, num_params, (double *) I_matrix, num_params);

    cblas_dtrsm(
            CblasColMajor,
            CblasLeft,
            CblasLower,
            CblasTrans,
            CblasNonUnit,
            num_params, num_params, 1.0, (double *) A, num_params, (double *) I_matrix, num_params);

    I_matrix[num_params * num_params] = time;

    results::Results::GetInstance()->SetTotalFisherTime(time);
    results::Results::GetInstance()->SetFisher00((T) I_matrix[0]);
    results::Results::GetInstance()->SetFisher11((T) I_matrix[4]);
    results::Results::GetInstance()->SetFisher22((T) I_matrix[8]);

    delete[] A;
    return I_matrix;
}

template<typename T>
void
LinearAlgebraMethods<T>::ExaGeoStatGetZObs(Configurations &aConfigurations, T *apZ, const int &aSize,
                                           DescriptorData<T> &aDescData, T *apMeasurementsMatrix, const int &aP) {

    auto z_desc = (CHAM_desc_t *) aDescData.GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
    if (!z_desc) {
        int full_problem_size = aConfigurations.GetProblemSize() * aP;
        int dts = aConfigurations.GetDenseTileSize();
        int p_grid = aConfigurations.GetPGrid();
        int q_grid = aConfigurations.GetQGrid();
        bool is_OOC = aConfigurations.GetIsOOC();


        // Set the floating point precision based on the template type
        FloatPoint float_point;
        if (sizeof(T) == SIZE_OF_FLOAT) {
            float_point = EXAGEOSTAT_REAL_FLOAT;
        } else if (sizeof(T) == SIZE_OF_DOUBLE) {
            float_point = EXAGEOSTAT_REAL_DOUBLE;
        } else {
            throw runtime_error("Unsupported for now!");
        }
        aDescData.SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, apMeasurementsMatrix, float_point, dts,
                                dts, dts * dts, full_problem_size, 1, 0, 0, full_problem_size, 1, p_grid, q_grid);
        z_desc = (CHAM_desc_t *) aDescData.GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY).chameleon_desc;
    }
    this->ExaGeoStatDesc2Lap(apZ, aSize, z_desc, UpperLower::EXAGEOSTAT_UPPER_LOWER);
}

template<typename T>
void LinearAlgebraMethods<T>::CovarianceMatrixCodelet(DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                                      const int &aTriangularPart, Locations<T> *apLocation1,
                                                      Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                      T *apLocalTheta, const int &aDistanceMetric,
                                                      const kernels::Kernel<T> *apKernel) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    RUNTIME_option_t options;
    RUNTIME_options_init((RUNTIME_option_t *) &options, (CHAM_context_t *) this->mpContext,
                         (RUNTIME_sequence_t *) aDescriptorData.GetSequence(),
                         (RUNTIME_request_t *) aDescriptorData.GetRequest());

    int temp, tempnn;
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
            temp = m == CHAM_apDescriptor->mt - 1 ? CHAM_apDescriptor->m - m * CHAM_apDescriptor->mb
                                                  : CHAM_apDescriptor->mb;
            m0 = m * CHAM_apDescriptor->mb;
            n0 = n * CHAM_apDescriptor->nb;
            // Register the data with StarPU
            starpu_insert_task(cl,
                               STARPU_VALUE, &temp, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocation2, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocation3, sizeof(Locations<T> *),
                               STARPU_VALUE, &apLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(kernels::Kernel<T> *),
                               0);
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
}

template<typename T>
void
LinearAlgebraMethods<T>::ExaGeoStatMLETraceTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescNum,
                                                     void *apDescTrace) {
    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m, m0, n0;
    int tempmm;
    auto A = *((CHAM_desc_t *) apDescA);
    struct starpu_codelet *cl = &this->cl_dtrace;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescA, m, m),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescNum, 0, 0),
                           STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescTrace, m, 0),
                           0);
    }

    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, this->mpContext);
}

template<typename T>
void
LinearAlgebraMethods<T>::ExaGeoStatGaussianToNonTileAsync(DescriptorData<T> &aDescriptorData, void *apDesc,
                                                          T *apTheta) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, aDescriptorData.GetSequence(), aDescriptorData.GetRequest());

    int m, m0;
    int temp;
    auto desc_z = (CHAM_desc_t *) apDesc;
    struct starpu_codelet *cl = &this->cl_gaussian_to_non;


    for (m = 0; m < desc_z->mt; m++) {
        temp = m == desc_z->mt - 1 ? desc_z->m - m * desc_z->mb : desc_z->mb;
        m0 = m * desc_z->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
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

    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, this->mpContext);
}

template<typename T>
int LinearAlgebraMethods<T>::ExaGeoStatMLEMSPETileAsync(void *apDescZPredict, void *apDescZMiss,
                                                        void *apDescError, void *apSequence,
                                                        void *apRequest) {
    // Check for Initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m, m0;
    int temp;
    auto Zpre = (CHAM_desc_t *) apDescZPredict;
    struct starpu_codelet *cl = &this->cl_dmse;

    for (m = 0; m < Zpre->mt; m++) {
        temp = m == Zpre->mt - 1 ? Zpre->m - m * Zpre->mb : Zpre->mb;

        m0 = m * Zpre->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescError, 0, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZPredict, m, 0),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescZMiss, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dmse",
#endif
                           0);
    }
    ExaGeoStatOptionsFree(&options);
    ExaGeoStatOptionsFinalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
int LinearAlgebraMethods<T>::ExaGeoStatMLETileAsyncMLOEMMOM(void *apDescExpr2, void *apDescExpr3, void *apDescExpr4,
                                                            void *apDescMLOE, void *apDescMMOM, void *apSequence,
                                                            void *apRequest) {

    // Check for Initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m, n, m0, n0;
    struct starpu_codelet *cl = &this->cl_dmloe_mmom;
    int temp, temp2;

    for (n = 0; n < ((CHAM_desc_t *) apDescExpr2)->nt; n++) {
        temp2 = n == ((CHAM_desc_t *) apDescExpr2)->nt - 1 ? ((CHAM_desc_t *) apDescExpr2)->n -
                                                             n * ((CHAM_desc_t *) apDescExpr2)->nb
                                                           : ((CHAM_desc_t *) apDescExpr2)->nb;
        for (m = 0; m < ((CHAM_desc_t *) apDescExpr2)->mt; m++) {

            temp = m == ((CHAM_desc_t *) apDescExpr2)->mt - 1 ? ((CHAM_desc_t *) apDescExpr2)->m -
                                                                m * ((CHAM_desc_t *) apDescExpr2)->mb
                                                              : ((CHAM_desc_t *) apDescExpr2)->mb;
            m0 = m * ((CHAM_desc_t *) apDescExpr2)->mb;
            n0 = n * ((CHAM_desc_t *) apDescExpr2)->nb;
            starpu_insert_task(cl,
                               STARPU_VALUE, &temp, sizeof(int),
                               STARPU_VALUE, &temp2, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescExpr2, m, n),
                               STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescExpr3, m, n),
                               STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescExpr4, m, n),
                               STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescMLOE, m, n),
                               STARPU_RW, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescMMOM, m, n),
                               0);
        }
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
void
LinearAlgebraMethods<T>::CopyDescriptorZ(DescriptorData<T> &aDescriptorData, void *apDescriptor, T *apDoubleVector) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (chameleon_context_s *) this->mpContext,
                         (RUNTIME_sequence_t *) aDescriptorData.GetSequence(),
                         (RUNTIME_request_t *) aDescriptorData.GetRequest());

    int m, m0;
    int temp;
    auto A = (CHAM_desc_t *) apDescriptor;
    struct starpu_codelet *cl = &this->cl_dzcpy;

    for (m = 0; m < A->mt; m++) {
        temp = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
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
int
LinearAlgebraMethods<T>::ExaGeoStatMLEMSPEBivariateTileAsync(void *apDescZPre, void *apDescZMiss, void *apDescsError1,
                                                             void *apDescsError2, void *apDescsError,
                                                             void *apSequence, void *apRequest) {
    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (chameleon_context_s *) this->mpContext,
                         (RUNTIME_sequence_t *) apSequence,
                         (RUNTIME_request_t *) apRequest);

    int m, m0;
    int tempmm;
    auto Zpre = (CHAM_desc_t *) apDescZPre;

    struct starpu_codelet *cl = &cl_dmse_bivariate;


    for (m = 0; m < Zpre->mt; m++) {
        tempmm = m == Zpre->mt - 1 ? Zpre->m - m * Zpre->mb : Zpre->mb;
        m0 = m * Zpre->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, ExaGeoStatDataGetAddr(apDescsError1, 0, 0),
                           STARPU_RW, ExaGeoStatDataGetAddr(apDescsError2, 0, 0),
                           STARPU_RW, ExaGeoStatDataGetAddr(apDescsError, 0, 0),
                           STARPU_R, ExaGeoStatDataGetAddr(apDescZPre, m, 0),
                           STARPU_R, ExaGeoStatDataGetAddr(apDescZMiss, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dmse_bivariate",
#endif
                           0);
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}

#ifdef USE_HICMA

template<typename T>
void LinearAlgebraMethods<T>::CopyDescriptors(void *apSourceDesc, void *apDestinationDesc, const int &aSize,
                                              const common::CopyDirection &aDirection) {
    auto *z = new T[aSize];
    int status;
    if (aDirection == common::CHAMELEON_TO_HICMA) {
        status = CHAMELEON_Desc2Lap((cham_uplo_t) EXAGEOSTAT_UPPER_LOWER, (CHAM_desc_t *) apSourceDesc, z, aSize);
        if (status != CHAMELEON_SUCCESS) {
            throw std::runtime_error("CHAMELEON_Desc2Lap Failed!");
        }
        status = HICMA_Lapack_to_Tile(z, aSize, (HICMA_desc_t *) apDestinationDesc);
        if (status != HICMA_SUCCESS) {
            throw std::runtime_error("HICMA_Lapack_to_Tile Failed!");
        }
    } else if (aDirection == common::HICMA_TO_CHAMELEON) {
        status = HICMA_Tile_to_Lapack((HICMA_desc_t *) apSourceDesc, z, aSize);
        if (status != HICMA_SUCCESS) {
            throw std::runtime_error("HICMA_Tile_to_Lapack Failed!");
        }
        status = CHAMELEON_Lap2Desc((cham_uplo_t) EXAGEOSTAT_UPPER_LOWER, z, aSize, (CHAM_desc_t *) apDestinationDesc);
        if (status != CHAMELEON_SUCCESS) {
            throw std::runtime_error("CHAMELEON_Lap2Desc Failed!");
        }
    }
    delete[] z;
}

#endif

template<typename T>
void
LinearAlgebraMethods<T>::ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA, const UpperLower &aUpperLower) {
    int status = CHAMELEON_Desc2Lap((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apDescA, apA, aLDA);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Desc2Lap Failed!");
    }
}


template<typename T>
void
LinearAlgebraMethods<T>::ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const UpperLower &aUpperLower) {
    int status = CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, apA, aLDA, (CHAM_desc_t *) apDescA);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Lap2Desc Failed!");
    }
}

template<typename T>
void LinearAlgebraMethods<T>::ExaGeoStatLaSetTile(const common::UpperLower &aUpperLower, T alpha, T beta,
                                                  void *apDescriptor) {
    int status = CHAMELEON_dlaset_Tile((cham_uplo_t) aUpperLower, alpha, beta, (CHAM_desc_t *) apDescriptor);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dlaset_Tile Failed!");
    }
}

template<typename T>
void LinearAlgebraMethods<T>::ExaGeoStatTrmmTile(const Side &aSide, const UpperLower &aUpperLower, const Trans &aTrans,
                                                 const Diag &aDiag, const T &alpha, void *apDescA, void *apDescB) {

    int status = CHAMELEON_dtrmm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans,
                                      (cham_diag_t) aDiag, alpha, (CHAM_desc_t *) apDescA, (CHAM_desc_t *) apDescB);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dtrmm_Tile Failed!");
    }
}

template<typename T>
void LinearAlgebraMethods<T>::ExaGeoStatGeaddTile(const common::Trans &aTrans, const T &aAlpha, void *apDescA,
                                                  const T &aBeta, void *apDescB) {
    int status = CHAMELEON_dgeadd_Tile((cham_trans_t) aTrans, aAlpha, (CHAM_desc_t *) apDescA, aBeta,
                                       (CHAM_desc_t *) apDescB);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dgeadd_Tile Failed!");
    }
}

template<typename T>
void LinearAlgebraMethods<T>::ExaGeoStatPosvTile(const common::UpperLower &aUpperLower, void *apA, void *apB) {
    int status = CHAMELEON_dposv_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dposv_Tile Failed!");
    }
}