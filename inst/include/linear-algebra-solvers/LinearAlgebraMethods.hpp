
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraMethods.hpp
 * @brief Header file for the LinearAlgebraMethods class, which defines the interface for linear algebra solvers.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
 * @details This header file defines the abstract class LinearAlgebraMethods, which provides an interface for linear algebra solvers.
 * The purpose of this interface is to allow different concrete linear algebra solvers to be interchangeable,
 * so that they can be used interchangeably by other parts of the software system that rely on linear algebra.
 *
**/

#ifndef EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP
#define EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP

#include <vector>

extern "C" {
#include <gsl/gsl_errno.h>
}

#include <common/Utils.hpp>
#include <helpers/DiskWriter.hpp>
#include <kernels/Kernel.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <results/Results.hpp>

namespace exageostat::linearAlgebra {

    /**
     * @class LinearAlgebraMethods
     * @brief A class that defines the interface for linear algebra solvers.
     * @tparam T Data Type: float or double.
     *
     */
    template<typename T>
    class LinearAlgebraMethods {
    public:

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        virtual ~LinearAlgebraMethods() = default;

        /**
         * @brief Initializes the descriptors necessary for the linear algebra solver.
         * @details This method initializes the descriptors necessary for the linear algebra solver.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aDescriptorData Descriptor Data object to be populated with descriptors and data.
         * @param[in] apMeasurementsMatrix Pointer to the measurement matrix.
         * @return void
         *
         */
        void InitiateDescriptors(configurations::Configurations &aConfigurations,
                                 dataunits::DescriptorData<T> &aDescriptorData, T *apMeasurementsMatrix = nullptr);

        /**
         * @brief Initializes the descriptors necessary for the Fisher prediction function.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aDescriptorData Descriptor Data object to be populated with descriptors and data.
         */
        void InitiateFisherDescriptors(configurations::Configurations &aConfigurations,
                                       dataunits::DescriptorData<T> &aDescriptorData);

            /**
         * @brief Initializes the descriptors necessary for the Prediction.
         * @details This method initializes the descriptors necessary for the linear algebra solver.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData DescriptorData object to be populated with descriptors and data.
         * @return void
         *
         */
        void InitiatePredictionDescriptors(configurations::Configurations &aConfigurations,
                                           dataunits::ExaGeoStatData<T> &aData);

        /**
         * @brief Initializes the descriptors necessary for the Prediction Auxiliary function MLE-MLOE-MMOM.
         * @details This method initializes the descriptors necessary for the linear algebra solver.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData DescriptorData object to be populated with descriptors and data.
         * @return void
         *
         */
        void InitiateMLOEMMOMDescriptors(configurations::Configurations &aConfigurations,
                                         dataunits::ExaGeoStatData<T> &aData);

        /**
         * @brief Generates synthetic data.
         * @param[in] aConfigurations The configurations object containing relevant settings.
         * @param[in] aHardware  ExaGeoStatHardware object representing the hardware.
         * @param[in,out] aData ExaGeoStatData object to be populated with synthetic data.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return None.
         *
         */
        void GenerateSyntheticData(configurations::Configurations &aConfigurations,
                                   const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &aData,
                                   const kernels::Kernel<T> &aKernel);

        /**
         * @brief Computes the covariance matrix.
         * @param[in] aDescriptorData pointer to the DescriptorData object holding descriptors and data.
         * @param[out] apDescriptor Pointer to the descriptor for the covariance matrix.
         * @param[in] aTriangularPart Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @param[in] apLocation1 Pointer to the first set of locations.
         * @param[in] apLocation2 Pointer to the second set of locations.
         * @param[in] apLocation3 Pointer to the third set of locations.
         * @param[in] apLocalTheta Pointer to the local theta values.
         * @param[in] aDistanceMetric Specifies the distance metric to use.
         * @param[in] apKernel Pointer to the kernel object to use.
         * @return void
         *
         */
        void CovarianceMatrixCodelet(dataunits::DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                     const int &aTriangularPart, dataunits::Locations<T> *apLocation1,
                                     dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                     T *apLocalTheta, const int &aDistanceMetric, const kernels::Kernel<T> *apKernel);

        /**
         * @brief Copies the descriptor data to a double vector.
         * @param[in] aDescriptorData pointer to the DescriptorData object holding descriptors and data.
         * @param[in] apDescriptor Pointer to the descriptor data.
         * @param[in,out] apDoubleVector Pointer to the double vector to copy the descriptor data to.
         * @return void
         *
         */
        void
        CopyDescriptorZ(dataunits::DescriptorData<T> &aDescriptorData, void *apDescriptor, T *apDoubleVector);

        /**
         * @brief Generates the observations vector.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] aDescriptorData pointer to the DescriptorData object holding descriptors and data.
         * @param[in] apLocation1 Pointer to the first set of locations.
         * @param[in] apLocation2 Pointer to the second set of locations.
         * @param[in] apLocation3 Pointer to the third set of locations.
         * @param[in] aDistanceMetric Specifies the distance metric to use.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return void
         *
         */
        void
        GenerateObservationsVector(configurations::Configurations &aConfigurations, dataunits::ExaGeoStatData<T> &aData,
                                   dataunits::Locations<T> *apLocation1, dataunits::Locations<T> *apLocation2,
                                   dataunits::Locations<T> *apLocation3, const int &aDistanceMetric,
                                   const kernels::Kernel<T> &aKernel);

        /**
         * @brief Calculates the log likelihood value of a given value theta.
         * @param[in] aHardware  ExaGeoStatHardware object representing the hardware.
         * @param[in,out] aData DescriptorData object to be populated with descriptors and data.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] apTheta Optimization parameter used by NLOPT.
         * @param[in] apMeasurementsMatrix measurements matrix to be stored in DescZ.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return log likelihood value
         *
         */
        virtual T ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &aData,
                                    configurations::Configurations &aConfigurations, const double *apTheta,
                                    T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) = 0;

        /**
         * @brief Converts a Gaussian descriptor to a non-tiled descriptor.
         * @param[in] aDescriptorData DescriptorData struct with the Gaussian descriptor.
         * @param[in] apDesc Pointer to the non-tiled descriptor.
         * @param[in] apTheta Theta vector.
         * @return void
         */
        void ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData<T> &aDescriptorData, void *apDesc, T *apTheta);

        /**
         * @brief Copies a matrix in the tile layout from source to destination
         * @param[in] aUpperLower Specifies the part of the matrix A to be copied to B.
         * @param[in] apA Source matrix A.
         * @param[in,out] apB Destination matrix B. On exit, B = A in the locations specified by Upper Lower.
         * @return void
         *
         */
        virtual void ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) = 0;

        /**
         * @brief Initialize the runtime option structure for either HiCMA or CHAMELEON.
         * @param[in, out] apOptions The options structure that needs to be initialized.
         * @param[in] apContext The runtime context in which to initialize the runtime support.
         * @param[in] apSequence The sequence structure to associate in the options.
         * @param[in] apRequest The request structure to associate in the options.
         * @return void
         */
        virtual void
        ExaGeoStatOptionsInit(void *apOptions, void *apContext, void *apSequence, void *apRequest) = 0;

        /**
         * @brief Submit the release of the workspaces associated to the options structure.
         * @param[in,out] apOptions The options structure for which to workspaces will be released
         * @return void
         */
        virtual void ExaGeoStatOptionsFree(void *apOptions) = 0;

        /**
         * @brief Finalize the runtime option structure for either HiCMA or CHAMELEON.
         * @param[in,out] apOptions The options structure that needs to be finalized.
         * @param[in] apContext The runtime context in which to finalize the runtime support.
         * @return void
         *
         */
        virtual void ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) = 0;

        /**
         * @brief Wait for the completion of a sequence.
         * @param[in] apSequence apSequence A pointer to either CHAMELEON or HiCMA sequence.
         * @return void
         *
         */
        virtual void ExaGeoStatSequenceWait(void *apSequence) = 0;

        /**
         * @brief Create Sequence.
         * @param[out] apSequence A pointer to either CHAMELEON or HiCMA sequence.
         * @return void
         *
         */
        virtual void
        ExaGeoStatCreateSequence(void *apSequence) = 0;

        /**
        * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
        * @param[in] aUpperLower Whether upper or lower part of the matrix A.
        * @param[in, out] apA Symmetric matrix A.
        * @param[in] aBand Diagonal thickness parameter.
        * @param[in] apCD Additional matrix CD.
        * @param[in] apCrk Additional matrix Crk.
        * @param[in] aMaxRank Maximum rank parameter.
        * @param[in] aAcc Accuracy parameter.
        * @return void
        */
        virtual void
        ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apCD, void *apCrk,
                            const int &aMaxRank, const int &aAcc) = 0;

        /**
         * @brief Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
         * @param[in] aSide Specifies whether op(A) appears on the left or on the right of X.
         * @param[in] aUpperLower Specifies whether the matrix A is upper triangular or lower triangular.
         * @param[in] aTrans Specifies the form of op( A ) to be used in the matrix multiplication.
         * @param[in] aDiag Specifies whether or not A is unit triangular.
         * @param[in] aAlpha Specifies the scalar alpha. When alpha is zero, A is not referenced and B need not be set before entry.
         * @param[in] apA The triangular matrix A.
         * @param[in] apCD Additional matrix CD.
         * @param[in] apCrk Additional matrix Crk.
         * @param[in, out] apZ The matrix B of dimension, on exit is overwritten by the solution matrix X.
         * @param[in] aMaxRank Maximum rank parameter.
         * @return void
         */
        virtual void ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                        const common::Trans &aTrans, const common::Diag &aDiag, const T &aAlpha,
                                        void *apA, void *apCD, void *apCrk, void *apZ, const int &aMaxRank) = 0;

        /**
         * @brief Calculate determinant for triangular matrix.
         * @param[in] apDescA Pointer to the descriptor of the matrix 'descA'.
         * @param[in] apSequence Pointer to a sequence structure for managing asynchronous execution.
         * @param[in] apRequest Pointer to a request structure for tracking the operation's status.
         * @param[out] apDescNum Pointer to the descriptor of the matrix to store the sum of elements.
         * @param[out] apDescTrace Pointer to the descriptor of the matrix to store the trace.
         * @return void
         */
        void
        ExaGeoStatMLETraceTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescNum,
                                    void *apDescTrace);

            /**
         * @brief Calculate determinant for triangular matrix.
         * @param[in] apDescA Exageostat descriptor.
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[in] apRequest Identifies this function call (for exception handling purposes).
         * @param[in] apDescDet determinant value
         * @return Returns 0 for success, error code otherwise.
         *
         */
        virtual int
        ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescDet) = 0;
        /**
         * @brief Solve a positive definite linear system of equations AX = B using tiled algorithms.
         * @param[in] aUpperLower Specifies whether the matrix A is upper triangular or lower triangular.
         * @param [in] apA coefficient matrix of the system of linear equations. This matrix is expected to be positive definite.
         * @param [in] apB Pointer to coefficient matrix of the system of linear equations. This matrix is expected to be positive definite.
         * @return void
         */
        void ExaGeoStatPosvTile(const common::UpperLower &aUpperLower, void *apA, void *apB);

        /**
         * @brief Calculate mean square prediction error (MSPE) scalar value of the prediction.
         * @param[in] apDescZPredict Observed measurements.
         * @param[in] apDescZMiss Missing measurements.
         * @param[out] apDescError Mean Square Prediction Error (MSPE).
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[out] apRequest Identifies this function call (for exception handling purposes).
         * @return Returns 0 for success, error code otherwise.
         */
        int ExaGeoStatMLEMSPETileAsync(void *apDescZPredict, void *apDescZMiss, void *apDescError, void *apSequence,
                                      void *apRequest);

        /**
         * Predict missing values base on a set of given values and covariance matrix/
         * @param[in] aData Reference to Data containing different MLE inputs.
         * @param[in] apTheta theta Vector with three parameter (Variance, Range, Smoothness) that is used to to generate the Covariance Matrix.
         * @param[in] aZMissNumber number of missing values (unknown observations).
         * @param[in] aZObsNumber number of observed values (known observations).
         * @param[in] apZObs observed values vector (known observations).
         * @param[in] apZActual actual missing values vector (in the case of testing MSPE).
         * @param[in] apZMiss missing values vector (unknown observations).
         * @param[in] aHardware  ExaGeoStatHardware object representing the hardware.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] aMissLocations Reference to Locations object containing missed locations.
         * @param[in] aObsLocations Reference to Locations object containing observed locations.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return the prediction Mean Square Error (MSPE).
         */
        T * ExaGeoStatMLEPredictTile(exageostat::dataunits::ExaGeoStatData<T> &aData, T *apTheta, const int &aZMissNumber,
                                     const int &aZObsNumber, T *apZObs, T *apZActual, T *apZMiss,
                                     const hardware::ExaGeoStatHardware &aHardware,
                                     configurations::Configurations &aConfiguration,
                                     exageostat::dataunits::Locations<T> &aMissLocations,
                                     exageostat::dataunits::Locations<T> &aObsLocations, const kernels::Kernel<T> &aKernel);

        /**
         * @brief Copy Lapack matrix to Descriptor Matrix.
         * @param[in] apA Lapack Matrix.
         * @param[in] aLDA Size.
         * @param[out] apDescA Matrix Descriptor.
         * @param[in] aUpperLower Specifies Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @return void
         */
        virtual void
        ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const common::UpperLower &aUpperLower) = 0;

        /**
         * @brief Copy Descriptor Matrix to Lapack matrix.
         * @param[out] apA Lapack Matrix.
         * @param[in] aLDA Size.
         * @param[in] apDescA Matrix Descriptor
         * @param[in] aUpperLower Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @return void
         */
        void ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA, const common::UpperLower &aUpperLower);

        /**
         * @brief Sets the values of all or part of a two-dimensional Tile.
         * @param[in] aUpperLower Specifies Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @param[in] aAlpha All the off diagonal array elements are set to aAlpha.
         * @param[in] aBeta All the diagonal array elements are set to aBeta.
         * @param[out] apDescriptor Pointer to matrix descriptor to be set with aAlpha and aBeta.
         * @return void
         */
        void ExaGeoStatLaSetTile(const common::UpperLower &aUpperLower, T aAlpha, T aBeta, void *apDescriptor);

        /**
         * @brief Copy the Z matrix into a pointer.
         * @param[out] apZ Pointer to an array to copy Z matrix into.
         * @param[in] aSize Size of the matrix.
         * @param[in] aDescData Descriptor data containing required Z matrix Descriptor.
         */
        void ExaGeoStatGetZObs(exageostat::configurations::Configurations &aConfigurations, T *apZ, const int &aSize,
                               exageostat::dataunits::DescriptorData<T> &aDescData, T *apMeasurementsMatrix);

        /**
         * @brief Predict missing values based on a set of given values and covariance matrix.
         * @details This function predicts missing values using the maximum likelihood estimation (MLE),
         * maximum likelihood on the empirical orthogonal functions (MLOE), and method of moments (MMOM).
         * @param[in] aConfigurations Configurations for the prediction.
         * @param[in, out] aData Data for prediction (input and output).
         * @param[in] aHardware Hardware specifications for the prediction.
         * @param[in] apTruthTheta Pointer to the true theta values.
         * @param[in] apEstimatedTheta Pointer to the estimated theta values.
         * @param[in] aMissLocations Locations of missing values.
         * @param[in] aObsLocations Locations of observed values.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return void
         */
        void ExaGeoStatMLETileMLOEMMOM(exageostat::configurations::Configurations &aConfigurations,
                                       exageostat::dataunits::ExaGeoStatData<T> &aData,
                                       const exageostat::hardware::ExaGeoStatHardware &aHardware, T *apTruthTheta,
                                       T *apEstimatedTheta, dataunits::Locations<T> &aMissLocations,
                                       dataunits::Locations<T> &aObsLocations, const kernels::Kernel<T> &aKernel);

        /**
         * @brief Maximum Likelihood Evaluation (MLE) Fisher method.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         * @param[in] aHardware  ExaGeoStatHardware object representing the hardware.
         * @param[in] apTheta Pointer containing three parameter (Variance, Range, Smoothness) that is used to to generate the Covariance Matrix.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return Fisher Matrix
         */
        T *
        ExaGeoStatFisherTile(configurations::Configurations &aConfigurations, dataunits::ExaGeoStatData<T> &aData,
                             const hardware::ExaGeoStatHardware &aHardware, T *apTheta,
                             const kernels::Kernel <T> &aKernel);

        /**
         * @brief Perform an asynchronous computation of MLE, MLOE, and MMOM for a tile.
         * @details his function performs the computation of Maximum Likelihood Estimation (MLE),
         * Maximum Likelihood on the Empirical Orthogonal Functions (MLOE), and
         * Method of Moments (MMOM) for a tile asynchronously.
         * @param[in] apDescExpr2 Descriptor for expression 2.
         * @param[in] apDescExpr3 Descriptor for expression 3.
         * @param[in] apDescExpr4 Descriptor for expression 4.
         * @param[in] apDescMLOE Descriptor for MLOE.
         * @param[in] apDescMMOM Descriptor for MMOM.
         * @param[in] apSequence Sequence for the computation.
         * @param[in] apRequest Request for the computation.
         * @return Result of the asynchronous operation.
         */
        int ExaGeoStatMLETileAsyncMLOEMMOM(void *apDescExpr2, void *apDescExpr3, void *apDescExpr4, void *apDescMLOE,
                                           void *apDescMMOM, void *apSequence, void *apRequest);

        /**
         * @brief Perform a matrix addition with scaling.
         * @details This function performs a matrix addition with scaling, given the matrices A and B.
         * @param[in] aTrans Specifies whether to transpose matrix A.
         * @param[in] aAlpha Scaling factor for matrix A.
         * @param[in] apDescA Descriptor for matrix A.
         * @param[in] aBeta Scaling factor for matrix B.
         * @param[in] apDescB Descriptor for matrix B.
         * @return void
         */
        void ExaGeoStatGeaddTile(const common::Trans &aTrans, const T &aAlpha, void *apDescA, const T &aBeta,
                                 void *apDescB);

        /**
         * @brief Perform a triangular matrix multiplication.
         * @details This function performs triangular matrix multiplication on two matrices A and B.
         * @param[in] aSide Specifies whether the multiplication is performed on the left or right side.
         * @param[in] aUpperLower Specifies whether the matrix is upper or lower triangular.
         * @param[in] aTrans Specifies whether to transpose the matrix.
         * @param[in] aDiag Specifies whether the diagonal elements are unitary or non-unitary.
         * @param[in] alpha Scaling factor for the multiplication.
         * @param[in] apDescA Descriptor for matrix A.
         * @param[in] apDescB Descriptor for matrix B.
         */
        void ExaGeoStatTrmmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                const common::Trans &aTrans, const common::Diag &aDiag, const T &alpha, void *apDescA,
                                void *apDescB);

        /**
         * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
         * @param[in] apA The descriptor to which belongs the piece of data
         * @param[in] aAm The row coordinate of the piece of data in the matrix
         * @param[in] aAn he column coordinate of the piece of data in the matrix
         * @return The runtime handler address of the piece of data.
         */
        virtual void *ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) = 0;

        /**
         * @brief Sets the context.
         * @param[in] apContext The context.
         *
         */
        void SetContext(void *apContext) {
            this->mpContext = apContext;
        }

        //// TODO: Create a Factory for Runtime system. This is a solution fix for now
        struct starpu_codelet cl_gaussian_to_non =
                {
                        .where        = STARPU_CPU,
                        .cpu_funcs    = {CORE_gaussian_to_non_starpu},
                        .nbuffers    = 1,
                        .modes        = {STARPU_RW},
                        .name        = "gaussian_to_non"
                };

        struct starpu_codelet cl_dzcpy =
                {
                        .where        = STARPU_CPU,
                        .cpu_funcs    = {CORE_dzcpy_Starpu},
                        .nbuffers    = 1,
                        .modes        = {STARPU_W},
                        .name        = "dzcpy"
                };

        struct starpu_codelet cl_dmdet =
                {
                        .where        = STARPU_CPU,
                        .cpu_funcs    = {CORE_dmdet_starpu},
                        .nbuffers    = 2,
                        .modes        = {STARPU_R, STARPU_RW},
                        .name        = "dmdet"
                };

        struct starpu_codelet cl_stride_vec =
                {
                        .where          = STARPU_CPU,
                        .cpu_funcs      = {CORE_stride_vector_starpu},
                        .nbuffers       = 3,
                        .modes          = {STARPU_R, STARPU_W, STARPU_W},
                        .name           = "stride_vec"
                };

        struct starpu_codelet cl_tri_stride_vec =
                {
                        .where          = STARPU_CPU,
                        .cpu_funcs      = {CORE_tri_stride_vector_starpu},
                        .nbuffers       = 4,
                        .modes          = {STARPU_R, STARPU_W, STARPU_W, STARPU_W},
                        .name           = "tri_stride_vec"
                };

        struct starpu_codelet cl_dmse =
                {
                        .where        = STARPU_CPU,
                        .cpu_funcs    = {CORE_dmse_starpu},
                        .nbuffers    = 3,
                        .modes        = {STARPU_RW, STARPU_R, STARPU_R},
                        .name        = "dmse"
                };

        struct starpu_codelet cl_dmloe_mmom =
                {
                        .where        = STARPU_CPU,
                        .cpu_funcs    = {CORE_dmloe_mmom_starpu},
                        .nbuffers    = 5,
                        .modes        = {STARPU_R, STARPU_R, STARPU_R, STARPU_RW, STARPU_RW},
                        .name        = "dmloe_mmom"
                };

        struct starpu_codelet cl_dcmg =
                {
                        .where        = STARPU_CPU,
                        .cpu_funcs     = {CORE_dcmg_starpu},
                        .nbuffers     = 1,
                        .modes        = {STARPU_W},
                        .name         = "dcmg"
                };

        struct starpu_codelet cl_dtrace =
                {
                        .where		= STARPU_CPU,
                        .cpu_funcs	= {CORE_dtrace_starpu},
                        .nbuffers	= 3,
                        .modes		= {STARPU_R, STARPU_RW, STARPU_W},
                        .name		= "dtrace"
                };

        struct starpu_codelet cl_ddotp =
                {
                        .where		= STARPU_CPU,
                        .cpu_funcs	= {CORE_ddotp_starpu},
                        .nbuffers	= 2,
                        .modes		= {STARPU_RW,STARPU_R},
                        .name		= "ddotp"
                };



        static void CORE_dcmg_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m, n, m0, n0;
            exageostat::dataunits::Locations<T> *location1;
            exageostat::dataunits::Locations<T> *location2;
            exageostat::dataunits::Locations<T> *location3;
            T *theta;
            T *A;
            int distance_metric;
            kernels::Kernel<T> *kernel;

            A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            starpu_codelet_unpack_args(apCodeletArguments, &m, &n, &m0, &n0, &location1, &location2, &location3, &theta,
                                       &distance_metric, &kernel);
            kernel->GenerateCovarianceMatrix(A, m, n, m0, n0, *location1, *location2, *location3, theta,
                                             distance_metric);
        }

        static void CORE_gaussian_to_non_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m, m0;
            T *z;
            T *theta;
            theta = new T[6];
            z = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);

            starpu_codelet_unpack_args(apCodeletArguments, &m, &m0, &theta[0], &theta[1], &theta[2], &theta[3],
                                       &theta[4], &theta[5]);
            //core function to convert Z tile from Gaussian to non-Gaussian.
            core_gaussian_to_non(z, theta, m);
            delete[] theta;
        }

        static void core_gaussian_to_non(T *Z, T *apLocalTheta, int aSize) {

            double xi = apLocalTheta[2];
            double omega = apLocalTheta[3];
            double g = apLocalTheta[4];
            double h = apLocalTheta[5];

            int i;
            if (h < 0) {
                LOGGER("The kurtosis parameter cannot be negative")
                return;
            }
            if (g == 0) {
                for (i = 0; i < aSize; i++)
                    Z[i] = xi + omega * Z[i] * (exp(0.5 * h * pow(Z[i], 2)));
            } else {
                for (i = 0; i < aSize; i++)
                    Z[i] = xi + omega * (exp(g * Z[i]) - 1) * (exp(0.5 * h * pow(Z[i], 2))) / g;
            }
        }

        static void CORE_dzcpy_Starpu(void *apBuffers[], void *apCodeletArguments) {
            int m;
            T *A;
            int m0;
            T *r;

            A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            starpu_codelet_unpack_args(apCodeletArguments, &m, &m0, &r);
            memcpy(A, &r[m0], m * sizeof(T));
        }

        static void CORE_dmdet_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m;
            int n;
            T *A;
            int m0;
            int n0;
            T det = 0;
            T *determinant = &det;

            *determinant = 0;
            A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            determinant = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
            starpu_codelet_unpack_args(apCodeletArguments, &m, &n, &m0, &n0);
            T local_det = Core_dmdet(A, m);
            *determinant += local_det;
        }

        static T Core_dmdet(T *apA, int aSize) {

            int i;
            T res = 0.0;
            for (i = 0; i < aSize; i++) {
                if (apA[i + i * aSize] > 0)
                    res += log(apA[i + i * aSize]);
            }
            return res;
        }

        static void CORE_stride_vector_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m;
            int temp;
            T *A;
            T *B;
            T *C;
            int m0;
            int i;
            int j = 0;

            A = (T *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            B = (T *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
            C = (T *) STARPU_MATRIX_GET_PTR(apBuffers[2]);
            starpu_codelet_unpack_args(apCodeletArguments, &temp, &m0, &m);

            for (i = 0; i < temp - 1; i += 2) {
                B[j] = A[i];
                C[j] = A[i + 1];
                j++;
            }
        }

        static void CORE_tri_stride_vector_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m;
            int temp;
            double *A;
            double *B;
            double *C;
            double *D;
            int m0;
            int i;
            int j;

            A = (double *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            B = (double *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
            C = (double *) STARPU_MATRIX_GET_PTR(apBuffers[2]);
            D = (double *) STARPU_MATRIX_GET_PTR(apBuffers[3]);

            starpu_codelet_unpack_args(apCodeletArguments, &temp, &m0, &m);

            //accept only temp divided by three (should be optimized)
            j = 0;
            for (i = 0; i < temp - 1; i += 3) {
                B[j] = A[i];
                C[j] = A[i + 1];
                D[j] = A[i + 2];
                j++;
            }
        }

        static void CORE_dmse_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m, m0, i;
            double *pZPredict;
            double *pZMiss;
            double *pError;
            double local_error = 0.0;

            pError = (double *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            pZPredict = (double *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
            pZMiss = (double *) STARPU_MATRIX_GET_PTR(apBuffers[2]);

            starpu_codelet_unpack_args(apCodeletArguments, &m, &m0);
            for (i = 0; i < m; i++) {
                local_error += pow((pZPredict[i] - pZMiss[i]), 2);
            }
            *pError += local_error;
        }

        static void CORE_dmloe_mmom_starpu(void *apBuffers[], void *apCodeletArguments) {
            int m;
            int n;
            int i;
            double *expr2;
            double *expr3;
            double *expr4;
            double *mloe;
            double *mmom;
            int m0;
            int n0;

            expr2 = (double *) STARPU_MATRIX_GET_PTR(apBuffers[0]);
            expr3 = (double *) STARPU_MATRIX_GET_PTR(apBuffers[1]);
            expr4 = (double *) STARPU_MATRIX_GET_PTR(apBuffers[2]);
            mloe = (double *) STARPU_MATRIX_GET_PTR(apBuffers[3]);
            mmom = (double *) STARPU_MATRIX_GET_PTR(apBuffers[4]);
            starpu_codelet_unpack_args(apCodeletArguments, &m, &n, &m0, &n0);
            double expr2_ = 0, expr3_ = 0, expr4_ = 0;
            for (i = 0; i < m * n; i += 2) {
                expr2_ += expr2[i];
                expr3_ += expr3[i];
                expr4_ += expr4[i];
            }

            if (expr3_ == 0.0) {
                *mloe -= 1.0;
            } else {
                *mloe += (expr2_ / expr3_) - 1.0;
            }

            if (expr2_ == 0.0) {
                *mmom -= 1.0;
            } else {
                *mmom += (expr4_ / expr2_) - 1.0;
            }
        }

        static void CORE_dtrace_starpu(void *buffers[], void *cl_arg) {
            int m;
            double *A;
            double s = 0;
            double *sum = &s;
            double *trace;

            *sum = 0;
            A = (double *) STARPU_MATRIX_GET_PTR(buffers[0]);
            sum = (double *) STARPU_MATRIX_GET_PTR(buffers[1]);
            trace = (double *) STARPU_MATRIX_GET_PTR(buffers[2]);
            starpu_codelet_unpack_args(cl_arg, &m);

            double local_s = core_dtrace(A, m, trace);
            *sum+= local_s;
        }

        static double core_dtrace (double *A, int m, double* trace) {
            int i;
            double res = 0.0;
            for (i = 0; i < m; i++)
            {
                res += A[i + i * m];
                trace[i] = A[i + i * m];
            }
            return res;
        }

        static void CORE_ddotp_starpu(void *buffers[], void *cl_arg){
            int m, m0;
            double * A;
            double * dotproduct;

            dotproduct	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
            A		= (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
            starpu_codelet_unpack_args(cl_arg, &m, &m0);
            double local_dot=cblas_ddot(m, A, 1, A, 1);
            *dotproduct += local_dot;
        }

        bool recover(char *apPath, int aIterationCount, T *apTheta, T *apLogLik, int aNumParams) {

            FILE *fp;
            char *line = nullptr;
            size_t len = 0;
            int count;
            int i;
            char *pch;
            fp = fopen(apPath, "r");
            if (fp == nullptr) {
                LOGGER("cannot open observations file\n")
                exit(EXIT_FAILURE);
            }
            while (getline(&line, &len, fp) != -1) {
                pch = strtok(line, " ");
                count = (int) strtol(pch, nullptr, 10);
                if (count == aIterationCount) {
                    pch = strtok(nullptr, " ");
                    for (i = 0; i < aNumParams; i++) {
                        apTheta[i] = strtol(pch, nullptr, 10);
                        pch = strtok(nullptr, " ");
                    }
                    *apLogLik = strtol(pch, nullptr, 10);
                    fclose(fp);
                    free(line);
                    return true;
                }
            }

            fclose(fp);
            free(line);
            return false;
        }

    protected:
        /// Used Context.
        void *mpContext = nullptr;
    };

    /**
    * @brief Instantiates the Linear Algebra methods class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraMethods)

}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP