
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraMethods.hpp
 * @brief Header file for the LinearAlgebraMethods class, which defines the interface for linear algebra solvers.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
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

#include <utilities/Logger.hpp>
#include <results/Results.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <runtime/RuntimeFunctions.hpp>
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
         * @param[in] aP the P value of the kernel multiplied by time slot.
         * @param[in] apMeasurementsMatrix Pointer to the measurement matrix.
         * @return void
         *
         */
        void InitiateDescriptors(configurations::Configurations &aConfigurations,
                                 dataunits::DescriptorData<T> &aDescriptorData,
                                 const int &aP, T *apMeasurementsMatrix = nullptr);

        /**
         * @brief Initializes the descriptors necessary for the Fisher prediction function.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aDescriptorData Descriptor Data object to be populated with descriptors and data.
         * @return void
         *
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
                                           std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Initializes the descriptors necessary for the Prediction Auxiliary function MLE-MLOE-MMOM.
         * @details This method initializes the descriptors necessary for the linear algebra solver.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData DescriptorData object to be populated with descriptors and data.
         * @param[in] aP the P value of the kernel multiplied by time slot.
         * @return void
         *
         */
        void InitiateMLOEMMOMDescriptors(configurations::Configurations &aConfigurations,
                                         std::unique_ptr<ExaGeoStatData<T>> &aData, const int &aP);

        /**
         * @brief Generates synthetic data.
         * @param[in] aConfigurations The configurations object containing relevant settings.
         * @param[in,out] aData ExaGeoStatData object to be populated with synthetic data.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return void.
         *
         */
        void GenerateSyntheticData(configurations::Configurations &aConfigurations,
                                   std::unique_ptr<ExaGeoStatData<T>> &aData,
                                   const kernels::Kernel<T> &aKernel);

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
        GenerateObservationsVector(configurations::Configurations &aConfigurations,
                                   std::unique_ptr<ExaGeoStatData<T>> &aData,
                                   dataunits::Locations<T> *apLocation1, dataunits::Locations<T> *apLocation2,
                                   dataunits::Locations<T> *apLocation3, const int &aDistanceMetric,
                                   const kernels::Kernel<T> &aKernel);

        /**
         * @brief Calculates the log likelihood value of a given value theta.
         * @param[in,out] aData DescriptorData object to be populated with descriptors and data.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] apTheta Optimization parameter used by NLOPT.
         * @param[in] apMeasurementsMatrix measurements matrix to be stored in DescZ.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return log likelihood value
         *
         */
        virtual T ExaGeoStatMLETile(std::unique_ptr<ExaGeoStatData<T>> &aData,
                                    configurations::Configurations &aConfigurations, const double *apTheta,
                                    T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) = 0;


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
        *
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
         *
         */
        virtual void ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                        const common::Trans &aTrans, const common::Diag &aDiag, const T &aAlpha,
                                        void *apA, void *apCD, void *apCrk, void *apZ, const int &aMaxRank) = 0;


        /**
         * @brief Solve a positive definite linear system of equations AX = B using tiled algorithms.
         * @param[in] aUpperLower Specifies whether the matrix A is upper triangular or lower triangular.
         * @param [in] apA coefficient matrix of the system of linear equations. This matrix is expected to be positive definite.
         * @param [in] apB Pointer to coefficient matrix of the system of linear equations. This matrix is expected to be positive definite.
         * @return void
         *
         */
        void ExaGeoStatPosvTile(const common::UpperLower &aUpperLower, void *apA, void *apB);

        /**
         * @brief Predict missing values base on a set of given values and covariance matrix/
         * @param[in] aData Reference to Data containing different MLE inputs.
         * @param[in] apTheta theta Vector with three parameter (Variance, Range, Smoothness) that is used to to generate the Covariance Matrix.
         * @param[in] aZMissNumber number of missing values (unknown observations).
         * @param[in] aZObsNumber number of observed values (known observations).
         * @param[in] apZObs observed values vector (known observations).
         * @param[in] apZActual actual missing values vector (in the case of testing MSPE).
         * @param[in] apZMiss missing values vector (unknown observations).
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] aMissLocations Reference to Locations object containing missed locations.
         * @param[in] aObsLocations Reference to Locations object containing observed locations.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return the prediction Mean Square Error (MSPE).
         *
         */
        T *ExaGeoStatMLEPredictTile(std::unique_ptr<ExaGeoStatData<T>> &aData, T *apTheta,
                                    const int &aZMissNumber,
                                    const int &aZObsNumber, T *apZObs, T *apZActual, T *apZMiss,
                                    configurations::Configurations &aConfiguration,
                                    exageostat::dataunits::Locations<T> &aMissLocations,
                                    exageostat::dataunits::Locations<T> &aObsLocations,
                                    const kernels::Kernel<T> &aKernel);

        /**
         * @brief Predict missing values base on a set of given values and Non-Gaussian covariance matrix/
         * @param[in] aData Reference to Data containing different MLE inputs.
         * @param[in] apTheta theta Vector with three parameter (Variance, Range, Smoothness) that is used to to generate the Covariance Matrix.
         * @param[in] aZMissNumber number of missing values (unknown observations).
         * @param[in] aZObsNumber number of observed values (known observations).
         * @param[in] apZObs observed values vector (known observations).
         * @param[in] apZActual actual missing values vector (in the case of testing MSPE).
         * @param[in] apZMiss missing values vector (unknown observations).
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] aMissLocations Reference to Locations object containing missed locations.
         * @param[in] aObsLocations Reference to Locations object containing observed locations.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return the prediction Mean Square Error (MSPE).
         *
         */
        T *ExaGeoStatMLENonGaussianPredictTile(std::unique_ptr<ExaGeoStatData<T>> &aData, T *apTheta,
                                               const int &aZMissNumber,
                                               const int &aZObsNumber, T *apZObs, T *apZActual, T *apZMiss,
                                               configurations::Configurations &aConfiguration,
                                               exageostat::dataunits::Locations<T> &aMissLocations,
                                               exageostat::dataunits::Locations<T> &aObsLocations,
                                               const kernels::Kernel<T> &aKernel);

        /**
         * @brief Copy Lapack matrix to Descriptor Matrix.
         * @param[in] apA Lapack Matrix.
         * @param[in] aLDA Size.
         * @param[out] apDescA Matrix Descriptor.
         * @param[in] aUpperLower Specifies Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @return void
         *
         */
        void ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const common::UpperLower &aUpperLower);

        /**
         * @brief Copy Descriptor Matrix to Lapack matrix.
         * @param[out] apA Lapack Matrix.
         * @param[in] aLDA Size.
         * @param[in] apDescA Matrix Descriptor
         * @param[in] aUpperLower Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @return void
         *
         */
        void ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA, const common::UpperLower &aUpperLower);

#ifdef USE_HICMA

        /**
         * @brief Copy Descriptor Matrix to another Descriptor matrix.
         * @param[out] apSourceDesc Descriptor matrix to be copied
         * @param[out] apDestinationDesc Descriptor matrix to be copied to.
         * @param[in] aSize Size of matrix to be copied.
         * @param[in] aDirection Specifies the type of Descriptors to be copied.
         * @return void
         *
         */
        void CopyDescriptors(void *apSourceDesc, void *apDestinationDesc, const int &aSize,
                             const common::CopyDirection &aDirection);

#endif

        /**
         * @brief Sets the values of all or part of a two-dimensional Tile.
         * @param[in] aUpperLower Specifies Specifies whether the upper or lower triangular part of the covariance matrix is stored.
         * @param[in] aAlpha All the off diagonal array elements are set to aAlpha.
         * @param[in] aBeta All the diagonal array elements are set to aBeta.
         * @param[out] apDescriptor Pointer to matrix descriptor to be set with aAlpha and aBeta.
         * @return void
         *
         */
        void ExaGeoStatLaSetTile(const common::UpperLower &aUpperLower, T aAlpha, T aBeta, void *apDescriptor);

        /**
         * @brief Copy the Z matrix into a pointer.
         * @param[out] apZ Pointer to an array to copy Z matrix into.
         * @param[in] aSize Size of the matrix.
         * @param[in] aDescData Descriptor data containing required Z matrix Descriptor.
         * @param[in] aP the P value of the kernel multiplied by time slot.
         * @return void
         *
         */
        void ExaGeoStatGetZObs(configurations::Configurations &aConfigurations, T *apZ, const int &aSize,
                               exageostat::dataunits::DescriptorData<T> &aDescData, T *apMeasurementsMatrix,
                               const int &aP);

        /**
         * @brief Predict missing values based on a set of given values and covariance matrix.
         * @details This function predicts missing values using the maximum likelihood estimation (MLE),
         * maximum likelihood on the empirical orthogonal functions (MLOE), and method of moments (MMOM).
         * @param[in] aConfigurations Configurations for the prediction.
         * @param[in, out] aData Data for prediction (input and output).
         * @param[in] apTruthTheta Pointer to the true theta values.
         * @param[in] apEstimatedTheta Pointer to the estimated theta values.
         * @param[in] aMissLocations Locations of missing values.
         * @param[in] aObsLocations Locations of observed values.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return void
         *
         */
        void ExaGeoStatMLETileMLOEMMOM(configurations::Configurations &aConfigurations,
                                       std::unique_ptr<ExaGeoStatData<T>> &aData,
                                       T *apTruthTheta, T *apEstimatedTheta, dataunits::Locations<T> &aMissLocations,
                                       dataunits::Locations<T> &aObsLocations, const kernels::Kernel<T> &aKernel);

        /**
         * @brief Maximum Likelihood Evaluation (MLE) Fisher method.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         * @param[in] apTheta Pointer containing three parameter (Variance, Range, Smoothness) that is used to to generate the Covariance Matrix.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return Fisher Matrix
         *
         */
        T *
        ExaGeoStatFisherTile(configurations::Configurations &aConfigurations,
                             std::unique_ptr<ExaGeoStatData<T>> &aData,
                             T *apTheta, const kernels::Kernel<T> &aKernel);

        /**
         * @brief Perform a matrix addition with scaling.
         * @details This function performs a matrix addition with scaling, given the matrices A and B.
         * @param[in] aTrans Specifies whether to transpose matrix A.
         * @param[in] aAlpha Scaling factor for matrix A.
         * @param[in] apDescA Descriptor for matrix A.
         * @param[in] aBeta Scaling factor for matrix B.
         * @param[in] apDescB Descriptor for matrix B.
         * @return void
         *
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
         *
         */
        void ExaGeoStatTrmmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                const common::Trans &aTrans, const common::Diag &aDiag, const T &alpha, void *apDescA,
                                void *apDescB);

        /**
         * @brief Recovers theta and log-likelihood from a file.
         * @param[in] apPath A pointer to the path of the file from which to recover the data.
         * @param[in] aIterationCount The iteration count to look for in the file.
         * @param[in,out] apTheta A pointer to the array where the theta values will be stored.
         * @param[in,out] apLogLik A pointer to the variable where the log-likelihood value will be stored.
         * @param[in] aNumParams The number of parameters (elements) in the theta array.
         * @return bool `true` if the specified iteration count is found and successfully parsed, `false` otherwise.
         *
         */
        bool Recover(char *apPath, const int &aIterationCount, T *apTheta, T *apLogLik, const int &aNumParams);
    };

    /**
    * @brief Instantiates the Linear Algebra methods class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraMethods)

}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP
