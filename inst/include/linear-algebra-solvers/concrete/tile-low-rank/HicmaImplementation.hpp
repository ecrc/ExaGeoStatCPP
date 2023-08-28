
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.hpp
 * @brief This file contains the declaration of HicmaImplementation class.
 * @details HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace tileLowRank {

            /**
             * @brief HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
             * @tparam T Data Type: float or double
             * 
             */
            template<typename T>
            class HicmaImplementation : public LinearAlgebraMethods<T> {
            public:

                /**
                 * @brief Default constructor.
                 */
                explicit HicmaImplementation() = default;

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 */
                ~HicmaImplementation() override = default;

                /**
                 * @brief Initializes the descriptors necessary for the linear algebra solver.
                 * @copydoc LinearAlgebraMethods::InitiateDescriptors()
                 * 
                 */
                void InitiateDescriptors(configurations::Configurations &aConfigurations,
                                         dataunits::DescriptorData<T> &aDescriptorData,
                                         T *apMeasurementsMatrix) override;

                /**
                 * @brief Computes the covariance matrix.
                 * @copydoc LinearAlgebraMethods::CovarianceMatrixCodelet()
                 * 
                 */
                void
                CovarianceMatrixCodelet(dataunits::DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                        int &aTriangularPart, dataunits::Locations<T> *apLocation1,
                                        dataunits::Locations<T> *apLocation2,
                                        dataunits::Locations<T> *apLocation3, T *aLocalTheta, int aDistanceMetric,
                                        exageostat::kernels::Kernel<T> *apKernel) override;

                /**
                 * @brief Generates the observations vector.
                 * @copydoc LinearAlgebraMethods::GenerateObservationsVector()
                 * 
                 */
                void GenerateObservationsVector(configurations::Configurations &apConfigurations,
                                                dataunits::DescriptorData<T> *apDescriptorData,
                                                dataunits::BaseDescriptor aDescriptor,
                                                dataunits::Locations<T> *apLocation1,
                                                dataunits::Locations<T> *apLocation2,
                                                dataunits::Locations<T> *apLocation3, int aDistanceMetric) override;

                /**
                 * @brief Copies the descriptor data to a double vector.
                 * @copydoc LinearAlgebraMethods::CopyDescriptorZ()
                 *
                 */
                void CopyDescriptorZ(dataunits::DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                     T *apDoubleVector) override;

                /**
                 * @brief Calculates the log likelihood value of a given value theta.
                 * @param aN unsigned variable used by NLOPT library.
                 * @param apTheta theta Vector with three parameter (Variance, Range, Smoothness)
                 * that is used to to generate the Covariance Matrix.
                 * @param apGrad double variable used by NLOPT library.
                 * @param apData MLE_data struct with different MLE inputs.
                 * @return log likelihood value
                */
                T
                ExaGeoStatMleTile(const hardware::ExaGeoStatHardware &apHardware, dataunits::ExaGeoStatData<T> &apData,
                                  configurations::Configurations &apConfigurations, const double *theta,
                                  T *apMeasurementsMatrix) override;

                /**
                 * @brief Calculates the log likelihood value of a given value theta.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGaussianToNonTileAsync()
                 *
                */
                void ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData<T> *apDescriptorData, void *apDesc,
                                                      T *apTheta) override;

                /**
                 * @brief Copies a matrix in the tile layout from source to destination
                 * @param aUpperLower Specifies the part of the matrix A to be copied to B.
                 * @param apA Source matrix A.
                 * @param apB Destination matrix B. On exit, B = A in the locations specified by UPLO.
                 * @return Successful exit
                 */
                int ExaGeoStatLapackCopyTile(common::UpperLower aUpperLower, void *apA, void *apB) override;

                /**
                 * @brief Conversion from LAPACK layout to HiCMA descriptor.
                 * @param aUpperLower Specifies the shape of the matrix A.
                 * @param apAf77 LAPACK matrix.
                 * @param aLda The leading dimension of the matrix Af77.
                 * @param apA Descriptor of the CHAMELEON matrix initialized with data from Af77.
                 * @return Successful exit
                 */
                int ExaGeoStatLapackToDescriptor(common::UpperLower aUpperLower, void *apAf77, int aLda,
                                                 void *apA) override;

                /**
                 * @brief Wait for the completion of a sequence.
                 * @param Identifies a set of routines sharing common exception handling.
                 * @return successful exit
                 */
                int
                ExaGeoStatSequenceWait(void *apSequence) override;

                /**
                 * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
                 * @param aUpperLower Whether upper or lower part of the matrix A
                 * @param apA Symmetric matrix A
                 * @return
                 */
                int
                ExaGeoStatPotrfTile(common::UpperLower aUpperLower, void *apA) override;

                /**
                * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
                * @param aSide Specifies whether op(A) appears on the left or on the right of X
                * @param aUpperLower Specifies whether the matrix A is upper triangular or lower triangular.
                * @param aTrans Specifies the form of op( A ) to be used in the matrix multiplication.
                * @param aDiag Specifies whether or not A is unit triangular.
                * @param aAlpha Specifies the scalar alpha. When alpha is zero then A is not referenced and B need not be set before entry.
                * @param apA The triangular matrix A
                * @param apB The matrix B of dimension ,on exit is overwritten by the solution matrix X.
                * @return successful exit
                */
                int
                ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans,
                                   common::Diag aDiag, T aAlpha, void *apA, void *apB) override;

                /**
                 * @brief Performs matrix multiplication.
                 * @param aTransA  Specifies whether the matrix A is transposed.
                 * @param aTransB Specifies whether the matrix B is transposed.
                 * @param aAlpha Specifies the scalar alpha.
                 * @param apA Matrix A.
                 * @param apB Matrix B.
                 * @param aBeta Specifies the scalar beta.
                 * @param apC On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
                 * @return successful exit.
                 */
                int
                ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha, void *apA, void *apB,
                                   T aBeta, void *apC) override;

                /**
                 * @brief Calculate determinant for triangular matrix.
                 * @param apDescA Exageostat descriptor.
                 * @param apSequence Identifies the sequence of function calls that this call belongs to.
                 * @param apRequest Identifies this function call (for exception handling purposes).
                 * @param apDescDet determinant value
                 * @return
                 */
                int
                ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                              void *apDescDet) override;

                /**
                 * @brief opy Chameleon descriptor to vector float*.
                 * @param apDescA Exageostat descriptor A.
                 * @param apDescB Exageostat descriptor B.
                 * @param apDescC Exageostat descriptor C.
                 * @param apSequence Identifies the sequence of function calls that this call belongs to.
                 * @param apRequest Identifies this function call (for exception handling purposes).
                 * @return
                 */
                int ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                   void *apSequence, void *apRequest) override;
            };

            /**
            * @brief Instantiates the Hicma TLR class for float and double types.
            * @tparam T Data Type: float or double
            *
            */
            EXAGEOSTAT_INSTANTIATE_CLASS(HicmaImplementation);

        }//namespace tileLowRank
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP