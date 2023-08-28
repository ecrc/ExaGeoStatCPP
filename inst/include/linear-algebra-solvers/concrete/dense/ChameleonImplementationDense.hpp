
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDense.hpp
 * @brief This file contains the declaration of ChameleonImplementationDense class.
 * @details ChameleonImplementationDense is a concrete implementation of LinearAlgebraMethods class for dense matrices.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>
#include <data-generators/DataGenerator.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace dense {

            /**
             * @brief ChameleonImplementationDense is a concrete implementation of LinearAlgebraMethods class for dense matrices.
             * @tparam T Data Type: float or double
             * 
             */
            template<typename T>
            class ChameleonImplementationDense : public LinearAlgebraMethods<T> {
            public:

                /**
                 * @brief Default constructor.
                 */
                explicit ChameleonImplementationDense() = default;

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 */
                ~ChameleonImplementationDense() override = default;

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
                void GenerateObservationsVector(configurations::Configurations &aConfigurations,
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
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMleTile()
                */
                T ExaGeoStatMleTile(const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &apData,
                                    configurations::Configurations &apConfigurations, const double *theta,
                                    T *apMeasurementsMatrix) override;

                /**
                 * @brief Converts Gaussian to non-Gaussian distributed random numbers for a matrix descriptor asynchronously.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGaussianToNonTileAsync()
                */
                void ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData<T> *apDescriptorData, void *apDesc,
                                                      T *apTheta) override;

                /**
                 * @brief Copies a matrix in the tile layout from source to destination
                 * @param[in] aUpperLower Specifies the part of the matrix A to be copied to B.
                 * @param[in] apA Source matrix A.
                 * @param[in,out] apB Destination matrix B. On exit, B = A in the locations specified by UPLO.
                 * @return Successful exit
                 */
                int ExaGeoStatLapackCopyTile(common::UpperLower aUpperLower, void *apA, void *apB) override;

                /**
                * @brief Conversion from LAPACK layout to CHAM_desct_t.
                * @param[in] aUpperLower Specifies the shape of the matrix A.
                * @param[in] apAf77 LAPACK matrix.
                * @param[in] aLda The leading dimension of the matrix Af77.
                * @param[in] apA Descriptor of the CHAMELEON matrix initialized with data from Af77.
                * @return successful exit
                */
                int ExaGeoStatLapackToDescriptor(common::UpperLower aUpperLower, void *apAf77, int aLda,
                                                 void *apA) override;

                /**
                 * @brief Wait for the completion of a sequence.
                 * @param[in] apSequence apSequence A pointer to either CHAMELEON or HiCMA sequence.
                 * @return successful exit
                 */
                int
                ExaGeoStatSequenceWait(void *apSequence) override;

                /**
                 * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
                 * @param[in] aUpperLower Whether upper or lower part of the matrix A
                 * @param[in] apA Symmetric matrix A
                 * @return successful exit
                 */
                int
                ExaGeoStatPotrfTile(common::UpperLower aUpperLower, void *apA) override;

                /**
                 * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
                 * @param[in] aSide Specifies whether op(A) appears on the left or on the right of X
                 * @param[in] aUpperLower Specifies whether the matrix A is upper triangular or lower triangular.
                 * @param[in] aTrans Specifies the form of op( A ) to be used in the matrix multiplication.
                 * @param[in] aDiag Specifies whether or not A is unit triangular.
                 * @param[in] aAlpha Specifies the scalar alpha. When alpha is zero then A is not referenced and B need not be set before entry.
                 * @param[in] apA The triangular matrix A
                 * @param[in,out] apB The matrix B of dimension ,on exit is overwritten by the solution matrix X.
                 * @return successful exit
                 */
                int
                ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans,
                                   common::Diag aDiag, T aAlpha, void *apA, void *apB) override;

                /**
                 * @brief Performs matrix multiplication.
                 * @param[in] aTransA  Specifies whether the matrix A is transposed.
                 * @param[in] aTransB Specifies whether the matrix B is transposed.
                 * @param[in] aAlpha Specifies the scalar alpha.
                 * @param[in] apA Matrix A.
                 * @param[in] apB Matrix B.
                 * @param[in] aBeta Specifies the scalar beta.
                 * @param[in,out] apC On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
                 * @return successful exit.
                 */
                int
                ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha, void *apA, void *apB,
                                   T aBeta, void *apC) override;

                /**
                 * @brief Calculate determinant for triangular matrix.
                 * @param[in] apDescA Exageostat descriptor.
                 * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
                 * @param[in] apRequest Identifies this function call (for exception handling purposes).
                 * @param[in] apDescDet determinant value
                 * @return
                 */
                int
                ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                              void *apDescDet) override;

                /**
                 * @brief opy Chameleon descriptor to vector float*.
                 * @param[in] apDescA Exageostat descriptor A.
                 * @param[in] apDescB Exageostat descriptor B.
                 * @param[in] apDescC Exageostat descriptor C.
                 * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
                 * @param[in] apRequest Identifies this function call (for exception handling purposes).
                 * @return
                 */
                int ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                   void *apSequence, void *apRequest) override;

            };

            /**
            * @brief Instantiates the chameleon dense class for float and double types.
            * @tparam T Data Type: float or double
            *
            */
            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDense)

        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP