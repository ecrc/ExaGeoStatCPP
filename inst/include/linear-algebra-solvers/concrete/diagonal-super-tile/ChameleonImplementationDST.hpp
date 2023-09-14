
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDST.hpp
 * @brief This file contains the declaration of ChameleonImplementationDST class.
 * @details ChameleonImplementationDST is a concrete implementation of LinearAlgebraMethods class for diagonal super tile matrices.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace diagonalSuperTile {

            /**
             * @brief ChameleonImplementationDST is a concrete implementation of LinearAlgebraMethods class for diagonal super tile matrices.
             * @tparam T Data Type: float or double
             * 
             */
            template<typename T>
            class ChameleonImplementationDST : public LinearAlgebraMethods<T> {

            public:

                /**
                 * @brief Default constructor.
                 */
                explicit ChameleonImplementationDST() = default;

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 */
                ~ChameleonImplementationDST() override = default;

                /**
                 * @brief Initializes the descriptors necessary for the linear algebra solver.
                 * @copydoc LinearAlgebraMethods::InitiateDescriptors()
                 * 
                 */
                void InitiateDescriptors(configurations::Configurations &aConfigurations,
                                         dataunits::DescriptorData<T> &aDescriptorData,
                                         T *apMeasurementsMatrix) override;

                /**
                 * @brief Initializes the chameleon descriptors necessary for the Prediction.
                 * @copydoc LinearAlgebraMethods::InitiateDescriptors()
                 */
                void InitiatePredictionDescriptors(configurations::Configurations &aConfigurations,
                                                   dataunits::ExaGeoStatData<T> &aData) override;

                /**
                 * @brief Initializes the descriptors necessary for the Prediction Auxiliary function MLE-MLOE-MMOM.
                 * @copydoc LinearAlgebraMethods::InitiateMloeMmomDescriptors()
                 */
                void InitiateMloeMmomDescriptors(configurations::Configurations &aConfigurations,
                                                 dataunits::ExaGeoStatData<T> &aData) override;

                /**
                 * @brief Computes the covariance matrix.
                 * @copydoc LinearAlgebraMethods::CovarianceMatrixCodelet()
                 * 
                 */
                void
                CovarianceMatrixCodelet(dataunits::DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                        const int &aTriangularPart, dataunits::Locations<T> *apLocation1,
                                        dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                        T *aLocalTheta, const int &aDistanceMetric,
                                        const std::string &aKernelName) override;

                /**
                 * @brief Computes the covariance matrix.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGaussianToNonTileAsync()
                 *
                 */
                void ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData<T> &aDescriptorData, void *apDesc,
                                                      T *apTheta) override;

                /**
                 * @brief Generates the observations vector.
                 * @copydoc LinearAlgebraMethods::GenerateObservationsVector()
                 * 
                 */
                void GenerateObservationsVector(configurations::Configurations &aConfigurations,
                                                dataunits::DescriptorData<T> &aDescriptorData,
                                                const dataunits::BaseDescriptor &aDescriptor,
                                                dataunits::Locations<T> *apLocation1,
                                                dataunits::Locations<T> *apLocation2,
                                                dataunits::Locations<T> *apLocation3,
                                                const int &aDistanceMetric) override;

                /**
                 * @brief Copies the descriptor data to a double vector.
                 * @copydoc LinearAlgebraMethods::CopyDescriptorZ()
                 *
                 */
                void CopyDescriptorZ(dataunits::DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                     T *apDoubleVector) override;

                /**
                 * @brief Calculates the log likelihood value of a given value theta.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMLETile()
                */
                T
                ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &aData,
                                  configurations::Configurations &aConfigurations, const double *theta,
                                  T *apMeasurementsMatrix) override;

                /**
                 * @brief Copies a matrix in the tile layout from source to destination
                 * @copydoc LinearAlgebraMethods::ExaGeoStatLapackCopyTile()
                 */
                int ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

                /**
                * @brief Conversion from LAPACK layout to CHAM_desct_t.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatLapackToDescriptor()
                */
                int ExaGeoStatLapackToDescriptor(const common::UpperLower &aUpperLower, void *apAf77, const int &aLda,
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
                 * @copydoc LinearAlgebraMethods::ExaGeoStatPotrfTile()
                 */
                int
                ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA) override;

                /**
                 * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatTrsmTile()
                 */
                int
                ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                   const common::Trans &aTrans,
                                   const common::Diag &aDiag, const T &aAlpha, void *apA, void *apB) override;

                /**
                 * @brief Performs matrix multiplication.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGemmTile()
                 */
                int
                ExaGeoStatGemmTile(const common::Trans &aTransA, const common::Trans &aTransB, const T &aAlpha,
                                   void *apA, void *apB,
                                   const T &aBeta, void *apC) override;

                /**
                 * @brief Calculate determinant for triangular matrix.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMeasureDetTileAsync()
                 */
                int
                ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                              void *apDescDet) override;

                /**
                 * @brief opy Chameleon descriptor to vector float*.
                 * @copydoc LinearAlgebraMethods::ExaGeoStaStrideVectorTileAsync()
                 */
                int ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                   void *apSequence, void *apRequest) override;

                /**
                 * @brief Solve a positive definite linear system of equations AX = B using tiled algorithms.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatPosvTile()
                 */
                int ExaGeoStatPosvTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

                /**
                 * @brief Calculate mean square error (MSE) scalar value the prediction.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMLEMseTileAsync()
                 */
                int ExaGeoStatMLEMseTileAsync(void *apDescZPredict, void *apDescZMiss, void *apDescError,
                                              void *apSequence, void *apRequest) override;

                /**
                 * Predict missing values base on a set of given values and covariance matrix/
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMLEPredictTile()
                 */
                T *
                ExaGeoStatMLEPredictTile(exageostat::dataunits::ExaGeoStatData<T> &aData, T *apTheta,
                                         const int &aZMissNumber,
                                         const int &aZObsNumber,
                                         T *apZObs, T *apZActual, T *apZMiss,
                                         const hardware::ExaGeoStatHardware &aHardware,
                                         configurations::Configurations &aConfiguration,
                                         exageostat::dataunits::Locations<T> &aMissLocations,
                                         exageostat::dataunits::Locations<T> &aObsLocations) override;

                /**
                 * @brief Copy Lapack matrix to Descriptor Matrix
                 * @copydoc LinearAlgebraMethods::ExaGeoStatLap2Desc()
                 */
                void ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA,
                                        const common::UpperLower &aUpperLower) override;

                /**
               * @brief Copy Descriptor Matrix to Lapack matrix.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatDesc2Lap()
               */
                void ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA,
                                        const common::UpperLower &aUpperLower) override;

                /**
                 * @brief Copy the Z matrix into a pointer.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGetZObs()
                 */
                void
                ExaGeoStatGetZObs(exageostat::configurations::Configurations &aConfigurations, T *apZ, const int &aSize,
                                  exageostat::dataunits::DescriptorData<T> &aDescData,
                                  T *apMeasurementsMatrix) override;

                /**
                 * @brief Predict missing values base on a set of given values and covariance matrix.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMLEMloeMmomTile()
                 */
                void ExaGeoStatMLEMloeMmomTile(exageostat::configurations::Configurations &aConfigurations,
                                               exageostat::dataunits::ExaGeoStatData<T> &aData,
                                               const exageostat::hardware::ExaGeoStatHardware &aHardware,
                                               T *apTruthTheta,
                                               T *apEstimatedTheta, dataunits::Locations<T> &aMissLocations,
                                               dataunits::Locations<T> &aObsLocations) override;

                /**
                 * @brief Perform an asynchronous computation of MLE, MLOE, and MMOM for a tile.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMLEMloeMmomTileAsync()
                 */
                int ExaGeoStatMLEMloeMmomTileAsync(void *apDescExpr2, void *apDescExpr3, void *apDescExpr4,
                                                   void *apDescMloe, void *apDescMmom,
                                                   void *apSequence, void *apRequest) override;

                /**
                 * @brief Perform a matrix addition with scaling.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGeaddTile()
                 */
                int ExaGeoStatGeaddTile(const common::Trans &aTrans, const T &aAlpha, void *apDescA, const T &aBeta,
                                        void *apDescB) override;

                /**
                 * @brief Perform a triangular matrix multiplication.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatTrmmTile()
                */
                void ExaGeoStatTrmmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                        const common::Trans &aTrans, const common::Diag &aDiag,
                                        const T &alpha, void *apDescA, void *apDescB) override;
            };

            /**
            * @brief Instantiates the chameleon DST class for float and double types.
            * @tparam T Data Type: float or double
            *
            */
            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDST)

        }//namespace diagonalSuperTile
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP