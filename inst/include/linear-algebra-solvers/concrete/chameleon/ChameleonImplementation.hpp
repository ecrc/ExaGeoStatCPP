
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.hpp
 * @brief This file contains the declaration of ChameleonImplementation class.
 * @details ChameleonImplementation is a concrete implementation of the LinearAlgebraMethods class for the common functionality implementation shared between dense and diagonal-super tile matrices.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {

            /**
             * @brief ChameleonImplementation is a concrete implementation of LinearAlgebraMethods class for dense or diagonal-super tile matrices.
             * @tparam T Data Type: float or double
             *
             */
            template<typename T>
            class ChameleonImplementation : public LinearAlgebraMethods<T> {
            public:

                /**
                 * @brief Calculates the log likelihood value of a given value theta.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatMLETile()
                 */
                T
                ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &aData,
                                  configurations::Configurations &aConfigurations, const double *theta,
                                  T *apMeasurementsMatrix) override;

                /**
                 * @brief Copy Lapack matrix to Descriptor Matrix
                 * @copydoc LinearAlgebraMethods::ExaGeoStatLap2Desc()
                 */
                void
                ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA,
                                        const common::UpperLower &aUpperLower) override;

                /**
                 * @brief Copy Descriptor Matrix to Lapack matrix.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatDesc2Lap()
                 */
                void
                ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA,
                                        const common::UpperLower &aUpperLower) override;

                /**
                 * @brief Copies a matrix in the tile layout from source to destination
                 * @copydoc LinearAlgebraMethods::ExaGeoStatLapackCopyTile()
                 */
                void
                ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

                /**
                 * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatDataGetAddr()
                 */
                void * ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) override;

                /**
                  * @brief Sets the values of all or part of a two-dimensional Tile.
                  * @copydoc LinearAlgebraMethods::ExaGeoStatLaSetTile()
                  */
                void ExaGeoStatLaSetTile(const common::UpperLower &aUpperLower, T alpha, T beta, void *apDescriptor) override;

                /**
                 * @brief Solve a positive definite linear system of equations AX = B using tiled algorithms.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatPosvTile()
                 */
                void
                ExaGeoStatPosvTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

                /**
                 * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatTrsmTile()
                 */
                void
                ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                   const common::Trans &aTrans,
                                   const common::Diag &aDiag, const T &aAlpha, void *apA, void *apB) override;


                /**
                 * @brief Perform a triangular matrix multiplication.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatTrmmTile()
                */
                void
                ExaGeoStatTrmmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                        const common::Trans &aTrans, const common::Diag &aDiag,
                                        const T &alpha, void *apDescA, void *apDescB) override;


                /**
                 * @brief Perform a matrix addition with scaling.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatGeaddTile()
                 */
                void
                ExaGeoStatGeaddTile(const common::Trans &aTrans, const T &aAlpha, void *apDescA, const T &aBeta,
                                    void *apDescB) override;

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
                int
                ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apSequence,
                                                   void *apRequest) override;

                /**
                 * @brief Copy Chameleon descriptor to vector float*.
                 * @copydoc LinearAlgebraMethods::ExaGeoStaStrideVectorTileAsync()
                 */
                int
                ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apDescD,
                                                   void *apSequence, void *apRequest) override;

                /**
                 * @brief Wait for the completion of a sequence.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatSequenceWait()
                 */
                void
                ExaGeoStatSequenceWait(void *apSequence) override;

                /**
                 * @brief Create CHAMELEON Sequence.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatCreateSequence()
                 */
                void
                ExaGeoStatCreateSequence(void * apSequence) override;

                /**
                 * @brief Initialize the runtime option structure for CHAMELEON
                 * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsInit()
                 */
                void
                ExaGeoStatOptionsInit(void *apOptoins, void * apContext, void * apSequence, void * apRequest) override;

                /**
                 * @brief Submit the release of the workspaces associated to the options structure.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsFree()
                 */
                void
                ExaGeoStatOptionsFree(void *apOptions) override;

                /**
                 * @brief Finalize the runtime option structure for CHAMELEON.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsFinalize()
                 */
                void
                ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) override;

            };
            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementation)
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP