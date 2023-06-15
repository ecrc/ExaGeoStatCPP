
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.hpp
 * @brief This file contains the declaration of HicmaImplementation class.
 * HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
 * @version 1.0.0
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace tileLowRank {

            /**
             * @brief
             * HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
             *
             * @tparam T Type of matrix elements.
             */
            template<typename T>
            class HicmaImplementation : public LinearAlgebraMethods<T>{
            public:

                /**
                 * @brief
                 * Default constructor.
                 */
                explicit HicmaImplementation() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 */
                ~HicmaImplementation() override = default;

                /**
                 * @brief Initializes the descriptors necessary for the linear algebra solver.
                 *
                 * This method initializes the descriptors necessary for the linear algebra solver.
                 */
                void InitiateDescriptors() override;

                /**
                 * @brief Destroys the descriptors used by the linear algebra solver.
                 *
                 * This method destroys the descriptors used by the linear algebra solver.
                 */
                void DestoryDescriptors() override;

                /**
                 * @brief Computes the covariance matrix.
                 *
                 * @param[in] descA Pointer to the descriptor for the covariance matrix.
                 * @param[in] uplo Specifies whether the upper or lower triangular part of the covariance matrix is stored.
                 * @param[in] apLocation1 Pointer to the first set of locations.
                 * @param[in] apLocation2 Pointer to the second set of locations.
                 * @param[in] apLocation3 Pointer to the third set of locations.
                 * @param[in] aLocalTheta Pointer to the local theta values.
                 * @param[in] aDistanceMetric Specifies the distance metric to use.
                 * @param[in] apKernel Pointer to the kernel function to use.
                 */
                void
                CovarianceMatrixCodelet(void *descA, int &uplo, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                        dataunits::Locations *apLocation3, double* apLocalTheta, int aDistanceMetric,
                                        exageostat::kernels::Kernel *apKernel) override;

                void CopyDescriptorZ(void *apDescA, double *apDoubleVector) override;
                /**
                 * @brief Generates the observations vector.
                 *
                 * @param[in] descA Pointer to the descriptor for the observations vector.
                 * @param[in] apLocation1 Pointer to the first set of locations.
                 * @param[in] apLocation2 Pointer to the second set of locations.
                 * @param[in] apLocation3 Pointer to the third set of locations.
                 * @param[in] aLocalTheta Pointer to the local theta values.
                 * @param[in] aDistanceMetric Specifies the distance metric to use.
                 * @param[in] apKernel Pointer to the kernel function to use.
                 */
                void GenerateObservationsVector(void *descA, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                                dataunits::Locations *apLocation3, std::vector<double> aLocalTheta, int aDistanceMetric, exageostat::kernels::Kernel * apKernel) override;

                /**
                 * @brief
                 * Initializes the context needed for the HICMA solver.
                 *
                 * @param apCoresNumber Number of cores to allocate.
                 * @param apGPUs Number of GPUs to allocate.
                 */
                void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) override;

                /**
                 * @brief
                 * Finalizes the context needed for the HICMA solver.
                 */
                void ExaGeoStatFinalizeContext() override;

                /**
                * @brief allocates approx matrix tile.
                *
                * @param[in,out] apDescriptor The descriptor for the tile.
                * @param[in] aIsOOC Whether the matrix is out-of-core.
                * @param[in] apMemSpace The memory space to use for the tile.
                * @param[in] aType2 The data type of the tile.
                * @param[in] aMB The row block size of the tile.
                * @param[in] aNB The column block size of the tile.
                * @param[in] aMBxNB The product of row and column block sizes.
                * @param[in] aLda The leading dimension of the tile.
                * @param[in] aN The total number of columns of the matrix.
                * @param[in] aSMB The row block size for the matrix distribution.
                * @param[in] aSNB The column block size for the matrix distribution.
                * @param[in] aM The total number of rows of the matrix.
                * @param[in] aN2 The total number of columns of the matrix after padding.
                * @param[in] aP The row coordinate of the tile.
                * @param[in] aQ The column coordinate of the tile.
                */
                void ExageostatAllocateMatrixTile(void ** apDescriptor, bool aIsOOC, T* apMemSpace, int aType2, int aMB,
                                                  int aNB, int aMBxNB, int aLda, int aN, int aSMB, int aSNB, int aM, int aN2, int aP, int aQ) override;


            private:
                //// Used context
                static void *apContext;

            };
            /**
             * @brief Instantiates the LinearAlgebraMethods class for float and double types.
             */
            EXAGEOSTAT_INSTANTIATE_CLASS(HicmaImplementation);

        }//namespace tileLowRank
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP