
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.hpp
 * @brief Header file for the LinearAlgebraMethods class, which defines the interface for linear algebra solvers.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
 *
 *  This header file defines the abstract class LinearAlgebraMethods, which provides an interface for linear algebra solvers.
 *  The purpose of this interface is to allow different concrete linear algebra solvers to be interchangeable,
 *  so that they can be used interchangeably by other parts of the software system that rely on linear algebra.
 *
**/

#ifndef EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP
#define EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP

#include <vector>

#include <common/Definitions.hpp>
#include <kernels/Kernel.hpp>
#include <configurations/Configurations.hpp>
#include <common/Utils.hpp>
#include <helpers/DiskWriter.hpp>

extern "C" {
#include <gsl/gsl_errno.h>
}

#define starpu_mpi_codelet(_codelet_) _codelet_


namespace exageostat {
    namespace linearAlgebra {

        /**
         * @class LinearAlgebraMethods
         * @brief A class that defines the interface for linear algebra solvers.
         * @tparam T The data type of the linear algebra solver.
         */
        template<typename T>
        class LinearAlgebraMethods {
        public:

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~LinearAlgebraMethods() = default;

            /**
             * @brief Initializes the descriptors necessary for the linear algebra solver.
             *
             * This method initializes the descriptors necessary for the linear algebra solver.
             */
            virtual void InitiateDescriptors() = 0;

            /**
             * @brief Destroys the descriptors used by the linear algebra solver.
             *
             * This method destroys the descriptors used by the linear algebra solver.
             */
            virtual void DestoryDescriptors() = 0;

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
            virtual void
            CovarianceMatrixCodelet(void *descA, int &uplo, dataunits::Locations *apLocation1,
                                    dataunits::Locations *apLocation2,
                                    dataunits::Locations *apLocation3, double *aLocalTheta, int aDistanceMetric,
                                    exageostat::kernels::Kernel *apKernel) = 0;

            virtual void CopyDescriptorZ(void *apDescA, double *apDoubleVector) = 0;

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
            virtual void GenerateObservationsVector(void *descA, dataunits::Locations *apLocation1,
                                                    dataunits::Locations *apLocation2,
                                                    dataunits::Locations *apLocation3, std::vector<double> aLocalTheta,
                                                    int aDistanceMetric, exageostat::kernels::Kernel *apKernel) = 0;

            /**
             * @brief Initializes the context for the linear algebra solver with the specified number of cores and GPUs.
             *
             * @param[in] apCoresNumber The number of cores to use for the solver.
             * @param[in] apGPUs The number of GPUs to use for the solver.
             */
            virtual void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) = 0;

            /**
             * @brief Finalizes the context for the linear algebra solver.
             */
            virtual void ExaGeoStatFinalizeContext() = 0;

            /**
             * @brief Sets the configurations for the linear algebra solver.
             *
             * @param[in] apConfigurations A pointer to the configurations for the solver.
             */
            void SetConfigurations(configurations::Configurations *apConfigurations) {
                this->mpConfigurations = apConfigurations;
            }

            /**
             * @brief Gets the matrix.
             *
             * @return Pointer to the matrix.
             */
            double *GetMatrix() {
                if (this->apMatrix == nullptr) {
                    throw std::runtime_error("Matrix is null");
                }
                return this->apMatrix;
            }

            /**
             * @brief allocates matrix tile.
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
            virtual void
            ExageostatAllocateMatrixTile(void **apDescriptor, bool aIsOOC, T *apMemSpace, int aType2, int aMB,
                                         int aNB, int aMBxNB, int aLda, int aN, int aSMB, int aSNB, int aM, int aN2,
                                         int aP, int aQ) = 0;

            static void cl_dcmg_cpu_func(void *buffers[], void *cl_arg) {

                int m, n, m0, n0;
                exageostat::dataunits::Locations *apLocation1;
                exageostat::dataunits::Locations *apLocation2;
                exageostat::dataunits::Locations *apLocation3;
                double *theta;
                double *A;
                int distance_metric;
                exageostat::kernels::Kernel *kernel;

                A = (double *) STARPU_MATRIX_GET_PTR(buffers[0]);

                starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &apLocation1, &apLocation2, &apLocation3, &theta,
                                           &distance_metric, &kernel);
                kernel->GenerateCovarianceMatrix(A, m, n, m0, n0, apLocation1,
                                                 apLocation2, apLocation3, theta, distance_metric);
            }

            //// These codlets and structs will be added to another level of abstraction and interface of runtime system. This is a quick fix for now.
            //// TODO: Create a Factory for Runtime system.
            struct starpu_codelet cl_dcmg =
                    {
                            .where        = STARPU_CPU,
                            .cpu_func     = cl_dcmg_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)
                            //    .cuda_func      = {cl_dcmg_cuda_func},
#endif
                            .nbuffers     = 1,
                            .modes        = {STARPU_W},
                            .name         = "dcmg"
                    };

            static void CORE_dzcpy_starpu(void *buffers[], void *cl_arg) {
                int m;
                double *A;
                int m0;
                double *r;

                A = (double *) STARPU_MATRIX_GET_PTR(buffers[0]);
                starpu_codelet_unpack_args(cl_arg, &m, &m0, &r);
                memcpy(A, &r[m0], m * sizeof(double));
            }

            struct starpu_codelet cl_dzcpy =
                    {
                            .where        = STARPU_CPU,
                            .cpu_funcs    = {CORE_dzcpy_starpu},
                            .nbuffers    = 1,
                            .modes        = {STARPU_W},
                            .name        = "dzcpy"
                    };

        protected:
            //// Used configurations map.
            configurations::Configurations *mpConfigurations = nullptr;
            //// Used Matrix
            double *apMatrix = nullptr;
        };

        /**
         * @brief Instantiates the LinearAlgebraMethods class for float and double types.
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraMethods)

    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP