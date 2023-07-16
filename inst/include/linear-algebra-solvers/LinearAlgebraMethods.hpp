
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.hpp
 * @brief Header file for the LinearAlgebraMethods class, which defines the interface for linear algebra solvers.
 * @version 1.0.0
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

#include <common/Definitions.hpp>
#include <kernels/Kernel.hpp>
#include <configurations/Configurations.hpp>
#include <configurations/data-modeling/DataModelingConfigurations.hpp>
#include <common/Utils.hpp>
#include <helpers/DiskWriter.hpp>

namespace exageostat {
    namespace linearAlgebra {

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
             * @return void
             *
             */
            virtual void InitiateDescriptors() = 0;

            /**
             * @brief Destroys the descriptors used by the linear algebra solver.
             * @details This method destroys the descriptors used by the linear algebra solver.
             * @return void
             *
             */
            virtual void DestoryDescriptors() = 0;

            /**
             * @brief Computes the covariance matrix.
             * @param[in] apDescriptor Pointer to the descriptor for the covariance matrix.
             * @param[in] aTriangularPart Specifies whether the upper or lower triangular part of the covariance matrix is stored.
             * @param[in] apLocation1 Pointer to the first set of locations.
             * @param[in] apLocation2 Pointer to the second set of locations.
             * @param[in] apLocation3 Pointer to the third set of locations.
             * @param[in] aLocalTheta Pointer to the local theta values.
             * @param[in] aDistanceMetric Specifies the distance metric to use.
             * @param[in] apKernel Pointer to the kernel function to use.
             * @return void
             *
             */
            virtual void
            CovarianceMatrixCodelet(void *apDescriptor, int &aTriangularPart, dataunits::Locations *apLocation1,
                                    dataunits::Locations *apLocation2,
                                    dataunits::Locations *apLocation3, double *aLocalTheta, int aDistanceMetric,
                                    exageostat::kernels::Kernel *apKernel) = 0;

            /**
             * @brief Copies the descriptor data to a double vector.
             * @param[in] apDescriptor Pointer to the descriptor data.
             * @param[in,out] apDoubleVector Pointer to the double vector to copy the descriptor data to.
             * @return void
             *
             */
            virtual void CopyDescriptorZ(void *apDescriptor, double *apDoubleVector) = 0;

            /**
             * @brief Generates the observations vector.
             * @param[in] apDescriptor Pointer to the descriptor for the observations vector.
             * @param[in] apLocation1 Pointer to the first set of locations.
             * @param[in] apLocation2 Pointer to the second set of locations.
             * @param[in] apLocation3 Pointer to the third set of locations.
             * @param[in] aLocalTheta Pointer to the local theta values.
             * @param[in] aDistanceMetric Specifies the distance metric to use.
             * @param[in] apKernel Pointer to the kernel function to use.
             * @return void
             *
             */
            virtual void GenerateObservationsVector(void *apDescriptor, dataunits::Locations *apLocation1,
                                                    dataunits::Locations *apLocation2,
                                                    dataunits::Locations *apLocation3, std::vector<double> aLocalTheta,
                                                    int aDistanceMetric, exageostat::kernels::Kernel *apKernel) = 0;

            /**
             * @brief Initializes the context for the linear algebra solver with the specified number of cores and GPUs.
             * @param[in] apCoresNumber The number of cores to use for the solver.
             * @param[in] apGPUs The number of GPUs to use for the solver.
             * @return void
             *
             */
            virtual void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) = 0;

            /**
             * @brief Finalizes the context for the linear algebra solver.
             *
             */
            virtual void ExaGeoStatFinalizeContext() = 0;

            /**
             * @brief Sets the configurations for the linear algebra solver.
             * @param[in] apConfigurations A pointer to the configurations for the solver.
             * @return void
             *
             */
            void SetConfigurations(configurations::Configurations *apConfigurations) {
                this->mpConfigurations = apConfigurations;
            }

            /**
             * @brief Gets the matrix.
             * @return Pointer to the matrix.
             *
             */
            double *GetMatrix() {
                if (this->apMatrix == nullptr) {
                    throw std::runtime_error("Matrix is null");
                }
                return this->apMatrix;
            }

            /**
             * @brief allocates matrix tile.
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
             * @return void
             *
             */
            virtual void
            ExageostatAllocateMatrixTile(void **apDescriptor, bool aIsOOC, T *apMemSpace, int aType2, int aMB,
                                         int aNB, int aMBxNB, int aLda, int aN, int aSMB, int aSNB, int aM, int aN2,
                                         int aP, int aQ) = 0;

            /**
             * @brief Calculates the log likelihood value of a given value theta.
             * @param aN unsigned variable used by NLOPT library.
             * @param apTheta theta Vector with three parameter (Variance, Range, Smoothness)
             * that is used to to generate the Covariance Matrix.
             * @param apGrad double variable used by NLOPT library.
             * @param apData MLE_data struct with different MLE inputs.
             * @return log likelihood value
             */
            virtual T
            ExageostatMleTile(std::vector<double> &apTheta, dataunits::Locations *apDataLocations, configurations::data_modeling::DataModelingConfigurations *apDataModelingConfiguration) = 0;

            /**
             * @brief Copies a matrix in the tile layout from source to destination
             * @param aUpperLower Specifies the part of the matrix A to be copied to B.
             * @param apA Source matrix A.
             * @param apB Destination matrix B. On exit, B = A in the locations specified by UPLO.
             * @return Successful exit
             */
            virtual int
            ExageostatLacpyTile(common::UpperLower aUpperLower, void *apA, void *apB) = 0;

            /**
             * @brief Conversion from LAPACK layout to Exageostat descriptor.
             * @param aUpperLower Specifies the shape of the matrix A.
             * @param apAf77 LAPACK matrix.
             * @param aLda The leading dimension of the matrix Af77.
             * @param apA Descriptor of the CHAMELEON matrix initialized with data from Af77.
             * @return
             */
            virtual int
            ExageostatLap2Desc(common::UpperLower aUpperLower, void *apAf77, int aLda, void * apA) = 0;

            /**
             * @brief Wait for the completion of a sequence.
             * @param Identifies a set of routines sharing common exception handling.
             * @return successful exit
             */
            virtual int
            ExageostatSequenceWait(void * apSequence) = 0;

            /**
             * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
             * @param aUpperLower Whether upper or lower part of the matrix A
             * @param apA Symmetric matrix A
             * @return
             */
            virtual int
            ExageostatDpotrfTile(common::UpperLower aUpperLower, void *apA) = 0;

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
            virtual int
            ExageostatDtrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA, void *apB) = 0;

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
            virtual int
            ExageostatDgemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha, void *apA, void *apB, T aBeta, void * apC) = 0;

            /**
             * @brief Calculate determinant for triangular matrix.
             * @param apDescA Exageostat descriptor.
             * @param apSequence Identifies the sequence of function calls that this call belongs to.
             * @param apRequest Identifies this function call (for exception handling purposes).
             * @param apDescDet determinant value
             * @return
             */
            virtual int
            ExageostatMleMdetTileAsync(void *apDescA, void * apSequence, void *apRequest, void *apDescDet) = 0;

            /**
             * @brief opy Chameleon descriptor to vector float*.
             * @param apDescA Exageostat descriptor A.
             * @param apDescB Exageostat descriptor B.
             * @param apDescC Exageostat descriptor C.
             * @param apSequence Identifies the sequence of function calls that this call belongs to.
             * @param apRequest Identifies this function call (for exception handling purposes).
             * @return
             */
            virtual int ExageostaStrideVecTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                 void * apSequence, void *apRequest) = 0;

            /**
             * @briefCodelet to generate covariance matrix in buffersiptor
             * descA in  dense format between two sets of locations
             * @param aUpperLower Upper or lower fill of the matrix.
             * @param apDescA Chameleon buffersiptor that handles the generated covariance matrix.
             * @param apL1 Location struct of the first input.
             * @param apL2 Location struct of the second input.
             * @param apLm
             * @param apTheta
             * @param apDm Distance metric "euclidean Distance ("ED" -->0) or "Great Circle Distance (GCD) -->1".
             * @param apKernelFun
             * @param sequence  Identifies the sequence of function calls that this call belongs to
                                (for completion checks and exception handling purposes).
             * @param request Identifies this function call (for exception handling purposes).
             * @return
             */
            virtual int ExageostatMleDcmgTileAsync(common::UpperLower aUpperLower, void *apDescA, dataunits::Locations *apL1,
                                                   dataunits::Locations *apL2, dataunits::Locations *apLm, T* apTheta,
                                                   std::string &aDm, std::string &aKernelFun,
                                                   void *apSequence, void *apRequest) = 0;

            //// These codlets and structs will be added to another level of abstraction and interface of runtime system. This is a quick fix for now.
            //// TODO: Create a Factory for Runtime system.

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

            static void CORE_dmdet_starpu(void *buffers[], void *cl_arg) {
                int m;
                int n;
                double* A;
                int m0;
                int n0;
                double det = 0;
                double* determinant = &det;

                *determinant = 0;
                A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
                determinant = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
                starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);
                double local_det = core_dmdet(A, m, n, m0, n0);
                *determinant += local_det;
            }

            static double core_dmdet(double* A, int m, int n, int m0, int n0) {

                int i;
                double res = 0.0;
                for (i = 0; i < m; i++) {
                    if (A[i + i * m] > 0)
                        res += log(A[i + i * m]);
                }
                return res;
            }

            static void CORE_stride_vecstarpu(void *buffers[], void *cl_arg) {
                int m;
                int tempmm;
                double* A;
                double* B;
                double* C;
                int m0;
                int i = 0;
                int j = 0;

                A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
                B = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
                C = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
                starpu_codelet_unpack_args(cl_arg, &tempmm, &m0, &m);
                j = 0;
                for (i = 0; i < tempmm - 1; i += 2) {
                    B[j] = A[i];
                    C[j] = A[i + 1];
                    j++;
                }
            }

            struct starpu_codelet cl_dzcpy =
                    {
                            .where        = STARPU_CPU,
                            .cpu_funcs    = {CORE_dzcpy_starpu},
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
                            .cpu_funcs      = {CORE_stride_vecstarpu},
                            .nbuffers       = 3,
                            .modes          = {STARPU_R, STARPU_W, STARPU_W},
                            .name           = "stride_vec"
                    };

            bool recover(char *path, int iter_count, T* theta, T* loglik, int num_params) {

                FILE *fp;
                char *line = NULL;
                size_t len = 0;
                ssize_t read;
                int count = 0;
                int i = 0;
                char *pch;
                fp = fopen(path, "r");
                if (fp == NULL) {
                    printf("cannot open observations file\n");
                    exit(EXIT_FAILURE);
                }

                while ((read = getline(&line, &len, fp)) != -1) {
                    pch = strtok(line, " ");
                    count = atoi(pch);
                    if (count == iter_count) {
                        pch = strtok(NULL, " ");
                        for (i = 0; i < num_params; i++) {
                            theta[i] = atof(pch);
                            pch = strtok(NULL, " ");
                        }
                        *loglik = atof(pch);
                        fclose(fp);
                        free(line);

                        return true;
                    }
                    count++;
                }

                fclose(fp);
                free(line);

                return false;
            }

        protected:
            //// Used configurations map.
            configurations::Configurations *mpConfigurations = nullptr;
            //// Used Matrix
            double *apMatrix = nullptr;
        };

        /**
        * @brief Instantiates the Linear Algebra methods class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraMethods)

    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP