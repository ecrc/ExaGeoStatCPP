
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

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>

#include <common/Definitions.hpp>
#include <common/Utils.hpp>
#include <configurations/Configurations.hpp>
#include <helpers/DiskWriter.hpp>
#include <kernels/Kernel.hpp>
#include <data-units/DescriptorData.hpp>
#include <data-units/ExaGeoStatData.hpp>

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
             * @param[in] aConfigurations Configurations object containing relevant settings.
             * @param[in,out] aDescriptorData DescriptorData object to be populated with descriptors and data.
             * @param[in] apMeasurementsMatrix Pointer to the measurement matrix.
             * @return void
             *
             */
            virtual void InitiateDescriptors(configurations::Configurations &aConfigurations,
                                             dataunits::DescriptorData<T> &aDescriptorData,
                                             T *apMeasurementsMatrix = nullptr) = 0;

            /**
             * @brief Generates synthetic data.
             * @param[in] apConfigurations The configurations object containing relevant settings.
             * @param[in] apHardware  ExaGeoStatHardware object representing the hardwar.
             * @param[in,out] apData ExaGeoStatData object to be populated with synthetic data.
             * @param[in] aType The descriptor type.
             * @return None.
             *
             */
            void GenerateSyntheticData(configurations::Configurations &aConfigurations,
                                       const hardware::ExaGeoStatHardware &aHardware,
                                       dataunits::ExaGeoStatData<T> &apData, const common::DescriptorType aType) {
                this->mpContext = aHardware.GetContext();
                this->InitiateDescriptors(aConfigurations, *apData.GetDescriptorData());
                auto descC = apData.GetDescriptorData()->GetDescriptor(aType, common::DESCRIPTOR_C);
                this->GenerateObservationsVector(aConfigurations, apData.GetDescriptorData(), descC,
                                                 apData.GetLocations(), apData.GetLocations(), nullptr, 0);
            }

            /**
             * @brief Computes the covariance matrix.
             * @param[in] apDescriptorData pointer to the DescriptorData object holding descriptors and data.
             * @param[out] apDescriptor Pointer to the descriptor for the covariance matrix.
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
            virtual void CovarianceMatrixCodelet(dataunits::DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                 int &aTriangularPart, dataunits::Locations<T> *apLocation1,
                                                 dataunits::Locations<T> *apLocation2,
                                                 dataunits::Locations<T> *apLocation3, T *aLocalTheta,
                                                 int aDistanceMetric, const std::string &aKernelName) = 0;

            /**
             * @brief Copies the descriptor data to a double vector.
             * @param[in] apDescriptorData pointer to the DescriptorData object holding descriptors and data.
             * @param[in] apDescriptor Pointer to the descriptor data.
             * @param[in,out] apDoubleVector Pointer to the double vector to copy the descriptor data to.
             * @return void
             *
             */
            virtual void
            CopyDescriptorZ(dataunits::DescriptorData<T> *apDescriptorData, void *apDescriptor, T *apDoubleVector) = 0;

            /**
             * @brief Generates the observations vector.
             * @param[in] aConfigurations Configurations object containing relevant settings.
             * @param[in] apDescriptorData pointer to the DescriptorData object holding descriptors and data.
             * @param[in] aDescriptor Pointer to the descriptor for the observations vector.
             * @param[in] apLocation1 Pointer to the first set of locations.
             * @param[in] apLocation2 Pointer to the second set of locations.
             * @param[in] apLocation3 Pointer to the third set of locations.
             * @param[in] aDistanceMetric Specifies the distance metric to use.
             * @return void
             *
             */
            virtual void GenerateObservationsVector(configurations::Configurations &aConfigurations,
                                                    dataunits::DescriptorData<T> *apDescriptorData,
                                                    dataunits::BaseDescriptor aDescriptor,
                                                    dataunits::Locations<T> *apLocation1,
                                                    dataunits::Locations<T> *apLocation2,
                                                    dataunits::Locations<T> *apLocation3, int aDistanceMetric) = 0;

            /**
             * @brief Calculates the log likelihood value of a given value theta.
             * @param[in] aHardware  ExaGeoStatHardware object representing the hardware.
             * @param[in] apData MLE_data struct with different MLE inputs.
             * @param[in] aConfigurations Configurations object containing relevant settings.
             * @param[in] apTheta Optimization parameter used by NLOPT.
             * @return log likelihood value
             *
             */
            virtual T
            ExaGeoStatMleTile(const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &apData,
                              configurations::Configurations &apConfigurations, const double *apTheta,
                              T *apMeasurementsMatrix) = 0;

            /**
             * @brief Converts a Gaussian descriptor to a non-tiled descriptor.
             * @param[in] apDescriptorData DescriptorData struct with the Gaussian descriptor.
             * @param[in] apDesc Pointer to the non-tiled descriptor.
             * @param[in] apTheta Theta vector.
             *
             */
            virtual void ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData<T> *apDescriptorData, void *apDesc,
                                                          T *apTheta) = 0;

            /**
             * @brief Copies a matrix in the tile layout from source to destination
             * @param[in] aUpperLower Specifies the part of the matrix A to be copied to B.
             * @param[in] apA Source matrix A.
             * @param[in,out] apB Destination matrix B. On exit, B = A in the locations specified by Upper Lower.
             * @return Successful exit
             *
             */
            virtual int ExaGeoStatLapackCopyTile(common::UpperLower aUpperLower, void *apA, void *apB) = 0;

            /**
             * @brief Conversion from LAPACK layout to Exageostat descriptor.
             * @param[in] aUpperLower Specifies the shape of the matrix A.
             * @param[in] apAf77 LAPACK matrix.
             * @param[in] aLda The leading dimension of the matrix Af77.
             * @param[in] apA Descriptor of the CHAMELEON matrix initialized with data from Af77.
             * @return Successful exit
             *
             */
            virtual int
            ExaGeoStatLapackToDescriptor(common::UpperLower aUpperLower, void *apAf77, int aLda, void *apA) = 0;

            /**
             * @brief Wait for the completion of a sequence.
             * @param[in] apSequence apSequence A pointer to either CHAMELEON or HiCMA sequence.
             * @return successful exit
             *
             */
            virtual int ExaGeoStatSequenceWait(void *apSequence) = 0;

            /**
             * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
             * @param[in] aUpperLower Whether upper or lower part of the matrix A
             * @param[in] apA Symmetric matrix A
             * @return
             *
             */
            virtual int ExaGeoStatPotrfTile(common::UpperLower aUpperLower, void *apA) = 0;

            /**
             * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
             * @param[in] aSide Specifies whether op(A) appears on the left or on the right of X
             * @param[in] aUpperLower Specifies whether the matrix A is upper triangular or lower triangular.
             * @param[in] aTrans Specifies the form of op( A ) to be used in the matrix multiplication.
             * @param[in] aDiag Specifies whether or not A is unit triangular.
             * @param[in] aAlpha Specifies the scalar alpha. When alpha is zero then A is not referenced and B need not be set before entry.
             * @param[in] apA The triangular matrix A
             * @param[in, out] apB The matrix B of dimension ,on exit is overwritten by the solution matrix X.
             * @return successful exit
             *
             */
            virtual int ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans,
                                           common::Diag aDiag, T aAlpha, void *apA, void *apB) = 0;

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
             *
             */
            virtual int
            ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha, void *apA, void *apB, T aBeta,
                               void *apC) = 0;

            /**
             * @brief Calculate determinant for triangular matrix.
             * @param[in] apDescA Exageostat descriptor.
             * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
             * @param[in] apRequest Identifies this function call (for exception handling purposes).
             * @param[in] apDescDet determinant value
             * @return successful exit.
             *
             */
            virtual int
            ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescDet) = 0;

            /**
             * @brief copy Chameleon descriptor to vector float*.
             * @param[in] apDescA Exageostat descriptor A.
             * @param[in] apDescB Exageostat descriptor B.
             * @param[in] apDescC Exageostat descriptor C.
             * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
             * @param[in] apRequest Identifies this function call (for exception handling purposes).
             * @return successful exit.
             *
             */
            virtual int ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apSequence,
                                                       void *apRequest) = 0;

            /**
             * @brief Sets the context.
             * @param[in] apContext The context.
             *
             */
            void SetContext(void *apContext) {
                this->mpContext = apContext;
            }

            void *GetContext() {
                return this->mpContext;
            }

            //// These codlets and structs will be added to another level of abstraction and interface of runtime system. This is a quick fix for now.
            //// TODO: Create a Factory for Runtime system.

            static void cl_dcmg_cpu_func(void *buffers[], void *cl_arg) {
                int m, n, m0, n0;
                exageostat::dataunits::Locations<T> *location1;
                exageostat::dataunits::Locations<T> *location2;
                exageostat::dataunits::Locations<T> *location3;
                T *theta;
                T *A;
                int distance_metric;
                exageostat::kernels::Kernel<T> *kernel;

                A = (T *) STARPU_MATRIX_GET_PTR(buffers[0]);
                starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &location1, &location2, &location3, &theta,
                                           &distance_metric, &kernel);
                kernel->GenerateCovarianceMatrix(A, m, n, m0, n0, *location1, *location2, *location3, theta,
                                                 distance_metric);
            }

            struct starpu_codelet cl_gaussian_to_non =
                    {
                            .where        = STARPU_CPU,
                            .cpu_funcs    = {CORE_gaussian_to_non_starpu},
                            .nbuffers    = 1,
                            .modes        = {STARPU_RW},
                            .name        = "gaussian_to_non"
                    };

            static void CORE_gaussian_to_non_starpu(void *buffers[], void *cl_arg) {
                int m, m0;
                T *z;
                T *theta;
                theta = new T[6];
                z = (T *) STARPU_MATRIX_GET_PTR(buffers[0]);

                starpu_codelet_unpack_args(cl_arg, &m, &m0,
                                           &theta[0], &theta[1], &theta[2],
                                           &theta[3], &theta[4], &theta[5]);

                //core function to convert Z tile from Gaussian to non-Gaussian.
                core_gaussian_to_non(z, theta, m);
                delete[] theta;
            }

            static void core_gaussian_to_non(T *Z, T *localtheta, int m) {

                double xi = localtheta[2];
                double omega = localtheta[3];
                double g = localtheta[4];
                double h = localtheta[5];

                int i;
                if (h < 0) {
                    printf("The kurtosis parameter cannot be negative");
                    return;
                }
                if (g == 0) {
                    for (i = 0; i < m; i++)
                        Z[i] = xi + omega * Z[i] * (exp(0.5 * h * pow(Z[i], 2)));
                } else {
                    for (i = 0; i < m; i++)
                        Z[i] = xi + omega * (exp(g * Z[i]) - 1) * (exp(0.5 * h * pow(Z[i], 2))) / g;
                }
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

            static void CORE_dzcpy_Starpu(void *buffers[], void *cl_arg) {
                int m;
                T *A;
                int m0;
                T *r;

                A = (T *) STARPU_MATRIX_GET_PTR(buffers[0]);
                starpu_codelet_unpack_args(cl_arg, &m, &m0, &r);
                memcpy(A, &r[m0], m * sizeof(T));
            }

            static void CORE_dmdet_starpu(void *buffers[], void *cl_arg) {
                int m;
                int n;
                T *A;
                int m0;
                int n0;
                T det = 0;
                T *determinant = &det;

                *determinant = 0;
                A = (T *) STARPU_MATRIX_GET_PTR(buffers[0]);
                determinant = (T *) STARPU_MATRIX_GET_PTR(buffers[1]);
                starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);
                T local_det = Core_dmdet(A, m);
                *determinant += local_det;
            }

            static T Core_dmdet(T *A, int m) {

                int i;
                T res = 0.0;
                for (i = 0; i < m; i++) {
                    if (A[i + i * m] > 0)
                        res += log(A[i + i * m]);
                }
                return res;
            }

            static void CORE_stride_vecstarpu(void *buffers[], void *cl_arg) {
                int m;
                int tempmm;
                T *A;
                T *B;
                T *C;
                int m0;
                int i;
                int j = 0;

                A = (T *) STARPU_MATRIX_GET_PTR(buffers[0]);
                B = (T *) STARPU_MATRIX_GET_PTR(buffers[1]);
                C = (T *) STARPU_MATRIX_GET_PTR(buffers[2]);
                starpu_codelet_unpack_args(cl_arg, &tempmm, &m0, &m);

                for (i = 0; i < tempmm - 1; i += 2) {
                    B[j] = A[i];
                    C[j] = A[i + 1];
                    j++;
                }
            }

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
                            .cpu_funcs      = {CORE_stride_vecstarpu},
                            .nbuffers       = 3,
                            .modes          = {STARPU_R, STARPU_W, STARPU_W},
                            .name           = "stride_vec"
                    };

            bool recover(char *path, int iter_count, T *theta, T *loglik, int num_params) {

                FILE *fp;
                char *line = nullptr;
                size_t len = 0;
                int count;
                int i;
                char *pch;
                fp = fopen(path, "r");
                if (fp == nullptr) {
                    printf("cannot open observations file\n");
                    exit(EXIT_FAILURE);
                }
                while (getline(&line, &len, fp) != -1) {
                    pch = strtok(line, " ");
                    count = (int) strtol(pch, nullptr, 10);
                    if (count == iter_count) {
                        pch = strtok(nullptr, " ");
                        for (i = 0; i < num_params; i++) {
                            theta[i] = strtol(pch, nullptr, 10);
                            pch = strtok(nullptr, " ");
                        }
                        *loglik = strtol(pch, nullptr, 10);
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

    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP