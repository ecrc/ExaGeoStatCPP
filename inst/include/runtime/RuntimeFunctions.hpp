
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RuntimeFunctions.hpp
 * @brief A class for runtime static functions.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-03-10
**/

#ifndef EXAGEOSTATCPP_RUNTIMEFUNCTIONS_HPP
#define EXAGEOSTATCPP_RUNTIMEFUNCTIONS_HPP

#include <kernels/Kernel.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::runtime {

    /**
    * @class RuntimeFunctions
    * @brief A class that defines runtime static functions.
    * @tparam T Data Type: float or double.
    *
    */
    template<typename T>
    class RuntimeFunctions {

    public:

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
        static void
        CovarianceMatrix(dataunits::DescriptorData <T> &aDescriptorData, void *apDescriptor, const int &aTriangularPart,
                         dataunits::Locations <T> *apLocation1, dataunits::Locations <T> *apLocation2,
                         dataunits::Locations <T> *apLocation3, T *apLocalTheta, const int &aDistanceMetric,
                         const kernels::Kernel <T> *apKernel);

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
         * @return void
         *
         */
        static void
        ExaGeoStatMLETileAsyncMLOEMMOM(void *apDescExpr2, void *apDescExpr3, void *apDescExpr4, void *apDescMLOE,
                                       void *apDescMMOM, void *apSequence, void *apRequest);

        /**
        * @brief Calculate mean square prediction error (MSPE) scalar value of the prediction.
        * @param[in] apDescZPredict Observed measurements.
        * @param[in] apDescZMiss Missing measurements.
        * @param[out] apDescError Mean Square Prediction Error (MSPE).
        * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
        * @param[out] apRequest Identifies this function call (for exception handling purposes).
        * @return void
        *
        */
        static void
        ExaGeoStatMLEMSPETileAsync(void *apDescZPredict, void *apDescZMiss, void *apDescError, void *apSequence,
                                   void *apRequest);

        /**
         * @brief Copies the descriptor data to a double vector.
         * @param[in] aComputation computation used in configuration.
         * @param[in] aDescriptorData pointer to the DescriptorData object holding descriptors and data.
         * @param[in] apDescriptor Pointer to the descriptor data.
         * @param[in,out] apDoubleVector Pointer to the double vector to copy the descriptor data to.
         * @return void
         *
         */
        static void
        CopyDescriptorZ(dataunits::DescriptorData <T> &aDescriptorData, void *apDescriptor, T *apDoubleVector);

        /**
        * @brief Converts a Gaussian descriptor to a non-tiled descriptor.
        * @param[in] aDescriptorData DescriptorData struct with the Gaussian descriptor.
        * @param[in] apDesc Pointer to the non-tiled descriptor.
        * @param[in] apTheta Theta vector.
        * @return void
        *
        */
        static void
        ExaGeoStatGaussianToNonTileAsync(dataunits::DescriptorData <T> &aDescriptorData, void *apDesc, T *apTheta);

        /**
        * @brief copy Chameleon descriptor to vector float*.
        * @param[in] apDescA Exageostat descriptor A.
        * @param[in] apDescB Exageostat descriptor B.
        * @param[in] apDescC Exageostat descriptor C.
        * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
        * @param[in] apRequest Identifies this function call (for exception handling purposes).
        * @return void
        *
        */
        static void
        ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apSequence, void *apRequest);

        /**
         * @brief Copy Chameleon descriptor to vector float*.
         * @param[in] apDescA Exageostat descriptor A.
         * @param[in] apDescB Exageostat descriptor B.
         * @param[in] apDescC Exageostat descriptor C.
         * @param[in] apDescD Exageostat descriptor D.
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[in] apRequest Identifies this function call (for exception handling purposes).
         * @return void
         *
         */
        static void
        ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apDescD, void *apSequence,
                                       void *apRequest);

        /**
         * @brief Calculate determinant for triangular matrix.
         * @param[in] aComputation computation used in configuration.
         * @param[in] apDescA Exageostat descriptor.
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[in] apRequest Identifies this function call (for exception handling purposes).
         * @param[in] apDescDet determinant value
         * @return void
         *
         */
        static void
        ExaGeoStatMeasureDetTileAsync(const common::Computation &aComputation, void *apDescA, void *apSequence,
                                      void *apRequest, void *apDescDet);

        /**
        * @brief Calculate determinant for triangular matrix.
        * @param[in] apDescA Pointer to the descriptor of the matrix 'descA'.
        * @param[in] apSequence Pointer to a sequence structure for managing asynchronous execution.
        * @param[in] apRequest Pointer to a request structure for tracking the operation's status.
        * @param[out] apDescNum Pointer to the descriptor of the matrix to store the sum of elements.
        * @param[out] apDescTrace Pointer to the descriptor of the matrix to store the trace.
        * @return void
        *
        */
        static void ExaGeoStatMLETraceTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescNum,
                                                void *apDescTrace);

        /**
        * @brief Computes dot product of A.A.
        * @param[in] apDescA  A Descriptor
        * @param[out] apDescProduct Stores the result of A.A.
        * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
        * @param[in] apRequest Identifies this function call (for exception handling purposes).
        * @return void
        *
        */
        static void ExaGeoStatDoubleDotProduct(void *apDescA, void *apDescProduct, void *apSequence, void *apRequest);

        /**
         * @brief Calculate mean square error (MSE) scalar value for Bivariate kernels.
         * @param[in] apDescZPre Observed measurements descZpre.
         * @param[in] apDescZMiss Missing measurements descZpre
         * @param[out] apDescError1 Mean Square Error (MSE) 1.
         * @param[out] apDescError2 Mean Square Error (MSE) 2.
         * @param[out] apDescError Mean Square Error (MSE).
         * @param[in] apSequence Sequence for the computation.
         * @param[in] apRequest Request for the computation.
         * @return void
         *
         */
        static void
        ExaGeoStatMLEMSPEBivariateTileAsync(void *apDescZPre, void *apDescZMiss, void *apDescError1, void *apDescError2,
                                            void *apDescError, void *apSequence, void *apRequest);

        /**
         * @brief Calculate the log likelihood of non-Gaussian MLE.
         * @param[in] aComputation computation used in configuration.
         * @param[in] apDescZ pointer to the Observed Measurements descriptor.
         * @param[in] apDescSum The log-likelihood Sum of descriptor Z.
         * @param[in] apTheta Pointer to Model parameters.
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[out] apRequest Identifies this function call (for exception handling purposes).
         * @return void
         *
         */
        static void
        ExaGeoStatNonGaussianLogLikeTileAsync(const common::Computation &aComputation, void *apDescZ, void *apDescSum,
                                              const T *apTheta, void *apSequence, void *apRequest);

        /**
        * @brief Transform the measurements vector inside the non-Gaussian MLE function.
        * @param[in] aComputation computation used in configuration.
        * @param[in] apDescZ pointer to the Observed Measurements descriptor.
        * @param[in] apTheta Pointer to Model parameters.
        * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
        * @param[in] apRequest Identifies this function call (for exception handling purposes).
        * @return void
        *
        */
        static void
        ExaGeoStatNonGaussianTransformTileAsync(const common::Computation &aComputation, void *apDescZ,
                                                const T *apTheta, void *apSequence, void *apRequest);
    };

    /**
     * @brief Instantiates the Runtime Functions class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(RuntimeFunctions)

}//namespace exageostat

#endif //EXAGEOSTATCPP_RUNTIMEFUNCTIONS_HPP
