
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file non-gaussian-transform-codelet.hpp
 * @brief A class for starpu codelet non-gaussian-transform.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-26
**/

#ifndef EXAGEOSTATCPP_NON_GAUSSIAN_TRANSFORM_CODELET_HPP
#define EXAGEOSTATCPP_NON_GAUSSIAN_TRANSFORM_CODELET_HPP

#include <runtime/starpu/helpers/StarPuHelpers.hpp>

namespace exageostat::runtime {

    /**
     * @class NonGaussianTransform Codelet
     * @brief A class for starpu codelet non gaussian transform.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_non_gaussian_transform and its CPU functions.
     *
     */
    template<typename T>
    class NonGaussianTransform {

    public:

        /**
         * @brief Default constructor
         *
         */
        NonGaussianTransform() = default;

        /**
         * @brief Default destructor
         *
         */
        ~NonGaussianTransform() = default;

        /**
         * @brief Inserts a task for Non-Gaussian transformation codelet processing.
         * @param[in,out] apDescZ A pointer to the descriptor for the dataset.
         * @param[in] apTheta A pointer to the transformation parameters.
         * @param[in] apStarPuHelpers A reference to a unique pointer of StarPuHelpers, used for accessing and managing data.
         * @return void
         *
         */
        void InsertTask(void *apDescZ, const T *apTheta, std::unique_ptr<StarPuHelpers> &apStarPuHelpers);

    private:

        /**
         * @brief Executes the Non-Gaussian transformation codelet function.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the dataset size,
         * offset, and the transformation parameters.
         * @return void
         *
         */
        static void cl_non_gaussian_transform_function(void **apBuffers, void *apCodeletArguments);

        /**
         * @brief Helper function for transforming a dataset to a non-Gaussian representation. It applies
         * the Non-Gaussian transformation to each element of a dataset using the
         * Newton-Raphson method for finding the root of a function.
         * @param[in,out] apDescripZ A pointer to the dataset to be transformed.
         * @param[in] apLocalTheta A pointer to the transformation parameters.
         * @param[in] aSize The size of the dataset.
         * @return void
         *
         */
        static void core_non_gaussian_transform_helper(T *apDescripZ, const T *apLocalTheta, const int &aSize);

        /**
         * @brief Implements the Newton-Raphson method for finding the root of a function.
         * @param[in] apDescriptorZ The initial guess for the root.
         * @param[in] aTransLocation The location parameter of the transformation.
         * @param[in] aTransScale The scale parameter of the transformation.
         * @param[in] aTransShape The shape parameter of the transformation.
         * @param[in] aTransKurtosis The kurtosis parameter of the transformation.
         * @param[in] aEpsilon The error threshold for the root-finding process.
         * @return T The calculated root of the function.
         *
         */
        static double newton_raphson(T apDescriptorZ, T aTransLocation, T aTransScale, T aTransShape, T aTransKurtosis, T aEpsilon);

        /**
         * @brief Calculates the Non-Gaussian transformation of a value.
         * @param[in] aOriginalValue The original value to be transformed.
         * @param[in] aCurrentValue The current value of the transformation variable.
         * @param[in] aTransLocation The location parameter of the transformation.
         * @param[in] aTransScale The scale parameter of the transformation.
         * @param[in] aTransShape The shape parameter of the transformation.
         * @param[in] aTransKurtosis The kurtosis parameter of the transformation.
         * @return T The transformed value.
         *
         */
        static double tukeyGHTransfor(T aOriginalValue, T aCurrentValue, T aTransLocation, T aTransScale, T aTransShape, T aTransKurtosis);

        /**
         * @brief Calculates the derivative of the Non-Gaussian transformation.
         * @param[in] aCurrentValue The current value of the transformation variable.
         * @param[in] aTransScale The scale parameter of the transformation.
         * @param[in] aTransShape The shape parameter of the transformation.
         * @param[in] aTransKurtosis The kurtosis parameter of the transformation.
         * @return T The derivative of the transformation at the given value.
         *
         */
        static double tukeyGHDiferencial(T aCurrentValue, T aTransScale, T aTransShape, T aTransKurtosis);

        /// starpu_codelet struct
        static struct starpu_codelet cl_non_gaussian_transform;
    };

    /**
     * @brief Instantiates the non-gaussian-transform codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(NonGaussianTransform)

}//namespace exageostat

#endif //EXAGEOSTATCPP_NON_GAUSSIAN_TRANSFORM_CODELET_HPP
