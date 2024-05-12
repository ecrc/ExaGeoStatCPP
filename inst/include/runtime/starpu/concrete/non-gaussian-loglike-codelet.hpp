
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file non-gaussian-loglike-codelet.hpp
 * @brief A class for starpu codelet non-gaussian-loglike.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-26
**/

#ifndef EXAGEOSTATCPP_NON_GAUSSIAN_LOGLIKE_CODELET_HPP
#define EXAGEOSTATCPP_NON_GAUSSIAN_LOGLIKE_CODELET_HPP

#include <runtime/starpu/helpers/StarPuHelpers.hpp>

namespace exageostat::runtime {

    /**
     * @class NonGaussianLoglike
     * @brief A class for starpu codelet non gaussian loglike.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_non_gaussian_loglike and its CPU functions.
     *
     */
    template<typename T>
    class NonGaussianLoglike {

    public:

        /**
         * @brief Constructor for NonGaussianLoglike
         *
         */
        NonGaussianLoglike() = default;

        /**
         * @brief Default destructor
         *
         */
        ~NonGaussianLoglike() = default;

        /**
         * @brief Inserts a task for Non-Gaussian log-likelihood codelet processing.
         * @param[in] apDescZ A pointer to the descriptor for the dataset.
         * @param[in,out] apDescSum A pointer to the descriptor for the sum.
         * @param[in] apTheta A pointer to the transformation parameters.
         * @param[in] aStarPuHelpers A reference to a unique pointer of StarPuHelpers, used for accessing and managing data.
         * @return void
         *
         */
        void
        InsertTask(void *apDescZ, void *apDescSum, const T *apTheta, std::unique_ptr<StarPuHelpers> &aStarPuHelper);

    private:

        /**
         * @brief Executes the Non-Gaussian log-likelihood codelet function.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the dataset size,
         * offset, and the transformation parameters.
         * @return void
         *
         */
        static void cl_non_gaussian_loglike_function(void **apBuffers, void *apCodeletArguments);

        /**
         * @brief Helper function for calculating the log-likelihood under a non-Gaussian distribution.
         * @param[in] apDescriptorZ A pointer to the dataset.
         * @param[in] apLocalTheta A pointer to the transformation parameters.
         * @param[in] aSize The size of the dataset.
         * @return T The calculated log-likelihood of the dataset under a non-Gaussian distribution.
         *
         */
        static double core_non_gaussian_loglike_helper(const T *apDescriptorZ, const T *apLocalTheta, const int &aSize);

        /// starpu_codelet struct
        static struct starpu_codelet cl_non_gaussian_loglike;
    };

    /**
     * @brief Instantiates the non-gaussian-loglike class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(NonGaussianLoglike)

}//namespace exageostat

#endif //EXAGEOSTATCPP_NON_GAUSSIAN_LOGLIKE_CODELET_HPP
