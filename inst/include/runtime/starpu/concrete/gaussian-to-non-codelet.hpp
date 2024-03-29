
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file gaussian-to-non-codelet.hpp
 * @brief A class for starpu codelet gaussian-to-non.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_GAUSSIAN_TO_NON_CODELET_H
#define EXAGEOSTATCPP_GAUSSIAN_TO_NON_CODELET_H

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
    * @class Gaussian-To-Non Codelet
    * @brief A class for starpu codelet gaussian-to-non.
    * @tparam T Data Type: float or double
    * @details This class encapsulates the struct cl_gaussian_to_non and its CPU functions.
    *
    */
    template<typename T>
    class GaussianCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        GaussianCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~GaussianCodelet() = default;

        /**
         * @brief Inserts a task for Gaussian to non-Gaussian conversion codelet processing.
         * @param[in,out] apDesc A pointer to the descriptor for the matrix tile.
         * @param[in] apTheta A pointer to the transformation parameters.
         * @return void
         *
         */
        void InsertTask(void *apDesc, T *apTheta);

    private:

        /**
         * @brief Executes the Gaussian to non-Gaussian conversion codelet function.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the matrix size,
         * offset, and the transformation parameters.
         * @return void
         *
         */
        static void cl_gaussian_to_non_function(void **apBuffers, void *apCodeletArguments);

        /**
         * @brief Transforms data from a Gaussian distribution to a non-Gaussian distribution.
         * @param[in,out] apDescriptorZ A pointer to the array of data to be transformed. This array is modified in place.
         * @param[in] apLocalTheta A pointer to the array of transformation parameters. The first element is the mean (`xi`),
         * the second element is the scale (`omega`), the third element is the skewness (`g`), and the fourth element is the kurtosis (`h`).
         * @param[in] aSize The size of the data array (`pZ`) and the transformation parameters array (`apLocalTheta`).
         * @throws std::runtime_error If the kurtosis parameter (`h`) is negative, indicating an invalid transformation parameter.
         * @return void
         *
         */
        static void core_gaussian_to_non(T *apDescriptorZ, const T *apLocalTheta, const int &aSize);

        /// starpu_codelet struct
        static struct starpu_codelet cl_gaussian_to_non;
    };

    /**
     * @brief Instantiates the gaussian-to-non codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(GaussianCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_GAUSSIAN_TO_NON_CODELET_H
