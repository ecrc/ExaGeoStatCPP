
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dcmg-codelet.hpp
 * @brief A class for starpu codelet dcmg.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-19
**/

#ifndef EXAGEOSTATCPP_DCMG_CODELET_HPP
#define EXAGEOSTATCPP_DCMG_CODELET_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::runtime {

    /**
     * @class DCMG Codelet
     * @brief A class for starpu codelet dcmg.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dcmg and its CPU functions.
     *
     */
    template<typename T>
    class DCMGCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DCMGCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DCMGCodelet() = default;

        /**
         * @brief Inserts a task for DCMG codelet processing.
         * @param[in,out] apDescriptor A pointer to the descriptor containing task information.
         * @param[in] aTriangularPart An integer specifying the triangular part of the matrix (upper or lower).
         * @param[in] apLocation1 A pointer to the first location object for the matrix elements.
         * @param[in] apLocation2 A pointer to the second location object for the matrix elements.
         * @param[in] apLocation3 A pointer to the third location object for the matrix elements.
         * @param[in] apLocalTheta A pointer to the local theta value.
         * @param[in] aDistanceMetric An integer specifying the distance metric to be used.
         * @param[in] apKernel A pointer to the kernel function to be applied during the task execution.
         * @return void
         *
         */
        void InsertTask(void *apDescriptor, const int &aTriangularPart, dataunits::Locations<T> *apLocation1,
                        dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3, T *apLocalTheta,
                        const int &aDistanceMetric, const kernels::Kernel<T> *apKernel);

    private:

        /**
         * @brief CPU Function used by starpu_codelet struct
         * @param[in] apBuffers An array of pointers to the buffers containing the matrix data.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure
         * @return void
         *
         */
        static void cl_dcmg_function(void **apBuffers, void *apCodeletArguments);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dcmg;
    };

    /**
     * @brief Instantiates the dcmg codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DCMGCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DCMG_CODELET_HPP