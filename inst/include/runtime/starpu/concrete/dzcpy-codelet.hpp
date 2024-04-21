
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dzcpy-codelet.hpp
 * @brief A class for starpu codelet dzcpy.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_DZCPY_CODELET_HPP
#define EXAGEOSTATCPP_DZCPY_CODELET_HPP

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
     * @class DZCPY Codelet
     * @brief A class for starpu codelet dzcpy.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dzcpy and its CPU functions.
     *
     */
    template<typename T>
    class DZCPYCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DZCPYCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DZCPYCodelet() = default;

        /**
         * @brief Inserts a task for DZCPY codelet processing.
         * @param[in,out] apDescriptor A pointer to the descriptor for the vector.
         * @param[in] apDoubleVector A pointer to the double vector to be copied.
         * @return void
         *
         */
        void InsertTask(void *apDescriptor, void *apDoubleVector);

    private:

        /**
          * @brief Executes the DZCPY codelet function for copying a double vector.
          * @param[in] apBuffers An array of pointers to the buffers.
          * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the vector size,
          * offset, and the pointer to the destination vector.
          * @return void
          *
          */
        static void cl_dzcpy_function(void **apBuffers, void *apCodeletArguments);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dzcpy;
    };

    /**
     * @brief Instantiates the dzcpy codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DZCPYCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DZCPY_CODELET_HPP


