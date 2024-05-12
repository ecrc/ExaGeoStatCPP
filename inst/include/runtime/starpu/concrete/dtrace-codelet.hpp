
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dtrace-codelet.hpp
 * @brief A class for starpu codelet dtrace.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_DTRACE_CODELET_HPP
#define EXAGEOSTATCPP_DTRACE_CODELET_HPP

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
     * @class DTRACE Codelet
     * @brief A class for starpu codelet dtrace.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dtrace and its CPU functions.
     *
     */
    template<typename T>
    class DTRACECodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DTRACECodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DTRACECodelet() = default;

        /**
         * @brief Inserts a task for DTRACE codelet processing.
         * @param[in] apDescA A pointer to the descriptor for the matrix.
         * @param[in,out] apDescNum A pointer to the descriptor for the sum.
         * @param[in,out] apDescTrace A pointer to the descriptor for the trace.
         * @return void
         *
         */
        void InsertTask(void *apDescA, void *apDescNum, void *apDescTrace);

    private:

        /**
         * @brief Executes the DTRACE codelet function for matrix trace calculation.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the matrix size.
         * @return void
         *
         */
        static void cl_dtrace_function(void **apBuffers, void *apCodeletArguments);

        /**
         * @brief Calculates the trace of a matrix.
         * @param[in] pDescriptor A pointer to the matrix data.
         * @param[in] aSize The size of the matrix (assumed to be square).
         * @param[in,out] pTrace A pointer to the buffer where the trace value will be stored.
         * @return The calculated trace of the matrix.
         *
         */
        static double core_dtrace(const T *pDescriptor, const int &aSize, T *pTrace);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dtrace;

    };

    /**
     * @brief Instantiates the dtrace codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DTRACECodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DTRACE_CODELET_HPP
