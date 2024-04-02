
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ddotp-codelet.hpp
 * @brief A class for starpu codelet ddotp.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_DDOTP_CODELET_HPP
#define EXAGEOSTATCPP_DDOTP_CODELET_HPP

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
     * @class DDOTP Codelet
     * @brief A class for starpu codelet ddotp.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_ddotp and its CPU functions.
     *
     */
    template<typename T>
    class DDOTPCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DDOTPCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DDOTPCodelet() = default;

        /**
         * @brief Inserts a task for DDOTP codelet processing.
         * @param[in] apDescA A pointer to the descriptor for the vector.
         * @param[in,out] apDescProduct A pointer to the descriptor for the dot product.
         * @return void
         *
         */
        void InsertTask(void *apDescA, void *apDescProduct);

    private:

        /**
         * @brief Executes the DDOTP codelet function for dot product calculation.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the vector size (m) and the offset (m0).
         * @return void
         *
         */
        static void cl_ddotp_function(void **apBuffers, void *apCodeletArguments);

        /// starpu_codelet struct
        static struct starpu_codelet cl_ddotp;

    };

    /**
     * @brief Instantiates the ddotp codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DDOTPCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DDOTP_CODELET_HPP
