
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmse-codelet.hpp
 * @brief A class for starpu codelet dmse.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-21
**/

#ifndef EXAGEOSTATCPP_DMSE_CODELET_HPP
#define EXAGEOSTATCPP_DMSE_CODELET_HPP

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
     * @class DMSE Codelet
     * @brief A class for starpu codelet dmse.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dmse and its CPU functions.
     *
     */
    template<typename T>
    class DMSECodelet {

    public:

        /**
         * @brief Constructor for DMSE codelet
         *
         */
        DMSECodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DMSECodelet() = default;

        /**
         * @brief Inserts a task for DMSE codelet processing.
         * @param[in,out] apDescError A pointer to the descriptor for the error sum.
         * @param[in] apDescZPredict A pointer to the descriptor for the predicted values.
         * @param[in] apDescZMiss A pointer to the descriptor for the observed values.
         * @return void
         *
         */
        void InsertTask(void *apDescError, void *apDescZPredict, void *apDescZMiss);

    private:

        /**
         * @brief Executes the DMSE codelet function for error calculation.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the vector size and offset.
         * @retur void
         *
         */
        static void cl_dmse_function(void **apBuffers, void *apCodeletArguments);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dmse;
    };

    /**
     * @brief Instantiates the dmse codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DMSECodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DMSE_CODELET_HPP
