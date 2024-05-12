
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmse-bivariate-codelet.hpp
 * @brief A class for starpu codelet dmse-bivariate.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_DMSE_BIVARIATE_CODELET_HPP
#define EXAGEOSTATCPP_DMSE_BIVARIATE_CODELET_HPP

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
     * @class DMSE Bivariate Codelet
     * @brief A class for starpu codelet dmse-bivariate.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dmse_bivariate and its CPU functions.
     *
     */
    template<typename T>
    class DMSEBivariateCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DMSEBivariateCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DMSEBivariateCodelet() = default;

        /**
         * @brief Inserts a task for DMSEBivariate codelet processing.
         * @param[in] apDescZMiss A pointer to the descriptor for the observed values.
         * @param[in] apDescZPre A pointer to the descriptor for the predicted values.
         * @param[in,out] apDescsError A pointer to the descriptor for the total error sum.
         * @param[in,out] apDescsError1 A pointer to the descriptor for the error sum for the first variable.
         * @param[in,out] apDescsError2 A pointer to the descriptor for the error sum for the second variable.
         * @return void
         *
         */
        void
        InsertTask(void *apDescZMiss, void *apDescZPre, void *apDescsError, void *apDescsError1, void *apDescsError2);

    private:

        /**
         * @brief Executes the DMSEBivariate codelet function for bivariate error calculation.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the vector size and offset.
         * @return void
         *
         */
        static void cl_dmse_bivariate_function(void **apBuffers, void *apCodeletArguments);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dmse_bivariate;
    };

    /**
     * @brief Instantiates the dmse-bivariate codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DMSEBivariateCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DMSE_BIVARIATE_CODELET_HPP
