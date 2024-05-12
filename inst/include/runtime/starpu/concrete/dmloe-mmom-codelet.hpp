
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmloe-mmom-codelet.hpp
 * @brief A class for starpu codelet dmloe-mmom.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-19
**/

#ifndef EXAGEOSTATCPP_DMLOE_MMOM_HPP
#define EXAGEOSTATCPP_DMLOE_MMOM_HPP

#include <common/Definitions.hpp>

namespace exageostat::runtime {

    /**
     * @class Dmloe-Mmom Codelet
     * @brief A class for starpu codelet dmloe-mmom.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dmloe_mmom and its CPU functions.
     *
     */
    template<typename T>
    class DmloeMmomCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DmloeMmomCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DmloeMmomCodelet() = default;

        /**
         * @brief Inserts a task for DmloeMmom codelet processing.
         * @param[in] apDescExpr1 A pointer to the descriptor for the first expression.
         * @param[in] apDescExpr2 A pointer to the descriptor for the second expression.
         * @param[in] apDescExpr3 A pointer to the descriptor for the third expression.
         * @param[in,out] apDescMLOE A pointer to the descriptor for the MLOE result.
         * @param[in,out] apDescMMOM A pointer to the descriptor for the MMOM result.
         * @return void
         *
         */
        void InsertTask(void *apDescExpr1, void *apDescExpr2, void *apDescExpr3, void *apDescMLOE,
                        void *apDescMMOM);

    private:

        /**
         * @brief Executes the DmloeMmom codelet function for MLOE and MMOM calculations.
         * @param[in] apBuffers An array of pointers to the buffers.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the matrix dimensions and offsets.
         * @return void
         *
         */
        static void cl_dmloe_mmom_function(void **apBuffers, void *apCodeletArguments);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dmloe_mmom;
    };

    /**
     * @brief Instantiates the dmloe-mmom codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DmloeMmomCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DMLOE_MMOM_HPP
