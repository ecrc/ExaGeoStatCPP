
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file dmdet-codelet.hpp
 * @brief A class for starpu codelet dmdet.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-21
**/

#ifndef EXAGEOSTATCPP_DMDET_CODELET_HPP
#define EXAGEOSTATCPP_DMDET_CODELET_HPP

#include <runtime/starpu/helpers/StarPuHelpers.hpp>

namespace exageostat::runtime {

    /**
     * @class DMDET Codelet
     * @brief A class for starpu codelet dmdet.
     * @tparam T Data Type: float or double
     * @details This class encapsulates the struct cl_dmdet and its CPU functions.
     *
     */
    template<typename T>
    class DMDETCodelet {

    public:

        /**
         * @brief Default constructor
         *
         */
        DMDETCodelet() = default;

        /**
         * @brief Default destructor
         *
         */
        ~DMDETCodelet() = default;

        /**
         * @brief Inserts a task for DMDET codelet processing.
         * @param[in] aComputation The type of computation to be performed, such as diagonal approximation or exact dense computation.
         * @param[in] apDescA A pointer to the descriptor for matrix A.
         * @param[in,out] apDescDet A pointer to the descriptor for the determinant.
         * @param[in] aStarPuHelpers A reference to a unique pointer of StarPuHelpers, used for accessing and managing data.
         * @return void
         *
         */
        void InsertTask(const common::Computation &aComputation, void *apDescA, void *apDescDet,
                        std::unique_ptr<StarPuHelpers> &aStarPuHelpers);

    private:

        /**
         * @brief Executes the DMDET codelet function for matrix determinant calculation.
         * @param[in] apBuffers An array of pointers to the buffers containing the matrix data and the determinant.
         * @param[in] apCodeletArguments A pointer to the codelet arguments structure, which includes the matrix size.
         * @return void
         *
         */
        static void cl_dmdet_function(void **apBuffers, void *apCodeletArguments);

        /**
         * @brief Calculates the determinant of a matrix.
         * @param[in] apDescriptor A pointer to the matrix data.
         * @param[in] aSize The size of the matrix (assumed to be square).
         * @return T The calculated determinant of the matrix.
         *
         */
        static T core_dmdet(const T *apDescriptor, const int &aSize);

        /// starpu_codelet struct
        static struct starpu_codelet cl_dmdet;

    };

    /**
     * @brief Instantiates the dmdet codelet class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DMDETCodelet)

}//namespace exageostat

#endif //EXAGEOSTATCPP_DMDET_CODELET_HPP
