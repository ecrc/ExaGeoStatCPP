
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaDescriptor.hpp
 * @brief Defines the Hicma Descriptor class for creating matrix descriptors using the HICMA library.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-08-15
**/

#ifndef EXAGEOSTATCPP_HICMADESCRIPTOR_HPP
#define EXAGEOSTATCPP_HICMADESCRIPTOR_HPP

#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>

#include <common/Definitions.hpp>

namespace exageostat::dataunits::descriptor {

    /**
     * @brief HicmaDescriptor is a class for creating matrix descriptors by HICMA library.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class HicmaDescriptor {

    public:
        /**
         * @brief Create a Hicma descriptor for a matrix with the given parameters.
         * @param[in] apDescriptor A pointer to the existing HICMA_desc_t descriptor. The new descriptor will be created based on this descriptor.
         * @param[in] aIsOOC A boolean value indicating whether the matrix is out-of-core or not.
         * @param[in] apMatrix A pointer to the beginning of the matrix.
         * @param[in] aFloatPoint The precision of the matrix.
         * @param[in] aMB The number of rows in a tile.
         * @param[in] aNB The number of columns in a tile.
         * @param[in] aSize The size of the matrix in elements including padding.
         * @param[in] aLM The number of rows of the entire matrix.
         * @param[in] aLN The number of columns of the entire matrix.
         * @param[in] aI The row index to the beginning of the sub-matrix.
         * @param[in] aJ The column index to the beginning of the sub-matrix.
         * @param[in] aM The number of rows of the sub-matrix.
         * @param[in] aN The number of columns of the sub-matrix.
         * @param[in] aP The number of rows of the 2D distribution grid.
         * @param[in] aQ The number of columns of the 2D distribution grid.
         * @param[in] aValidOOC Boolean refer to whether this descriptor can be created with OOC technology or not.
         * @return A pointer to the newly created HICMA_desc_t descriptor.
         *
         */
        static HICMA_desc_t *CreateHicmaDescriptor(void *apDescriptor, const bool &aIsOOC, void *apMatrix,
                                                   const common::FloatPoint &aFloatPoint, const int &aMB,
                                                   const int &aNB, const int &aSize, const int &aLM, const int &aLN,
                                                   const int &aI, const int &aJ, const int &aM, const int &aN,
                                                   const int &aP, const int &aQ, const bool &aValidOOC);

        /**
         * @brief destroys and finalize a descriptor
         * @param[in] apDescriptor A pointer to the existing HICMA_desc_t descriptor.
         * @return An error code or success code.
         *
         */
        static int DestroyHicmaDescriptor(void *apDescriptor);
    };

    /**
    * @brief Instantiates the Hicma descriptor methods class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(HicmaDescriptor)

}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMADESCRIPTOR_HPP
