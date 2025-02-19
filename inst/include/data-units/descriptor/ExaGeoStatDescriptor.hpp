
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatDescriptor.hpp
 * @brief Class for creating matrix descriptors used in CHAMELEON and HiCMA libraries.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-07-16
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP
#define EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP

#if DEFAULT_RUNTIME

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>

#else
#include <runtime/parsec/ParsecHeader.h>
#endif

#include <common/Definitions.hpp>

/**
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containing in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
namespace exageostat::dataunits::descriptor {

    /**
     * @brief ExaGeoStatDescriptor is a class for creating matrix descriptors used in CHAMELEON and HiCMA libraries.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ExaGeoStatDescriptor {

    public:

        /**
         * @brief Create a descriptor for a matrix with the given parameters.
         * @param[out] apDescriptor A pointer to the existing to the descriptor. The new descriptor will be created based on this descriptor Type.
         * @param[in] aDescriptorType The type of the descriptor.
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
         * @return A pointer to the newly created descriptor.
         *
         */
        void *CreateDescriptor(void *apDescriptor, const common::DescriptorType &aDescriptorType, const bool &aIsOOC,
                               void *apMatrix, const common::FloatPoint &aFloatPoint, const int &aMB, const int &aNB,
                               const int &aSize, const int &aLM, const int &aLN, const int &aI, const int &aJ,
                               const int &aM, const int &aN, const int &aP, const int &aQ, const bool &aValidOOC);

        /**
         * @brief destroys and finalize a descriptor
         * @param[in] aDescriptorType The type of the descriptor.
         * @param[in] apDescriptor A pointer to the existing descriptor.
         * @return An error code or success code.
         *
         */
        int DestroyDescriptor(const common::DescriptorType &aDescriptorType, void *apDescriptor);
    };

    /**
    * @brief Instantiates the ExaGeoStat Descriptor methods class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStatDescriptor)

}//namespace exageostat


#endif //EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP
