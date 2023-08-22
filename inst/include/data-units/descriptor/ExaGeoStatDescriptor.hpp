
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatDescriptor.hpp
 * @brief Class for creating matrix descriptors used in CHAMELEON and HiCMA libraries.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-07-16
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP
#define EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP

#include <iostream>

#ifdef EXAGEOSTAT_USE_CHAMELEON
extern "C" {
#include <chameleon/struct.h>
}
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
extern "C"{
#include <hicma_struct.h>
}
#endif

#include <common/Definitions.hpp>

/**
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containning in order
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
namespace exageostat {
    namespace dataunits {
        namespace descriptor {

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
                 * @param[in] aI The row index to the beginning of the submatrix.
                 * @param[in] aJ The column index to the beginning of the submatrix.
                 * @param[in] aM The number of rows of the submatrix.
                 * @param[in] aN The number of columns of the submatrix.
                 * @param[in] aP The number of rows of the 2D distribution grid.
                 * @param[in] aQ The number of columns of the 2D distribution grid.
                 * @return A pointer to the newly created descriptor.
                 *
                 */
                void *CreateDescriptor(void *apDescriptor, common::DescriptorType aDescriptorType, bool aIsOOC,
                                       void *apMatrix, common::FloatPoint aFloatPoint, int aMB, int aNB, int aSize,
                                       int aLM, int aLN, int aI, int aJ, int aM, int aN, int aP, int aQ);

                /**
                 * @brief destroys and finalize a descriptor
                 * @param[in] aDescriptorType The type of the descriptor.
                 * @param[in] apDescriptor A pointer to the existing descriptor.
                 * @return An error code or success code.
                 *
                 */
                int DestroyDescriptor(common::DescriptorType aDescriptorType, void *apDescriptor);
            };

            /**
            * @brief Instantiates the ExaGeoStat Descriptor methods class for float and double types.
            * @tparam T Data Type: float or double
            *
            */
            EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStatDescriptor)

        }//namespace descriptor
    }//namespace configurations
}//namespace exageostat


#endif //EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP
