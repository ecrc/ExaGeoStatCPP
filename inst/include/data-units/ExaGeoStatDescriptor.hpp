/**
 * @file ExaGeoStatDescriptor.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
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
};
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

        /**
         * @Class ExaGeoStatDescriptor
         * @brief Provides functions to create and manipulate descriptors (either chameleon or hicma) for geostatistical matrices.
         * @tparam T
         */
        template<typename T>
        class ExaGeoStatDescriptor {

        public:
#ifdef EXAGEOSTAT_USE_CHAMELEON

            /**
             * @brief Create a chameleon descriptor for a matrix with the given parameters.
             * @param apDescriptor A pointer to the existing CHAM_desc_t descriptor. The new descriptor will be created based on this descriptor.
             * @param aIsOOC A boolean value indicating whether the matrix is out-of-core or not.
             * @param apMatrix A pointer to the beginning of the matrix.
             * @param aFloatPoint The precision of the matrix.
             * @param aMB The number of rows in a tile.
             * @param aNB The number of columns in a tile.
             * @param aSize The size of the matrix in elements including padding.
             * @param aLM The number of rows of the entire matrix.
             * @param aLN The number of columns of the entire matrix.
             * @param aI The row index to the beginning of the submatrix.
             * @param aJ The column index to the beginning of the submatrix.
             * @param aM The number of rows of the submatrix.
             * @param aN The number of columns of the submatrix.
             * @param aP The number of rows of the 2D distribution grid.
             * @param aQ The number of columns of the 2D distribution grid.
             * @return A pointer to the newly created CHAM_desc_t descriptor.
             *
             */
            CHAM_desc_t *CreateChameleonDescriptor(CHAM_desc_t *apDescriptor, bool aIsOOC, void *apMatrix,
                                                   common::FloatPoint aFloatPoint, int aMB, int aNB, int aSize, int aLM,
                                                   int aLN, int aI, int aJ, int aM, int aN, int aP, int aQ);

            /**
             * @brief Create a submatrix for a chameleon descriptor with the given parameters.
             * @param apDescriptor The descriptor for which the submatrix will be created.
             * @param aI The row index to the beginning of the submatrix.
             * @param aJ The column index to the beginning of the submatrix.
             * @param aM The number of rows of the submatrix.
             * @param aN The number of columns of the submatrix.
             * @return A pointer to the newly created CHAM_desc_t submatrix.
             */
            CHAM_desc_t *CreateChameleonSubMatrixDescriptor(CHAM_desc_t *apDescriptor, int aI, int aJ, int aM, int aN);

#endif
#ifdef EXAGEOSTAT_USE_HiCMA

            /**
             * @brief Create a hicma descriptor for a matrix with the given parameters.
             * @param apDescriptor A pointer to the existing CHAM_desc_t descriptor. The new descriptor will be created based on this descriptor.
             * @param aIsOOC A boolean value indicating whether the matrix is out-of-core or not.
             * @param apMatrix A pointer to the beginning of the matrix.
             * @param aFloatPoint The precision of the matrix.
             * @param aMB The number of rows in a tile.
             * @param aNB The number of columns in a tile.
             * @param aSize The size of the matrix in elements including padding.
             * @param aLM The number of rows of the entire matrix.
             * @param aLN The number of columns of the entire matrix.
             * @param aI The row index to the beginning of the submatrix.
             * @param aJ The column index to the beginning of the submatrix.
             * @param aM The number of rows of the submatrix.
             * @param aN The number of columns of the submatrix.
             * @param aP The number of rows of the 2D distribution grid.
             * @param aQ The number of columns of the 2D distribution grid.
             * @return A pointer to the newly created HICMA_desc_t descriptor.
             *
             */
            HICMA_desc_t* CreateHicmaDescriptor(HICMA_desc_t *apDescriptor, bool aIsOOC, void *apMatrix,
                                                common::FloatPoint aFloatPoint, int aMB, int aNB, int aSize, int aLM,
                                                int aLN, int aI, int aJ, int aM, int aN, int aP, int aQ);

            HICMA_desc_t *CreateHicmaSubMatrixDescriptor(HICMA_desc_t *apDescriptor, int aI, int aJ, int aM, int aN);

#endif

        };

        /**
        * @brief Instantiates the Linear Algebra methods class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStatDescriptor)
    }//namespace configurations
}//namespace exageostat


#endif //EXAGEOSTATCPP_EXAGEOSTATDESCRIPTOR_HPP
