
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonDescriptor.hpp
 * @brief Defines the ChameleonDescriptor class for creating matrix descriptors using the CHAMELEON library.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-08-15
**/

#ifndef EXAGEOSTATCPP_CHAMELEONDESCRIPTOR_HPP
#define EXAGEOSTATCPP_CHAMELEONDESCRIPTOR_HPP

extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
}

#include <common/Definitions.hpp>

namespace exageostat {
    namespace dataunits {
        namespace descriptor {

            /**
             * @brief ChameleonDescriptor is a class for creating matrix descriptors by CHAMELEON library.
             * @tparam T Data Type: float or double
             *
             */
            template<typename T>
            class ChameleonDescriptor {

            public:
                /**
                 * @brief Create a chameleon descriptor for a matrix with the given parameters.
                 * @param[in] apDescriptor A pointer to the existing CHAM_desc_t descriptor. The new descriptor will be created based on this descriptor.
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
                 * @return A pointer to the newly created CHAM_desc_t descriptor.
                 *
                 */
                static CHAM_desc_t *CreateChameleonDescriptor(void *apDescriptor, const bool &aIsOOC, void *apMatrix,
                                                              const common::FloatPoint &aFloatPoint, const int &aMB,
                                                              const int &aNB, const int &aSize, const int &aLM,
                                                              const int &aLN, const int &aI, const int &aJ,
                                                              const int &aM, const int &aN, const int &aP,
                                                              const int &aQ);

                /**
                 * @brief destroys and finalize a descriptor
                 * @param[in] apDescriptor A pointer to the existing CHAM_desc_t descriptor.
                 * @return An error code or success code.
                 *
                 */
                static int DestroyChameleonDescriptor(void *apDescriptor);
            };

            /**
            * @brief Instantiates the Chameleon descriptor methods class for float and double types.
            * @tparam T Data Type: float or double
            *
            */
            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonDescriptor)

        }//namespace descriptor
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONDESCRIPTOR_HPP
