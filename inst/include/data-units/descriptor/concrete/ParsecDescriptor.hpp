
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecDescriptor.hpp
 * @brief Defines the ParsecDescriptor class for creating matrix descriptors using the PaRSEC library.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-10-18
**/

#ifndef EXAGEOSTATCPP_ParsecDESCRIPTOR_HPP
#define EXAGEOSTATCPP_ParsecDESCRIPTOR_HPP

#include <runtime/parsec/ParsecHeader.h>
#include <common/Definitions.hpp>

namespace exageostat::dataunits::descriptor {

    /**
     * @brief ParsecDescriptor is a class for creating matrix descriptors by Parsec library.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ParsecDescriptor {

    public:
        /**
         * @brief Create a Parsec descriptor for a matrix with the given parameters.
         * @param[in] apDescriptor A pointer to the existing parsec_matrix_block_cyclic_t descriptor. The new descriptor will be created based on this descriptor.
         * @return A pointer to the newly created parsec_matrix_block_cyclic_t descriptor.
         *
         */
        static parsec_matrix_block_cyclic_t *CreateParsecDescriptor(void *apDescriptor);

        /**
         * @brief destroys and finalize a descriptor
         * @param[in] apDescriptor A pointer to the existing parsec_matrix_block_cyclic_t descriptor.
         * @return An error code or success code.
         *
         */
        static int DestroyParsecDescriptor(void *apDescriptor);
    };

    /**
    * @brief Instantiates the Parsec descriptor methods class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(ParsecDescriptor)

}//namespace exageostat

#endif //EXAGEOSTATCPP_ParsecDESCRIPTOR_HPP
