
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecDescriptor.cpp
 * @brief Defines the ParsecDescriptor class for creating matrix descriptors using the PaRSEC library.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-10-18
**/

#include <data-units/descriptor/concrete/ParsecDescriptor.hpp>

using namespace exageostat::dataunits::descriptor;

template<typename T>
parsec_matrix_block_cyclic_t *
ParsecDescriptor<T>::CreateParsecDescriptor(void *apDescriptor) {

    parsec_matrix_block_cyclic_t *parsec_desc = new parsec_matrix_block_cyclic_t();
    return parsec_desc;
}

template<typename T>
int ParsecDescriptor<T>::DestroyParsecDescriptor(void *apDesc) {
    auto Parsec_desc = (parsec_matrix_block_cyclic_t *) apDesc;
    parsec_data_free(Parsec_desc->mat);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t *) Parsec_desc);
    return 0;
}