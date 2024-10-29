
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecHeader.hpp
 * @brief A header file to include hicma_parsec and undo str definition inside it.
 * @details Due to an error occuers when using hicma_parsec.h which set a defination called str() it conflicts with the standart C++ str() fn
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-10-08
**/

#include <hicma_parsec.h>
#undef str
