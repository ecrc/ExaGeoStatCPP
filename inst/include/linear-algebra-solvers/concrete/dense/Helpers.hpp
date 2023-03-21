
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Helpers.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-21
**/

#ifndef EXAGEOSTATCPP_HELPERS_HPP
#define EXAGEOSTATCPP_HELPERS_HPP

#include <chameleon.h>

//// TODO: What's data->ooc ????

#define EXAGEOSTAT_ALLOCATE_MATRIX_TILE(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
    if (_memspace_ == nullptr && _mb_ != 1  && _nb_ !=1)                                                         \
CHAMELEON_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
        _p_, _q_);             \
else                                                            \
CHAMELEON_Desc_Create(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
        _p_, _q_);                                                                                                                        \


#endif //EXAGEOSTATCPP_HELPERS_HPP
