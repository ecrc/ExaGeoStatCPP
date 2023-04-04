
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file MatrixAllocation.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-27
**/

#ifndef EXAGEOSTATCPP_MATRIXALLOCATION_HPP
#define EXAGEOSTATCPP_MATRIXALLOCATION_HPP

#ifdef EXAGEOSTAT_USE_CHAMELEON
extern "C"{
#include <chameleon.h>
}
#endif

#ifdef EXAGEOSTAT_USE_HiCMA
extern "C"{
#include <hicma.h>
}
#endif


#define EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(_desc_, _isOOC_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
    if (_isOOC_ && _memspace_ == nullptr && _mb_ != 1  && _nb_ !=1)                                                                               \
CHAMELEON_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
        _p_, _q_);             \
else                                                            \
CHAMELEON_Desc_Create(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
        _p_, _q_);                                                                                                                        \


#define EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(_desc_, _isOOC_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
    if (_isOOC_ && _memspace_ == nullptr && _mb_ != 1  && _nb_ !=1)                                                                               \
HICMA_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
        _p_, _q_);             \
else                                                            \
HICMA_Desc_Create(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
        _p_, _q_);                                                                                                                        \


#endif //EXAGEOSTATCPP_MATRIXALLOCATION_HPP
