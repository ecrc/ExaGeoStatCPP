
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file MatrixAllocation.hpp
 * @brief Header file for matrix allocation macros.
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

/**
 * @brief Macro for allocating a dense matrix tile.
 *
 * @param[in,out] _desc_ The descriptor for the tile.
 * @param[in] _isOOC_ Whether the matrix is out-of-core.
 * @param[in] _memspace_ The memory space to use for the tile.
 * @param[in] _type2_ The data type of the tile.
 * @param[in] _mb_ The row block size of the tile.
 * @param[in] _nb_ The column block size of the tile.
 * @param[in] _mbXnb_ The product of row and column block sizes.
 * @param[in] _lda_ The leading dimension of the tile.
 * @param[in] _n_ The total number of columns of the matrix.
 * @param[in] _smb_ The row block size for the matrix distribution.
 * @param[in] _snb_ The column block size for the matrix distribution.
 * @param[in] _m_ The total number of rows of the matrix.
 * @param[in] _n2_ The total number of columns of the matrix after padding.
 * @param[in] _p_ The row coordinate of the tile.
 * @param[in] _q_ The column coordinate of the tile.
 */
#define EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(_desc_, _isOOC_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
    if (_isOOC_ && _memspace_ == nullptr && _mb_ != 1  && _nb_ !=1)                                                                               \
        CHAMELEON_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_);             \
    else                                                            \
        CHAMELEON_Desc_Create(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, _p_, _q_);                                                                                                                        \


/**
 * @brief Macro for allocating an approximate matrix tile.
 *
 * @param[in,out] _desc_ The descriptor for the tile.
 * @param[in] _isOOC_ Whether the matrix is out-of-core.
 * @param[in] _memspace_ The memory space to use for the tile.
 * @param[in] _type2_ The data type of the tile.
 * @param[in] _mb_ The row block size of the tile.
 * @param[in] _nb_ The column block size of the tile.
 * @param[in] _mbXnb_ The product of row and column block sizes.
 * @param[in] _lda_ The leading dimension of the tile.
 * @param[in] _n_ The total number of columns of the matrix.
 * @param[in] _smb_ The row block size for the matrix distribution.
 * @param[in] _snb_ The column block size for the matrix distribution.
 * @param[in] _m_ The total number of rows of the matrix.
 * @param[in] _n2_ The total number of columns of the matrix after padding.
 * @param[in] _p_ The row coordinate of the tile.
 * @param[in] _q_ The column coordinate of the tile.
 */
#define EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(_desc_, _isOOC_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
    if (_isOOC_ && _memspace_ == nullptr && _mb_ != 1  && _nb_ !=1)                                                                               \
        HICMA_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
        _p_, _q_);             \
else                                                            \
HICMA_Desc_Create(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
        _p_, _q_);                                                                                                                        \

#endif //EXAGEOSTATCPP_MATRIXALLOCATION_HPP
