
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ByteHandler.hpp
 * @brief Implementation of byte manipulation functions for ExaGeoStat.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-01-24
**/

#include <helpers/ByteHandler.hpp>

namespace exageostat::helpers {

    uint64_t SpreadBits(uint64_t aInputByte) {

        aInputByte &= 0x000000000000ffff;
        // aInputByte = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
        aInputByte = (aInputByte ^ (aInputByte << 24)) & 0x000000ff000000ff;
        // aInputByte = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
        aInputByte = (aInputByte ^ (aInputByte << 12)) & 0x000f000f000f000f; //000 7000f000f000f
        // aInputByte = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
        aInputByte = (aInputByte ^ (aInputByte << 6)) & 0x0303030303030303; //0 0001 0 0011 0 0011 0 0011 0
        // aInputByte = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
        aInputByte = (aInputByte ^ (aInputByte << 3)) & 0x1111111111111111;
        // aInputByte = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
        return aInputByte;
    }

    uint64_t ReverseSpreadBits(uint64_t aInputByte) {

        aInputByte &= 0x1111111111111111;
        // aInputByte = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
        aInputByte = (aInputByte ^ (aInputByte >> 3)) & 0x0303030303030303;
        // aInputByte = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
        aInputByte = (aInputByte ^ (aInputByte >> 6)) & 0x000f000f000f000f;
        // aInputByte = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
        aInputByte = (aInputByte ^ (aInputByte >> 12)) & 0x000000ff000000ff;
        // aInputByte = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
        aInputByte = (aInputByte ^ (aInputByte >> 24)) & 0x000000000000ffff;
        // aInputByte = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
        return aInputByte;
    }

    bool CompareUint64(const uint64_t &aFirstValue, const uint64_t &aSecondValue) {
        return aFirstValue < aSecondValue;
    }
}