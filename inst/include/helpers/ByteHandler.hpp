
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

#ifndef EXAGEOSTATCPP_BYTEHANDLER_HPP
#define EXAGEOSTATCPP_BYTEHANDLER_HPP

#include <memory>

namespace exageostat::helpers {

    /**
     * @brief Spread bits by three spaces.
     * @param[in] aInputByte The input 64 bit to be spread.
     * @return The byte after being spread.
     *
     */
    uint64_t SpreadBits(uint64_t aInputByte);

    /**
     * @brief Reverse Spread bits operation.
     * @param[in] aInputByte The input spread 64 bit to be compacted.
     * @return The byte after being compacted.
     *
     */
    uint64_t ReverseSpreadBits(uint64_t aInputByte);

    /**
     * @brief Compares two Unit64 values
     * @param[in] aFirstValue Constant reference to the first input 64 bit value.
     * @param[in] aSecondValue Constant reference to the second input 64 bit value.
     * @return True if the second value is bigger than the first value, false otherwise.
     *
     */
    bool CompareUint64(const uint64_t &aFirstValue, const uint64_t &aSecondValue);

}
#endif //EXAGEOSTATCPP_BYTEHANDLER_HPP
