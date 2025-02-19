
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file EnumStringParser.hpp
 * @brief Provides utility functions for parsing enumeration values from strings.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-01-20
**/

#ifndef EXAGEOSTATCPP_ENUMSTRINGPARSER_HPP
#define EXAGEOSTATCPP_ENUMSTRINGPARSER_HPP

#ifdef min
#undef min
#endif

#include <algorithm>

#include <utilities/ErrorHandler.hpp>
#include <common/Definitions.hpp>

/**
 * @brief Convert a string representation of computation mode to its corresponding enum value.
 * @param[in] aComputation String representation of computation mode.
 * @return Computation enum value.
 */
inline exageostat::common::Computation GetInputComputation(std::string aComputation) {
    std::transform(aComputation.begin(), aComputation.end(),
                   aComputation.begin(), ::tolower);

    if (aComputation == "exact" || aComputation == "dense") {
        return exageostat::common::EXACT_DENSE;
    } else if (aComputation == "dst" || aComputation == "diag_approx") {
        return exageostat::common::DIAGONAL_APPROX;
    } else if (aComputation == "tlr" || aComputation == "tile_low_rank") {
        return exageostat::common::TILE_LOW_RANK;
    } else {
        const std::string msg = "Error in Initialization : Unknown computation Value" + std::string(aComputation);
        throw API_EXCEPTION(msg, INVALID_ARGUMENT_ERROR);
    }
}

/**
 * @brief Converts string to dimension enum.
 * @param[in] aDimension Dimension as a string.
 * @return Dimension as an enum.
 */
inline exageostat::common::Dimension GetInputDimension(std::string aDimension) {
    std::transform(aDimension.begin(), aDimension.end(),
                   aDimension.begin(), ::tolower);

    if (aDimension == "2d") {
        return exageostat::common::Dimension2D;
    } else if (aDimension == "3d") {
        return exageostat::common::Dimension3D;
    } else if (aDimension == "st") {
        return exageostat::common::DimensionST;
    } else {
        const std::string msg = "Error in Initialization : Unknown computation Value" + std::string(aDimension);
        throw API_EXCEPTION(msg, INVALID_ARGUMENT_ERROR);
    }
}

#endif //EXAGEOSTATCPP_ENUMSTRINGPARSER_HPP
