
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LocationGenerator.cpp
 * @brief Generates and manages spatial locations for ExaGeoStat.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <cmath>
#include <algorithm>

#include <data-generators/LocationGenerator.hpp>
#include <helpers/ByteHandler.hpp>

using namespace exageostat::generators;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T> void LocationGenerator<T>::GenerateLocations(const int &aN, const int &aTimeSlot, const Dimension &aDimension, Locations <T> &aLocations) {

    aLocations.SetSize(aN);
    int index = 0;
    aLocations.SetDimension(aDimension);

    int rootN;
    if (aDimension == Dimension3D) {
        //Cubic root.
        rootN = ceil(cbrt(aN));
    } else {
        //Square root.
        rootN = ceil(sqrt(aN));
    }

    int *grid = new int[rootN]();
    for (auto i = 0; i < rootN; i++) {
        grid[i] = i + 1;
    }
    T range_low = -0.4, range_high = 0.4;

    for (auto i = 0; i < rootN && index < aN; i++) {
        for (auto j = 0; j < rootN && index < aN; j++) {
            if (aDimension == Dimension3D) {
                for (auto k = 0; k < rootN && index < aN; k++) {
                    aLocations.GetLocationX()[index] =
                            (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    aLocations.GetLocationY()[index] =
                            (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    aLocations.GetLocationZ()[index] =
                            (grid[k] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    index++;
                }
            } else {
                aLocations.GetLocationX()[index] =
                        (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                aLocations.GetLocationY()[index] =
                        (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                if (aDimension == DimensionST) {
                    aLocations.GetLocationZ()[index] = 1.0;
                }
                index++;
            }
        }
    }
    delete[] grid;
    if (aDimension != DimensionST) {
        SortLocations(aN, aDimension, aLocations);
    } else {
        for (auto i = 0; i < aN; i++) {
            aLocations.GetLocationX()[i] = aLocations.GetLocationX()[i];
            aLocations.GetLocationY()[i] = aLocations.GetLocationY()[i];
            aLocations.GetLocationZ()[i] = (T) (i / aTimeSlot + 1);
        }
    }
}

template<typename T>
T LocationGenerator<T>::UniformDistribution(const T &aRangeLow, const T &aRangeHigh) {
    T myRand = (T) rand() / (T) (1.0 + RAND_MAX);
    T range = aRangeHigh - aRangeLow;
    return (myRand * range) + aRangeLow;
}

template<typename T>
void
LocationGenerator<T>::SortLocations(const int &aN, const Dimension &aDimension, Locations <T> &aLocations) {

    // Some sorting, required by spatial statistics code
    uint16_t x, y, z;
    uint64_t vectorZ[aN];

    // Encode data into vector z
    for (auto i = 0; i < aN; i++) {
        x = (uint16_t)(aLocations.GetLocationX()[i] * (double) UINT16_MAX + .5);
        y = (uint16_t)(aLocations.GetLocationY()[i] * (double) UINT16_MAX + .5);
        if (aDimension != Dimension2D) {
            z = (uint16_t)(aLocations.GetLocationZ()[i] * (double) UINT16_MAX + .5);
        } else {
            z = (uint16_t) 0.0;
        }
        vectorZ[i] = (SpreadBits(z) << 2) + (SpreadBits(y) << 1) + SpreadBits(x);
    }
    // Sort vector z
    std::sort(vectorZ, vectorZ + aN, CompareUint64);

    // Decode data from vector z
    for (auto i = 0; i < aN; i++) {
        x = ReverseSpreadBits(vectorZ[i] >> 0);
        y = ReverseSpreadBits(vectorZ[i] >> 1);
        z = ReverseSpreadBits(vectorZ[i] >> 2);
        aLocations.GetLocationX()[i] = (double) x / (double) UINT16_MAX;
        aLocations.GetLocationY()[i] = (double) y / (double) UINT16_MAX;
        if (aDimension == Dimension3D) {
            aLocations.GetLocationZ()[i] = (double) z / (double) UINT16_MAX;
        }
    }
}