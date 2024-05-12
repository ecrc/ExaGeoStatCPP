
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DistanceCalculationHelpers.cpp
 * @brief Contains the implementation of the DistanceCalculationHelpers class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <cmath>

#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::helpers;

template<typename T>
T DistanceCalculationHelpers<T>::DegreeToRadian(T aDegree) {
    return (aDegree * PI / 180);
}

template<typename T>
T DistanceCalculationHelpers<T>::DistanceEarth(T &aLatitude1, T &aLongitude1, T &aLatitude2, T &aLongitude2) {
    T lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = DegreeToRadian(aLatitude1);
    lon1r = DegreeToRadian(aLongitude1);
    lat2r = DegreeToRadian(aLatitude2);
    lon2r = DegreeToRadian(aLongitude2);
    u = sin((lat2r - lat1r) / 2);
    v = sin((lon2r - lon1r) / 2);
    return 2.0 * EARTH_RADIUS * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

template<typename T>
T DistanceCalculationHelpers<T>::CalculateDistance(exageostat::dataunits::Locations<T> &aLocations1,
                                                   exageostat::dataunits::Locations<T> &aLocations2,
                                                   const int &aIdxLocation1, const int &aIdxLocation2,
                                                   const int &aDistanceMetric, const int &aFlagZ) {

    T x1 = aLocations1.GetLocationX()[aIdxLocation1];
    T y1 = aLocations1.GetLocationY()[aIdxLocation1];
    T x2 = aLocations2.GetLocationX()[aIdxLocation2];
    T y2 = aLocations2.GetLocationY()[aIdxLocation2];
    T dx = x2 - x1;
    T dy = y2 - y1;
    T dz;

    if (aLocations1.GetLocationZ() == nullptr || aLocations2.GetLocationZ() == nullptr || aFlagZ == 0) {
        //if 2D
        if (aDistanceMetric == 1) {
            return DistanceEarth(x1, y1, x2, y2);
        }
        return sqrt(pow(dx, 2) + pow(dy, 2));
    } else {
        //if 3D
        if (aDistanceMetric == 1) {
            throw std::runtime_error("Great Circle (GC) distance is only valid for 2D!");
        }
        T z1 = aLocations1.GetLocationZ()[aIdxLocation1];
        T z2 = aLocations2.GetLocationZ()[aIdxLocation2];
        dz = z2 - z1;
        return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
    }
}