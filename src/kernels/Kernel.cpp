
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Kernel.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-12
**/
#include <kernels/Kernel.hpp>
#include <cmath>
#include <iostream>

using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace std;

double Kernel::CalculateDistance(Locations *apLocations1, Locations *apLocations2, int aIdxLocation1, int aIdxLocation2, int aDistanceMetric, int aFlagZ) {

    double x1 = apLocations1->GetLocationX()[aIdxLocation1];
    double y1 = apLocations1->GetLocationY()[aIdxLocation1];
    double x2 = apLocations2->GetLocationX()[aIdxLocation2];
    double y2 = apLocations2->GetLocationY()[aIdxLocation2];

    double dx = apLocations1->GetLocationX()[aIdxLocation1] - apLocations2->GetLocationX()[aIdxLocation2];
    double dy = apLocations1->GetLocationY()[aIdxLocation1] - apLocations2->GetLocationY()[aIdxLocation2];
    double dz = 0.0;

    if (aFlagZ == 1) {
        if (aDistanceMetric == 1) {
            throw std::runtime_error("Great Circle (GC) distance is only valid for 2D");
        }
        dz = apLocations1->GetLocationZ()[aIdxLocation1] - apLocations2->GetLocationZ()[aIdxLocation2];
    }
    if (aDistanceMetric == 1){
        return DistanceEarth(x1, y1, x2, y2);
    }
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double Kernel::DistanceEarth(double aLatitude1, double aLongitude1, double aLatitude2, double aLongitude2) {

    const double deg2rad = M_PI / 180.0;

    double latitude1 = aLatitude1 * deg2rad;
    double longitude1 = aLongitude1 * deg2rad;
    double latitude2 = aLatitude2 * deg2rad;
    double longitude2 = aLongitude2 * deg2rad;

    double dLat = latitude2 - latitude1;
    double dLon = longitude2 - longitude1;
    double a = sin(dLat/2) * sin(dLat/2) + cos(latitude1) * cos(latitude2) * sin(dLon/2) * sin(dLon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));

    return EARTH_RADIUS * c;
}


