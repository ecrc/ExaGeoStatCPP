
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Kernel.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-12
**/
#include <cmath>
#include <iostream>

extern "C" {
#include <gsl/gsl_sf_bessel.h>
}

#include <kernels/Kernel.hpp>

using namespace exageostat::dataunits;
using namespace exageostat::kernels;

using namespace std;

template<typename T>
T Kernel<T>::CalculateDistance(Locations<T> *apLocations1, Locations<T> *apLocations2, int &aIdxLocation1,
                                    int &aIdxLocation2, int &aDistanceMetric, int &aFlagZ) {

    T x1 = apLocations1->GetLocationX()[aIdxLocation1];
    T y1 = apLocations2->GetLocationY()[aIdxLocation1];
    T x2 = apLocations2->GetLocationX()[aIdxLocation2];
    T y2 = apLocations2->GetLocationY()[aIdxLocation2];

    T dx = x2 - x1;
    T dy = y2 - y1;
    T dz = 0;

    //TODO: I guess the condition should be if aFlagZ == 0 not 1, since in the C version documentation:
    // * @param z_flag  0--->2D or 1-->3D
    if (apLocations1->GetLocationZ() == NULL || apLocations2->GetLocationZ() == NULL || aFlagZ == 1) {
        //if 2D
        if (aDistanceMetric == 1){
            return DistanceEarth(x1, y1, x2, y2);
        }
        return sqrt(pow(dx, 2) + pow(dy, 2));
    } else {
        //if 3D
        if (aDistanceMetric == 1){
            printf("Great Circle (GC) distance is only valid for 2d\n");
            exit(0);
        }
        T z1 =  apLocations1->GetLocationZ()[aIdxLocation1];
        T z2 =  apLocations2->GetLocationZ()[aIdxLocation2];
        dz = z2 - z1;
        return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
    }

}

static double deg2rad(double deg) {
    return (deg * PI / 180);
}

template<typename T>
T Kernel<T>::DistanceEarth(T &aLatitude1, T &aLongitude1, T &aLatitude2, T &aLongitude2) {

//    const T deg2rad = M_PI / 180.0;
//
//    T latitude1 = aLatitude1 * deg2rad;
//    T longitude1 = aLongitude1 * deg2rad;
//    T latitude2 = aLatitude2 * deg2rad;
//    T longitude2 = aLongitude2 * deg2rad;
//
//    T dLat = latitude2 - latitude1;
//    T dLon = longitude2 - longitude1;
//    T a = sin(dLat / 2) * sin(dLat / 2) + cos(latitude1) * cos(latitude2) * sin(dLon / 2) * sin(dLon / 2);
//    T c = 2 * atan2(sqrt(a), sqrt(1 - a));

    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(aLatitude1);
    lon1r = deg2rad(aLongitude1);
    lat2r = deg2rad(aLatitude2);
    lon2r = deg2rad(aLongitude2);
    u = sin((lat2r - lat1r) / 2);
    v = sin((lon2r - lon1r) / 2);
    return 2.0 * EARTH_RADIUS * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

template<typename T>
T Kernel<T>::CalculateDerivativeBesselInputNu(const T &aOrder, const T &aInputValue) {
    if (aOrder < 1) {
        T nu_new = abs(aOrder - 1);
        return (-0.5 * (-CalculateDerivativeBesselNu(nu_new, aInputValue) +
                        CalculateDerivativeBesselNu(abs(aOrder + 1), aInputValue)));
    } else {
        return (-0.5 * (CalculateDerivativeBesselNu(aOrder - 1, aInputValue) +
                        CalculateDerivativeBesselNu(abs(aOrder + 1), aInputValue)));
    }
}

template<typename T>
T Kernel<T>::CalculateDerivativeBesselNu(const T &aOrder, const T &aInputValue) {
    if (aOrder == 0) {
        return 0;
    } else {
        // Use a small step size to calculate the derivative numerically
        const T step_size = 0.000000001;
        return (gsl_sf_bessel_Knu(aOrder + step_size, aInputValue) - gsl_sf_bessel_Knu(aOrder, aInputValue)) /
               step_size;
    }
}

template<typename T>
T Kernel<T>::CalculateSecondDerivativeBesselNu(const T &aOrder, const T &aInputValue) {
    return (-0.5 * (CalculateSecondDerivativeBesselNuInput(aOrder - 1, aInputValue) +
                    CalculateSecondDerivativeBesselNuInput(aOrder + 1, aInputValue)));
}

template<typename T>
T Kernel<T>::CalculateSecondDerivativeBesselNuInput(const T &aOrder, const T &aInputValue) {
    return (aOrder / aInputValue * gsl_sf_bessel_Knu(aOrder, aInputValue) - gsl_sf_bessel_Knu(aOrder + 1, aInputValue));
}

template<typename T>
int Kernel<T>::GetPValue() const {
    return this->mP;
}

template<typename T>
void Kernel<T>::SetPValue(int aP) {
    // Each kernel has its own initial P value, But in case of used spacetime kernels then aP won't be equal to 1.
    // In case of uni-variate spacetime P = 1 * time slot
    // In case of Bi-variate spacetime P = 2 * time slot
    this->mP = this->mP * aP;
}

template<typename T>
int Kernel<T>::GetParametersNumbers() const {
    return this->mParametersNumber;
}

