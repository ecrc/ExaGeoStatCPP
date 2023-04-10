
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Locations.cpp
 * @brief Implementation of the Locations class
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-27
**/
#include <data-units/Locations.hpp>

using namespace exageostat::dataunits;

void Locations::SetLocationX(double *apLocationX) {
    this->mpLocationX = apLocationX;
}

double *Locations::GetLocationX() {
    return this->mpLocationX;
}

void Locations::SetLocationY(double *apLocationY) {
    this->mpLocationY = apLocationY;
}

double *Locations::GetLocationY() {
    return this->mpLocationY;
}

void Locations::SetLocationZ(double *apLocationZ) {
    this->mpLocationZ = apLocationZ;
}

double *Locations::GetLocationZ() {
    return this->mpLocationZ;
}