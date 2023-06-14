
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
using namespace exageostat::common;

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

Locations::Locations(int aSize, Dimension aDimension) {

    this->mpLocationX = (double *) malloc(aSize * sizeof(double));
    this->mpLocationY = (double *) malloc(aSize * sizeof(double));
    if (aDimension != common::Dimension2D) {
        this->mpLocationZ = (double *) malloc(aSize * sizeof(double));
    }
}

Locations::~Locations() {

    if (this->mpLocationX != nullptr) {
        free(this->mpLocationX);
    }
    if (this->mpLocationY != nullptr) {
        free(this->mpLocationY);
    }
    if (this->mpLocationZ != nullptr) {
        free(this->mpLocationZ);
    }
}
