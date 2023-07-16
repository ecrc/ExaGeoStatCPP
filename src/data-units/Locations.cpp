
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
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
    if (this->mpLocationX == nullptr) {
        throw std::runtime_error("LocationX is null");
    }
    return this->mpLocationX;
}

void Locations::SetLocationY(double *apLocationY) {
    this->mpLocationY = apLocationY;
}

double *Locations::GetLocationY() {
    if (this->mpLocationY == nullptr) {
        throw std::runtime_error("LocationY is null");
    }
    return this->mpLocationY;
}

void Locations::SetLocationZ(double *apLocationZ) {
    this->mpLocationZ = apLocationZ;
}

double *Locations::GetLocationZ() {
    return this->mpLocationZ;
}

void Locations::SetSize(int aSize){
    this->mSize = aSize;
}

int Locations::GetSize(){
    return this->mSize;
}

void Locations::SetDimension(Dimension aDimension){
    this->mDimension = aDimension;
}

Dimension Locations::GetDimension(){
    return this->mDimension;
}

Locations * Locations::CalculateMedianLocations() {
    auto median_Locations = new Locations(this->mSize, this->mDimension);

    double x_min = this->GetLocationX()[0], x_max = this->GetLocationX()[0],
            y_min = this->GetLocationY()[0],y_max = this->GetLocationY()[0],
            z_min ,z_max;
    if(this->mDimension != common::Dimension2D){
        z_min = this->GetLocationZ()[0];
        z_max = this->GetLocationZ()[0];
    }

    for(int i = 0; i < this->mSize; i ++){
        double x = this->mpLocationX[i];
        double y = this->mpLocationX[i];

        x_min = (x < x_min) ? x : x_min;
        x_max = (x > x_max) ? x : x_max;
        y_min = (y < y_min) ? y : y_min;
        y_max = (y > y_max) ? y : y_max;

        if(this->mDimension != common::Dimension2D){
            double z = this->mpLocationX[i];
            z_min = (z < z_min) ? z : z_min;
            z_max = (z > z_max) ? z : z_max;
        }
    }
    median_Locations->GetLocationX()[0] = x_min + (x_max - x_min) / 2;
    median_Locations->GetLocationY()[0] = y_min + (y_max - y_min) / 2;
    if (this->mDimension != common::Dimension2D){
        median_Locations->GetLocationY()[0] = z_min + (z_max - z_min) / 2;
    }
    return median_Locations;
}


Locations::Locations(int aSize, Dimension aDimension) {
    this->mSize = aSize;
    this->mDimension = aDimension;
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
