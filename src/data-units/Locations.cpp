
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

template<typename T>
void Locations<T>::SetLocationX(T *apLocationX) {
    this->mpLocationX = apLocationX;
}

template<typename T>
T *Locations<T>::GetLocationX() {
    if (this->mpLocationX == nullptr) {
        throw std::runtime_error("LocationX is null");
    }
    return this->mpLocationX;
}

template<typename T>
void Locations<T>::SetLocationY(T *apLocationY) {
    this->mpLocationY = apLocationY;
}

template<typename T>
T *Locations<T>::GetLocationY() {
    if (this->mpLocationY == nullptr) {
        throw std::runtime_error("LocationY is null");
    }
    return this->mpLocationY;
}

template<typename T>
void Locations<T>::SetLocationZ(T *apLocationZ) {
    this->mpLocationZ = apLocationZ;
}

template<typename T>
T *Locations<T>::GetLocationZ() {
    return this->mpLocationZ;
}

template<typename T>
void Locations<T>::SetSize(int aSize) {
    this->mSize = aSize;
}

template<typename T>
int Locations<T>::GetSize() {
    return this->mSize;
}

template<typename T>
void Locations<T>::SetDimension(Dimension aDimension) {
    this->mDimension = aDimension;
}

template<typename T>
Dimension Locations<T>::GetDimension() {
    return this->mDimension;
}



template<typename T>
Locations<T>::Locations(int aSize, Dimension aDimension) {

    this->mSize = aSize;
    this->mDimension = aDimension;
    this->mpLocationX = (T *) malloc(aSize * sizeof(T));
    this->mpLocationY = (T *) malloc(aSize * sizeof(T));
    if (aDimension != common::Dimension2D) {
        this->mpLocationZ = (T *) malloc(aSize * sizeof(T));
    }
}

template<typename T>
Locations<T>::~Locations() {


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
