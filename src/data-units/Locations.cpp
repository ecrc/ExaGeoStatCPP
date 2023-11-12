
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Locations.cpp
 * @brief Implementation of the Locations class
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-27
**/

#include <iostream>
#include <cstring>

#include <data-units/Locations.hpp>

using namespace std;

using namespace exageostat::dataunits;
using namespace exageostat::common;

template<typename T>
void Locations<T>::SetLocationX(T &aLocationX, const int &aSize) {

    if (aLocationX && aSize == this->mSize) {
        memcpy(this->mpLocationX, &aLocationX, this->mSize * sizeof(T));
    } else {
        throw std::runtime_error("Invalid value for setting Locations X");
    }
}

template<typename T>
T *Locations<T>::GetLocationX() {

    if (this->mpLocationX == nullptr) {
        throw std::runtime_error("LocationX is null");
    }
    return this->mpLocationX;
}

template<typename T>
void Locations<T>::SetLocationY(T &aLocationY, const int &aSize) {

    if (aLocationY && aSize == this->mSize) {
        memcpy(this->mpLocationY, &aLocationY, this->mSize * sizeof(T));
    } else {
        throw std::runtime_error("Invalid value for setting Locations Y");
    }
}

template<typename T>
T *Locations<T>::GetLocationY() {
    if (this->mpLocationY == nullptr) {
        throw std::runtime_error("LocationY is null");
    }
    return this->mpLocationY;
}

template<typename T>
void Locations<T>::SetLocationZ(T &aLocationZ, const int &aSize) {

    if (aLocationZ && aSize == this->mSize) {
        memcpy(this->mpLocationZ, &aLocationZ, this->mSize * sizeof(T));
    } else {
        throw std::runtime_error("Invalid value for setting Locations Z");
    }
}

template<typename T>
T *Locations<T>::GetLocationZ() {
    return this->mpLocationZ;
}

template<typename T>
void Locations<T>::SetSize(const int &aSize) {
    this->mSize = aSize;
}

template<typename T>
int Locations<T>::GetSize() {
    return this->mSize;
}

template<typename T>
void Locations<T>::SetDimension(const Dimension &aDimension) {
    this->mDimension = aDimension;
}

template<typename T>
Dimension Locations<T>::GetDimension() {
    return this->mDimension;
}

template<typename T>
Locations<T>::Locations(const int &aSize, const Dimension &aDimension) {

    this->mSize = aSize;
    this->mDimension = aDimension;
    this->mpLocationX = new T[aSize];
    this->mpLocationY = new T[aSize];
    if (aDimension != common::Dimension2D) {
        this->mpLocationZ = new T[aSize];
    }
}

template<typename T>
Locations<T>::~Locations() {

    if (this->mpLocationX != nullptr) {
        delete[] this->mpLocationX;
    }
    if (this->mpLocationY != nullptr) {
        delete[] this->mpLocationY;
    }
    if (this->mpLocationZ != nullptr) {
        delete[] this->mpLocationZ;
    }
}
