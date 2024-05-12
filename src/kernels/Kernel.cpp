
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Kernel.cpp
 * @brief implementation file for the Kernels class, which contains the main kernel functions.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-12
**/

#include <kernels/Kernel.hpp>

using namespace std;

using namespace exageostat::kernels;

template<typename T>
int Kernel<T>::GetVariablesNumber() const {
    return this->mVariablesNumber;
}

template<typename T>
void Kernel<T>::SetPValue(int aTimeSlot) {
    // Each kernel has its own initial P value, But in case of used spacetime kernels then Time Slot won't be equal to 1.
    // In case of uni-variate spacetime P = 1 * time slot
    // In case of Bi-variate spacetime P = 2 * time slot
    // P and timeslot will be constant with each created kernel, overriding the Calculated P value to handle case of calling the function multiple times.
    this->mVariablesNumber = this->mP * aTimeSlot;
}

template<typename T>
int Kernel<T>::GetParametersNumbers() const {
    return this->mParametersNumber;
}