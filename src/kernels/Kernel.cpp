
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Kernel.cpp
 * @brief implementation file for the Kernels class, which contains the main kernel functions.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-12
**/

extern "C" {
#include <gsl/gsl_sf_bessel.h>
}

#include <kernels/Kernel.hpp>

using namespace std;

using namespace exageostat::dataunits;
using namespace exageostat::kernels;

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
int Kernel<T>::GetP() const {
    return this->mCalculatedP;
}

template<typename T>
void Kernel<T>::SetPValue(int aTimeSlot) {
    // Each kernel has its own initial P value, But in case of used spacetime kernels then Time Slot won't be equal to 1.
    // In case of uni-variate spacetime P = 1 * time slot
    // In case of Bi-variate spacetime P = 2 * time slot
    // P and timeslot will be constant with each created kernel, overriding the Calculated P value to handle case of calling the function multiple times.
    this->mCalculatedP = this->mP * aTimeSlot;
}

template<typename T>
int Kernel<T>::GetParametersNumbers() const {
    return this->mParametersNumber;
}