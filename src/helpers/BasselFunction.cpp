
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BasselFunction
 * @brief This file contains the BasselFunction class which provides methods for computing derivatives of the modified Bessel function of the second kind. These functions are crucial in statistical and mathematical computations, especially in fields such as geostatistics and spatial analysis.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-01-24
**/

extern "C" {
#include <gsl/gsl_sf_bessel.h>
}

#include <helpers/BasselFunction.hpp>

using namespace exageostat::helpers;

template<typename T>
T BasselFunction<T>::CalculateDerivativeBesselNu(const T &aOrder, const T &aInputValue) {
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
T BasselFunction<T>::CalculateSecondDerivativeBesselNu(const T &aOrder, const T &aInputValue) {
    return (-0.5 * (CalculateSecondDerivativeBesselNuInput(aOrder - 1, aInputValue) +
                    CalculateSecondDerivativeBesselNuInput(aOrder + 1, aInputValue)));
}

template<typename T>
T BasselFunction<T>::CalculateSecondDerivativeBesselNuInput(const T &aOrder, const T &aInputValue) {
    return (aOrder / aInputValue * gsl_sf_bessel_Knu(aOrder, aInputValue) - gsl_sf_bessel_Knu(aOrder + 1, aInputValue));
}
