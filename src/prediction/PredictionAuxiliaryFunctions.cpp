
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PredictionAuxiliaryFunctions.cpp
 * @brief Contains the implementation of the PredictionAuxiliaryFunctions class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <cmath>

#include <prediction/PredictionAuxiliaryFunctions.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>
#include <utilities/Logger.hpp>

using namespace exageostat::prediction;
using namespace exageostat::dataunits;

template<typename T>
void PredictionAuxiliaryFunctions<T>::PredictIDW(T *apZMiss, T *apZActual, T *apZObs, const int &aZMissNumber,
                                                 const int &aZObsNumber, Locations<T> &aMissLocation,
                                                 Locations<T> &aObsLocation, T *apMSPE) {

    int i, j;
    T sigma_1 = 0;
    T sigma_2 = 0;
    T error = 0, dij;
    T error1 = 0, error2 = 0;
    T x1, y1, x2, y2;
    int n = 2;

    for (j = 0; j < aZMissNumber; j++) {
        for (i = 0; i < aZObsNumber; i++) {
            dij = helpers::DistanceCalculationHelpers<T>::CalculateDistance(aObsLocation, aMissLocation, i, j,
                                                                            common::EUCLIDEAN_DISTANCE, 0);
            if (dij != 0) {
                sigma_1 += apZObs[i] / pow(dij, n);
                sigma_2 += 1.0 / pow(dij, n);
            }
        }
        apZMiss[j] = sigma_1 / sigma_2;
        if (j % 2 == 0)
            error1 += pow((apZActual[j] - apZMiss[j]), 2);
        else
            error2 += pow((apZActual[j] - apZMiss[j]), 2);
        error += pow((apZActual[j] - apZMiss[j]), 2);
        sigma_1 = 0;
        sigma_2 = 0;
    }

    apMSPE[0] = error / aZMissNumber;
    apMSPE[1] = error1 / (aZMissNumber / 2);
    apMSPE[2] = error2 / (aZMissNumber / 2);

    VERBOSE("- Z Actual .. Z Miss")
    for (int index = 0; index < aZMissNumber; index++) {
        VERBOSE(" (" << apZActual[index] << ", " << apZMiss[index] << ")")
    }
    LOGGER("\t\t- Prediction Error (IDW): " << apMSPE[0] << " - " << apMSPE[1] << " - " << apMSPE[2])
}