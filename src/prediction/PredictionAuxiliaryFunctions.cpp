
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PredictionAuxiliaryFunctions.cpp
 * @brief Contains the implementation of the PredictionAuxiliaryFunctions class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <cmath>

#include <prediction/PredictionAuxiliaryFunctions.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::common;
using namespace exageostat::prediction;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
void
PredictionAuxiliaryFunctions<T>::PredictTileIDW(T *apZMiss, T *apZActual, T *apZObs, int aZMissNumber, int aZObsNumber,
                                                Locations<T> &aMissLocation, Locations<T> &aObsLocation, T *apMSPE) {
    int i, j;
    T sigma_1 = 0;
    T sigma_2 = 0;
    T error = 0, dij;
    T error1 = 0, error2 = 0;
    T x1, y1, x2, y2;
    int n = 2;

    for (j = 0; j < aZMissNumber; j++) {
        for (i = 0; i < aZObsNumber; i++) {
            dij = DistanceCalculationHelpers<T>::CalculateDistance(aObsLocation, aMissLocation, i, j,
                                                                   DistanceMetric::EUCLIDIAN_DISTANCE, 0);
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

    int index;
    for (index = 0; index < aZMissNumber; index++)
        printf("(%3.6f, %3.6f)\n ", apZActual[index], apZMiss[index]);
    fprintf(stderr, "\n\nPrediction Error (IDW): %3.9f  -  %3.9f -  %3.9f\n", apMSPE[0], apMSPE[1], apMSPE[2]);
}