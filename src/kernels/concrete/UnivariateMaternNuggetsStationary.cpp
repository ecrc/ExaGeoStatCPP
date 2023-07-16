
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNuggetsStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternNuggetsStationary.hpp>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternNuggetsStationary::UnivariateMaternNuggetsStationary() {
    this->mP = 1;
    this->mParametersNumber = 4;
}

Kernel *UnivariateMaternNuggetsStationary::Create() {
    return new UnivariateMaternNuggetsStationary();
}

namespace exageostat::kernels {
    bool UnivariateMaternNuggetsStationary::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternNuggetsStationary", UnivariateMaternNuggetsStationary::Create);
}

void UnivariateMaternNuggetsStationary::GenerateCovarianceMatrix(double *apMatrixA, int &aRowsNumber, int &aColumnsNumber,
                                                       int &aRowOffset, int &aColumnOffset, Locations *apLocation1,
                                                       Locations *apLocation2, Locations *apLocation3,
                                                       double *aLocalTheta, int &aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double expr;
    double con;
    double sigma_square = aLocalTheta[0];// * aLocalTheta[0];

    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    con = sigma_square * con;
    int flag = 0;

    if (apLocation1->GetLocationZ() == nullptr || apLocation2->GetLocationZ() == nullptr) {
        for (i = 0; i < aRowsNumber; i++) {
            j0 = aColumnOffset;
            for (j = 0; j < aColumnsNumber; j++) {
                expr = CalculateDistance(apLocation1, apLocation2, j0, i0, aDistanceMetric, flag) / aLocalTheta[1];
                if (expr == 0)
                    apMatrixA[i + j * aRowsNumber] = sigma_square + aLocalTheta[3];
                else
                    apMatrixA[i + j * aRowsNumber] =
                            con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr); // Matern Function
                j0++;
            }
            i0++;
        }
    } else {
        for (i = 0; i < aRowsNumber; i++) {
            j0 = aColumnOffset;
            for (j = 0; j < aColumnsNumber; j++) {
                flag = 1;
                expr = CalculateDistance(apLocation1, apLocation2, j0, i0, aDistanceMetric, flag);
                if (expr == 0)
                    apMatrixA[i + j * aRowsNumber] = sigma_square + aLocalTheta[3];
                else
                    apMatrixA[i + j * aRowsNumber] =
                            con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr); // Matern Function
                j0++;
            }
            i0++;
        }
    }
}