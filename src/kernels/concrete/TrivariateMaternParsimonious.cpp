
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TrivariateMaternParsimonious.cpp
 * @brief Implementation of the BivariateMaternParsimonious kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#include <kernels/concrete/TrivariateMaternParsimonious.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
TrivariateMaternParsimonious<T>::TrivariateMaternParsimonious() {
    this->mP = 3;
    this->mParametersNumber = 10;
}

template<typename T>
Kernel<T> *TrivariateMaternParsimonious<T>::Create() {
    return new TrivariateMaternParsimonious();
}

namespace exageostat::kernels {
    template<typename T> bool TrivariateMaternParsimonious<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "TrivariateMaternParsimonious", TrivariateMaternParsimonious<T>::Create);
}

template<typename T>
void TrivariateMaternParsimonious<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                               const int &aColumnsNumber,
                                                               const int &aRowOffset, const int &aColumnOffset,
                                                               dataunits::Locations<T> &aLocation1,
                                                               dataunits::Locations<T> &aLocation2,
                                                               dataunits::Locations<T> &aLocation3, T *aLocalTheta,
                                                               const int &aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double con1, con2, con3, con12, con13, con23, rho12, rho13, rho23, nu12, nu13, nu23;

    con1 = pow(2, (aLocalTheta[4] - 1)) * tgamma(aLocalTheta[4]);
    con1 = 1.0 / con1;
    con1 = aLocalTheta[0] * con1;

    con2 = pow(2, (aLocalTheta[5] - 1)) * tgamma(aLocalTheta[5]);
    con2 = 1.0 / con2;
    con2 = aLocalTheta[1] * con2;

    con3 = pow(2, (aLocalTheta[6] - 1)) * tgamma(aLocalTheta[6]);
    con3 = 1.0 / con3;
    con3 = aLocalTheta[2] * con3;

    // The average
    nu12 = 0.5 * (aLocalTheta[4] + aLocalTheta[5]);
    nu13 = 0.5 * (aLocalTheta[4] + aLocalTheta[6]);
    nu23 = 0.5 * (aLocalTheta[5] + aLocalTheta[6]);

    rho12 = aLocalTheta[7] * sqrt((tgamma(aLocalTheta[4] + 1) * tgamma(aLocalTheta[5] + 1)) /
                                  (tgamma(aLocalTheta[4]) * tgamma(aLocalTheta[5]))) *
            tgamma(nu12) / tgamma(nu12 + 1);

    rho13 = aLocalTheta[8] * sqrt((tgamma(aLocalTheta[4] + 1) * tgamma(aLocalTheta[6] + 1)) /
                                  (tgamma(aLocalTheta[4]) * tgamma(aLocalTheta[6]))) *
            tgamma(nu13) / tgamma(nu13 + 1);

    rho23 = aLocalTheta[9] * sqrt((tgamma(aLocalTheta[5] + 1) * tgamma(aLocalTheta[6] + 1)) /
                                  (tgamma(aLocalTheta[5]) * tgamma(aLocalTheta[6]))) *
            tgamma(nu23) / tgamma(nu23 + 1);

    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho12 * sqrt(aLocalTheta[0] * aLocalTheta[1]) * con12;

    con13 = pow(2, (nu13 - 1)) * tgamma(nu13);
    con13 = 1.0 / con13;
    con13 = rho13 * sqrt(aLocalTheta[0] * aLocalTheta[2]) * con13;

    con23 = pow(2, (nu23 - 1)) * tgamma(nu23);
    con23 = 1.0 / con23;
    con23 = rho23 * sqrt(aLocalTheta[1] * aLocalTheta[2]) * con23;

    i0 /= 3;
    int matrix_size = aRowsNumber * aColumnsNumber;
    int index;
    int flag = 0;

    for (i = 0; i < aRowsNumber - 1; i += 3) {
        j0 = aColumnOffset / 3;
        for (j = 0; j < aColumnsNumber - 1; j += 3) {
            expr = this->CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[3];

            if (expr == 0) {

                apMatrixA[i + j * aRowsNumber] = aLocalTheta[0];
                index = (i + 1) + j * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = rho12 * sqrt(aLocalTheta[0] * aLocalTheta[1]);
                }
                index = i + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[i + (j + 1) * aRowsNumber] = rho12 * sqrt(aLocalTheta[0] * aLocalTheta[1]);
                }
                index = (i + 2) + j * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 2) + j * aRowsNumber] = rho13 * sqrt(aLocalTheta[0] * aLocalTheta[2]);
                }
                index = i + (j + 2) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[i + (j + 2) * aRowsNumber] = rho13 * sqrt(aLocalTheta[0] * aLocalTheta[2]);
                }

                index = (i + 1) + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = aLocalTheta[1];
                }

                index = (i + 1) + (j + 2) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + (j + 2) * aRowsNumber] = rho23 * sqrt(aLocalTheta[1] * aLocalTheta[2]);
                }
                index = (i + 2) + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 2) + (j + 1) * aRowsNumber] = rho23 * sqrt(aLocalTheta[1] * aLocalTheta[2]);
                }
                index = (i + 2) + (j + 2) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 2) + (j + 2) * aRowsNumber] = aLocalTheta[2];
                }
            } else {
                apMatrixA[i + j * aRowsNumber] =
                        con1 * pow(expr, aLocalTheta[4]) * gsl_sf_bessel_Knu(aLocalTheta[4], expr);

                index = (i + 1) + j * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12, expr);
                }
                index = i + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[i + (j + 1) * aRowsNumber] =
                            con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12, expr);
                }
                index = (i + 2) + j * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 2) + j * aRowsNumber] = con13 * pow(expr, nu13) * gsl_sf_bessel_Knu(nu13, expr);
                }
                index = i + (j + 2) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[i + (j + 2) * aRowsNumber] = con13 * pow(expr, nu13) * gsl_sf_bessel_Knu(nu13, expr);
                }
                index = (i + 1) + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] =
                            con2 * pow(expr, aLocalTheta[5]) * gsl_sf_bessel_Knu(aLocalTheta[5], expr);
                }

                index = (i + 1) + (j + 2) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + (j + 2) * aRowsNumber] =
                            con23 * pow(expr, nu23) * gsl_sf_bessel_Knu(nu23, expr);
                }
                index = (i + 2) + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 2) + (j + 1) * aRowsNumber] =
                            con23 * pow(expr, nu23) * gsl_sf_bessel_Knu(nu23, expr);
                }
                index = (i + 2) + (j + 2) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 2) + (j + 2) * aRowsNumber] =
                            con3 * pow(expr, aLocalTheta[6]) * gsl_sf_bessel_Knu(aLocalTheta[6], expr);
                }
            }
            j0++;
        }
        i0++;
    }
}
