
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateMaternFlexible.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/BivariateMaternFlexible.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

template<typename T>
BivariateMaternFlexible<T>::BivariateMaternFlexible() {
    this->mP = 2;
    this->mParametersNumber = 11;
}

template<typename T>
Kernel<T> *BivariateMaternFlexible<T>::Create() {
    return new BivariateMaternFlexible();
}

namespace exageostat::kernels {
    template<typename T> bool BivariateMaternFlexible<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "BivariateMaternFlexible", BivariateMaternFlexible<T>::Create);
}

template<typename T>
void BivariateMaternFlexible<T>::GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber,
                                                       int &aRowOffset, int &aColumnOffset, Locations<T> *apLocation1,
                                                       Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                       T *aLocalTheta, int &aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0;
    double expr1 = 0.0, expr2 = 0.0, expr12 = 0.0;
    double con1 = 0.0, con2 = 0.0, con12 = 0.0, scale12 = 0.0, rho = 0.0, nu12 = 0.0, sigma_square11 = 0.0, sigma_square22 = 0.0;
    double scale1 = aLocalTheta[0], scale2 = aLocalTheta[1], nu1 = aLocalTheta[4], nu2 = aLocalTheta[5];

    scale12 = pow(0.5 * (pow(scale1, -2) + pow(scale2, -2)) + aLocalTheta[2] * (1 - aLocalTheta[3]),
                  -0.5); //Remark 1 (c) of Apanasovich et al. (2012)

    nu12 = 0.5 * (nu1 + nu2) + aLocalTheta[6] * (1 - aLocalTheta[7]); //Theorem 1 (i) of Apanasovich et al. (2012).

    rho = aLocalTheta[8] * aLocalTheta[9] * aLocalTheta[10] *
          pow(scale12, 2 * aLocalTheta[6] + (nu1 + nu2))
          * tgamma(0.5 * (nu1 + nu2) + 1) * tgamma(nu12) /
          tgamma(nu12 + 1); //Equation (8) of Apanasovich et al. (2012).

    sigma_square11 = aLocalTheta[8] * aLocalTheta[8] *
                     pow(scale1, 2 * aLocalTheta[6] + nu1 + nu1) *
                     tgamma(nu1); //Equation (8) of Apanasovich et al. (2012).

    sigma_square22 = aLocalTheta[9] * aLocalTheta[9] *
                     pow(scale2, 2 * aLocalTheta[6] + nu2 + nu2) *
                     tgamma(nu2); //Equation (8) of Apanasovich et al. (2012).

    con1 = pow(2, (nu1 - 1)) * tgamma(nu1);
    con1 = 1.0 / con1;
    con1 = sigma_square11 * con1;

    con2 = pow(2, (nu2 - 1)) * tgamma(nu2);
    con2 = 1.0 / con2;
    con2 = sigma_square22 * con2;

    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho * con12;

    i0 /= 2;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i += 2) {
        j0 = aColumnOffset / 2;
        for (j = 0; j < aColumnsNumber; j += 2) {
            expr1 = this->CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, flag) / scale1;
            expr2 = this->CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, flag) / scale2;
            expr12 = this->CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, flag) / scale12;
            if (expr1 == 0) {
                apMatrixA[i + j * aRowsNumber] = aLocalTheta[0];
                if (((i + 1) + j * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = rho;
                }
                if ((i + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[i + (j + 1) * aRowsNumber] = rho;
                }
                if (((i + 1) + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = aLocalTheta[1];
                }
            } else {
                apMatrixA[i + j * aRowsNumber] = con1 * pow(expr1, nu1)
                                                 * gsl_sf_bessel_Knu(nu1, expr1);

                if (((i + 1) + j * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = con12 * pow(expr12, nu12) * gsl_sf_bessel_Knu(nu12, expr12);
                }
                if ((i + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[i + (j + 1) * aRowsNumber] = con12 * pow(expr12, nu12) * gsl_sf_bessel_Knu(nu12, expr12);
                }
                if (((i + 1) + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = con2 * pow(expr2, nu2) * gsl_sf_bessel_Knu(nu2, expr2);
                }
            }
            j0++;
        }
        i0++;
    }

}
