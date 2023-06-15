
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdnuNu.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdnuNu.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternDdnuNu::UnivariateMaternDdnuNu() {
    /// TODO: FIX THEIR VALUES
    this->mP = 1;
    this->mParametersNumber = 3;
}

Kernel *UnivariateMaternDdnuNu::Create() {
    return new UnivariateMaternDdnuNu();
}

namespace exageostat::kernels {
    bool UnivariateMaternDdnuNu::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternDdnuNu", UnivariateMaternDdnuNu::Create);
}

void UnivariateMaternDdnuNu::GenerateCovarianceMatrix(double *apMatrixA, int &aRowsNumber, int &aColumnsNumber,
                                                        int &aRowOffset, int &aColumnOffset, Locations *apLocation1,
                                                        Locations *apLocation2, Locations *apLocation3,
                                                        double *aLocalTheta, int &aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double expr = 0.0;
    double con = 0.0;
    double nu_expr = 0.0;
    double nu_expr_dprime = 0.0;
    double sigma_square = aLocalTheta[0];
    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
                nu_expr = (1 - aLocalTheta[2]) * 1 / pow(2, aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2])
                          * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                          pow(2, 1 - aLocalTheta[2])
                          * (-1 / tgamma(aLocalTheta[2]) * gsl_sf_psi(aLocalTheta[2]) * pow(expr, aLocalTheta[2])
                             * gsl_sf_bessel_Knu(aLocalTheta[2], expr) + 1 / tgamma(aLocalTheta[2])
                                                                          * (pow(expr, aLocalTheta[2]) * log(expr)
                                                                             * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                                                                             pow(expr, aLocalTheta[2])
                                                                             * CalculateDerivativeBesselNu(aLocalTheta[2], expr)));

                nu_expr_dprime = (1 - aLocalTheta[2]) * 1 / pow(2, aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2])
                                 * pow(expr, aLocalTheta[2]) * CalculateDerivativeBesselNu(aLocalTheta[2], expr) +
                                 pow(2, 1 - aLocalTheta[2])
                                 * (-1 / tgamma(aLocalTheta[2]) * gsl_sf_psi(aLocalTheta[2]) * pow(expr, aLocalTheta[2])
                                    * CalculateDerivativeBesselNu(aLocalTheta[2], expr) + 1 / tgamma(aLocalTheta[2])
                                                                          * (pow(expr, aLocalTheta[2]) * log(expr)
                                                                             * CalculateDerivativeBesselNu(aLocalTheta[2], expr) +
                                                                             pow(expr, aLocalTheta[2])
                                                                             * CalculateSecondDerivativeBesselNu(aLocalTheta[2], expr)));
                apMatrixA[i + j * aRowsNumber] = (-0.5 * con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                                                  (1 - aLocalTheta[2]) / 2 * nu_expr
                                                  - (1 / aLocalTheta[2] + 0.5 * pow(aLocalTheta[2], 2)) * con * pow(expr, aLocalTheta[2]) *
                                                    gsl_sf_bessel_Knu(aLocalTheta[2], expr) - gsl_sf_psi(aLocalTheta[2])
                                                                                               * nu_expr + log(expr) * nu_expr +
                                                  nu_expr_dprime) * sigma_square;
            }
            j0++;
        }
        i0++;
    }
}