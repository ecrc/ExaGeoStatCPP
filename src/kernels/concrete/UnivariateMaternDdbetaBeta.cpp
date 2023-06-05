
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdbetaBeta.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdbetaBeta.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternDdbetaBeta::UnivariateMaternDdbetaBeta() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

Kernel *UnivariateMaternDdbetaBeta::Create() {
    return new UnivariateMaternDdbetaBeta();
}

namespace exageostat::kernels {
    bool UnivariateMaternDdbetaBeta::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternDdbetaBeta", UnivariateMaternDdbetaBeta::Create);
}

void UnivariateMaternDdbetaBeta::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                          int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                          Locations *apLocation2, Locations *apLocation3,
                                                          double *aLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0;
    double expr = 0.0;
    double con = 0.0;
    double beta_expr = 0.0;
    double beta_expr_prime = 0.0;
    double sigma_square = aLocalTheta[0];
    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
                beta_expr = -aLocalTheta[2] / aLocalTheta[1] * pow(expr, aLocalTheta[2])
                            * gsl_sf_bessel_Knu(aLocalTheta[2], expr) - pow(expr, aLocalTheta[2])
                                                                         * (aLocalTheta[2] / expr *
                                                                            gsl_sf_bessel_Knu(aLocalTheta[2], expr) -
                                                                            gsl_sf_bessel_Knu(aLocalTheta[2] + 1,
                                                                                              expr)) * expr /
                                                                         aLocalTheta[1];

                beta_expr_prime = -aLocalTheta[2] / aLocalTheta[1] * pow(expr, aLocalTheta[2])
                                  * (aLocalTheta[2] / expr * gsl_sf_bessel_Knu(aLocalTheta[2], expr) -
                                     gsl_sf_bessel_Knu(aLocalTheta[2] + 1, expr)) - pow(expr, aLocalTheta[2])
                                                                                     * (-0.5 *
                                                                                        ((aLocalTheta[2] / expr *
                                                                                          gsl_sf_bessel_Knu(
                                                                                                  aLocalTheta[2],
                                                                                                  expr) -
                                                                                          gsl_sf_bessel_Knu(
                                                                                                  aLocalTheta[2] + 1,
                                                                                                  expr)) -
                                                                                         pow(expr, aLocalTheta[2])
                                                                                         + (aLocalTheta[2] / expr *
                                                                                            gsl_sf_bessel_Knu(
                                                                                                    aLocalTheta[2],
                                                                                                    expr) -
                                                                                            gsl_sf_bessel_Knu(
                                                                                                    aLocalTheta[2] + 1,
                                                                                                    expr)) -
                                                                                         pow(expr, aLocalTheta[2])))
                                                                                     * expr /
                                                                                     aLocalTheta[1];
                apMatrixA[i + j * aRowsNumber] =
                        (aLocalTheta[2] / pow(aLocalTheta[1], 2) * pow(expr, aLocalTheta[2]) *
                         gsl_sf_bessel_Knu(aLocalTheta[2], expr)
                         - aLocalTheta[2] / aLocalTheta[1] * beta_expr +
                         2 * expr / pow(aLocalTheta[1], 2) * pow(expr, aLocalTheta[2])
                         * (aLocalTheta[2] / expr * gsl_sf_bessel_Knu(aLocalTheta[2], expr) -
                            gsl_sf_bessel_Knu(aLocalTheta[2] + 1, expr)) - expr / aLocalTheta[1] * beta_expr_prime) *
                        sigma_square * con; // derivative with respect to beta
            }
            j0++;
        }
        i0++;
    }
}