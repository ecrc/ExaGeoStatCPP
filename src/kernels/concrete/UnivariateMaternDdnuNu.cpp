
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
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
#include<cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_psi.h>


using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;


void UnivariateMaternDdnuNu::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                        int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                        Locations *apLocation2, Locations *apLocation3,
                                                        double *apLocalTheta, int aDistanceMetric) {

}