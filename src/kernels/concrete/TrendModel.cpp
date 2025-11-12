
// Copyright (c) 2017-2025 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TrendModel.cpp
 * @brief Implementation of the TrendModel kernel.
 * @version 2.0.0
 * @author Ali Hakam
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-11-12
**/

#include <kernels/concrete/TrendModel.hpp>


using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
TrendModel<T>::TrendModel() {
    this->mP = 1;
    this->mParametersNumber = 1;
}

template<typename T>
Kernel<T> *TrendModel<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["TrendModel"] = 1;
    return new TrendModel();
}

namespace exageostat::kernels {
    template<typename T> bool TrendModel<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "TrendModel", TrendModel<T>::Create);
}

template<typename T>
void
TrendModel<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                        const int &aRowOffset, const int &aColumnOffset, Locations<T> &aLocation1,
                                        Locations<T> &aLocation2, Locations<T> &aLocation3, T *aLocalTheta,
                                        const int &aDistanceMetric) {
    int i, j;
    int row_offset = aRowOffset;
    int column_offset = aColumnOffset;
    T theta_1 = aLocalTheta[1];
    T theta_2 = aLocalTheta[2];
    T *forcing_theta = &aLocalTheta[3];

    double i_x = row_offset + 1;
    int ty;
    double theta_pow = 1.0;
    double sum = 0.0;

    for (i = 0; i < aRowsNumber; i++) {
        column_offset = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            if (column_offset == 0) {
                apMatrixA[i + j * aRowsNumber] = 1.0;
            }
            else if( column_offset==1 ) {
                apMatrixA[i + j * aRowsNumber] = forcing_theta[(int) (row_offset / theta_1) +
                                                               238]; // 190 for 1940 -- 238 for 1988
            } else if(column_offset==2) {
                ty = (row_offset) / theta_1;
                for (int k = 0; k < ty + 238; k++) { // 190 for 1940 -- 238 for 1988
                    for (int kk = k; kk < ty + 237; kk++) { // 189 for 1940 -- 237 for 1988
                        theta_pow *= aLocalTheta[0];
                    }
                    sum += theta_pow * forcing_theta[k];
                    theta_pow = 1;
                }
                apMatrixA[i + j * aRowsNumber] = (1 - aLocalTheta[0]) * sum;
                sum = 0;
                theta_pow = 1;
            } else {
			    if( j%2==0 ) {
				    apMatrixA[i + j * aRowsNumber]=sin(2.0 * PI * (i_x) * (floor((column_offset-3.0)/2.0)+1.0) / (theta_1));
			    } else{
				    apMatrixA[i + j * aRowsNumber]=cos(2.0 * PI * (i_x) * ((column_offset-3.0)/2.0+1.0) / (theta_1));
                }
		    }
		    column_offset++;
        }
        row_offset++;
        i_x++;
    }
}