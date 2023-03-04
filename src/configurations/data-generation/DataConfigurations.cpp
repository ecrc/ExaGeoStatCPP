
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataConfigurations.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-03
**/
#include <iostream>
#include <configurations/data-generation/DataConfigurations.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace std;

string DataConfigurations::GetKernel() {
    return this->mKernel;
}

void DataConfigurations::SetKernel(std::string aKernel) {
    this->mKernel = aKernel;
}

void DataConfigurations::SetIsSynthetic(bool aIsSynthetic) {
    this->mIsSynthetic = aIsSynthetic;
}

bool DataConfigurations::GetIsSynthetic() {
    return this->mIsSynthetic;
}

void DataConfigurations::CheckKernelValue(std::string aKernel) {
    if (aKernel != "univariate_matern_stationary"
        and aKernel != "univariate_matern_non_stationary"
        and aKernel != "bivariate_matern_flexible"
        and aKernel != "bivariate_matern_parsimonious"
        and aKernel != "bivariate_matern_parsimonious_profile"
        and aKernel != "univariate_matern_nuggets_stationary"
        and aKernel != "univariate_spacetime_matern_stationary"
        and aKernel != "univariate_matern_dsigma_square"
        and aKernel != "univariate_matern_dnu"
        and aKernel != "univariate_matern_dbeta"
        and aKernel != "univariate_matern_ddsigma_square"
        and aKernel != "univariate_matern_ddsigma_square_beta"
        and aKernel != "univariate_matern_ddsigma_square_nu"
        and aKernel != "univariate_matern_ddbeta_beta"
        and aKernel != "univariate_matern_ddbeta_nu"
        and aKernel != "univariate_matern_ddnu_nu"
        and aKernel != "bivariate_spacetime_matern_stationary"
        and aKernel != "univariate_matern_non_gaussian"
        and aKernel != "univariate_exp_non_gaussian"
        and aKernel != "bivariate_spacetime_matern_stationary"
        and aKernel != "trivariate_matern_parsimonious"
        and aKernel != "trivariate_matern_parsimonious_profile"
        and aKernel != "univariate_matern_non_stat") {

        throw range_error("Invalid value for Kernel. Please check manual.");
    }
}
