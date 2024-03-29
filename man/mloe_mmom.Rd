
% Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file mloe_mmom.Rd
% @brief roxygen2 documentation for the R wrapper of R_ExaGeoStatMLOE_MMOM function
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{mloe_mmom}
\alias{mloe_mmom}
\title{Mean Misspecification of the Mean Square Error (MMOM) and Mean Loss of Efficiency (MLOE) using exact method.}

\usage{
mloe_mmom(data=empty_data, matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, estimated_theta=estimated_theta, true_theta=estimated_theta, z_miss=3 )
}

\arguments{
 \item{kernel}{A string specifying the kernel to use - available kernels
  (        "BivariateMaternFlexible",
           "BivariateMaternParsimonious",
           "BivariateSpacetimeMaternStationary",
           "TrivariateMaternParsimonious",
           "UnivariateExpNonGaussian",
           "UnivariateMaternDbeta",
           "UnivariateMaternDdbetaBeta",
           "UnivariateMaternDdbetaNu",
           "UnivariateMaternDdnuNu",
           "UnivariateMaternDdsigmaSquare",
           "UnivariateMaternDdsigmaSquareBeta",
           "UnivariateMaternDdsigmaSquareNu",
           "UnivariateMaternDnu",
           "UnivariateMaternDsigmaSquare",
           "UnivariateMaternNonGaussian",
           "UnivariateMaternNuggetsStationary",
           "UnivariateMaternStationary",
           "UnivariatePowExpStationary",
           "UnivariateSpacetimeMaternStationary",
           "bivariate_matern_flexible",
           "bivariate_matern_parsimonious",
           "bivariate_spacetime_matern_stationary",
           "trivariate_matern_parsimonious",
           "univariate_exp_non_gaussian",
           "univariate_matern_dbeta",
           "univariate_matern_ddbeta_beta",
           "univariate_matern_ddbeta_nu",
           "univariate_matern_ddnu_nu",
           "univariate_matern_ddsigma_square",
           "univariate_matern_ddsigma_square_beta",
           "univariate_matern_ddsigma_square_nu",
           "univariate_matern_dnu",
           "univariate_matern_dsigma_square",
           "univariate_matern_non_gaussian",
           "univariate_matern_nuggets_stationary",
           "univariate_matern_stationary",
           "univariate_pow_exp_stationary",
           "univariate_spacetime_matern_stationary"
  )}
\item{distance_matrix}{ A string specifying the distance metric, either "euclidean" or "great_circle". Default is "euclidean"}
\item{estimated_theta}{A list of estimated theta parameters}
\item{dts}{A numeric value representing the time step size}
\item{lts}{A numeric value representing the length step size. Default is 0}
\item{dimension}{A string specifying the data dimension, either "2D" or "3D". Default is "2D"}
\item{z_miss}{ A numeric value representing the missing value indicator}
\item{data}{ A list of data vectors}
\item{matrix}{ A matrix object. Default is `R_NilValue`}
\item{x}{A numeric vector. Default is `R_NilValue`}
\item{y}{A numeric vector. Default is `R_NilValue`}
\item{z}{A numeric vector. Default is `R_NilValue`}
}

\value{
A vector of MLOE/MMOM values
}

\description{
This function calculates Mean Misspecification of the Mean Square Error (MMOM) and Mean Loss of Efficiency (MLOE).
}

\examples{
ncores <- 2
ngpus <- 0
computation <- "exact"

hardware <- new(Hardware, computation, ncores, ngpus)

problem_size <- 4
dimension = "2D"
dts <- 2
z_value <- c( -1.272336140360187606, -2.590699695867695773, 0.512142584178685967, -0.163880452049749520)
x_locations <- c(0.193041886015106440, 0.330556191348134576, 0.181612878614480805, 0.370473792629892440)
y_locations <- c(0.103883421072709245, 0.135790035858701447, 0.434683756771190977, 0.400778210116731537)

empty_data <- new(Data, problem_size, dimension)
mloe_mmom(data=empty_data, matrix=z_value, x=x_locations, y=y_locations, kernel="univariate_matern_stationary", dts = 2, estimated_theta=c(1,0.1,0.5), true_theta=c(1,0.1,0.5), z_miss=3 )
}
