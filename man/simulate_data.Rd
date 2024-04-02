
% Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file simulate_data.Rd
% @brief roxygen2 documentation for the R wrapper of R_ExaGeoStatLoadData function
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{simulate_data}
\alias{simulate_data}
\title{This function simulates data based on the provided computation method, kernel, distance matrix, and other parameters.}

\usage{
simulate_data(kernel, initial_theta, distance_matrix = "euclidean", problem_size, seed = 0, dts, lts = 0, dimension = "2D", log_path = "", data_path = "", observations_file = "", recovery_file = "")
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
\item{initial_theta}{A list of initial theta parameters.}
\item{distance_matrix}{A string specifying the distance metric, either "euclidean" or "great_circle". Default is "euclidean".}
\item{problem_size}{A numeric value representing the size of the problem to simulate.}
\item{seed}{ A numeric value specifying the seed for random number generation. Default is 0.}
\item{dts}{ A numeric value representing the time step size.}
\item{lts}{ A numeric  value representing the length step size. Default is 0.}
\item{dimension}{A string specifying the data dimension, either "2D" or "3D". Default is "2D".}
\item{log_path}{A string specifying the path for logging. Default is "".}
\item{data_path}{A string specifying the path for data storage. Default is "".}
\item{observations_file}{A string specifying the file name for observations. Default is "".}
\item{recovery_file}{A string specifying the file name for recovery. Default is "".}
}

\value{
A pointer to ExaGeoStatData object that contains the loaded data.
}

\description{
This function loads data into an ExaGeoStatData object using the provided configuration and computational settings.
}

\examples{
ncores <- 2
ngpus <- 0
computation <- "exact"
hardware <- new(Hardware, computation, ncores, ngpus)

dimension = "2D"
problem_size <- 4
empty_data <- new(Data, problem_size, dimension)

dts <- 2
kernel <- "univariate_matern_stationary"
initial_theta <- c(1,0.1,0.5)

exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta, problem_size=problem_size, dts=dts, dimension=dimension)
}
