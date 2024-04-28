
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file simulate_data.Rd
% @brief roxygen2 documentation for the R Interface of R_ExaGeoStatLoadData function
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{simulate_data}
\alias{simulate_data}
\title{Simulate Data function}

\usage{
simulate_data(kernel, initial_theta, distance_matrix = "euclidean", problem_size,
seed = 0, dts, lts = 0, dimension = "2D", log_path = "", data_path = "",
observations_file = "", recovery_file = "")
}

\arguments{
  \item{kernel}{A string specifies the kernel to use. Available kernels include:
      \itemize{
        \item "BivariateMaternFlexible"
        \item "BivariateMaternParsimonious"
        \item "BivariateSpacetimeMaternStationary"
        \item "TrivariateMaternParsimonious"
        \item "UnivariateExpNonGaussian"
        \item "UnivariateMaternDbeta"
        \item "UnivariateMaternDdbetaBeta"
        \item "UnivariateMaternDdbetaNu"
        \item "UnivariateMaternDdnuNu"
        \item "UnivariateMaternDdsigmaSquare"
        \item "UnivariateMaternDdsigmaSquareBeta"
        \item "UnivariateMaternDdsigmaSquareNu"
        \item "UnivariateMaternDnu"
        \item "UnivariateMaternDsigmaSquare"
        \item "UnivariateMaternNonGaussian"
        \item "UnivariateMaternNuggetsStationary"
        \item "UnivariateMaternStationary"
        \item "UnivariatePowExpStationary"
        \item "UnivariateSpacetimeMaternStationary"
      }
    }
\item{initial_theta}{A list of initial theta parameters.}
\item{distance_matrix}{A string specifies the distance metric, either "euclidean" or "great_circle". Default is "euclidean".}
\item{problem_size}{A numeric value represents the size of the problem to simulate.}
\item{seed}{ A numeric value specifies the seed for random number generation. Default is 0.}
\item{dts}{ A numeric value represents the time step size.}
\item{lts}{ A numeric  value represents the length step size. Default is 0.}
\item{dimension}{A string specifies the data dimension, either "2D" or "3D". Default is "2D".}
\item{log_path}{A string specifies the path for logging.}
\item{data_path}{A string specifies the path for data storage.}
\item{observations_file}{A string specifies the file name for observations.}
\item{recovery_file}{A string specifies the file name for recovery.}
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
hardware <- new(Hardware, computation, ncores, ngpus, 1, 1)

dimension = "2D"
problem_size <- 4
empty_data <- new(Data, problem_size, dimension)

dts <- 2
kernel <- "univariate_matern_stationary"
initial_theta <- c(1,0.1,0.5)

exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta,
problem_size=problem_size, dts=dts, dimension=dimension)
}
