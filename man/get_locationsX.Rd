
% Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file get_locationsX.Rd
% @brief roxygen2 documentation for the R wrapper of get_locationsX function for X coordinates.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{get_locationsX}
\alias{get_locationsX}
\title{Get X Locations}

\description{
Retrieves X coordinates of locations from ExaGeoStatData object.
}

\usage{
get_locationsX(data)
}

\arguments{
\item{data}{A list of ExaGeoStatData that contains the locations.}
}

\value{
A numeric vector of X locations.
}

\examples{
ncores <- 1
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

x <- get_locationsX(data=exageostat_data)
}
