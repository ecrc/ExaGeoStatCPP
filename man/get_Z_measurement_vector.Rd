// TODO: dots
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file get_Z_measurement_vector.Rd
% @brief roxygen2 documentation for the R Interface of get_Z_measurment_vector function.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{get_Z_measurement_vector}
\alias{get_Z_measurement_vector}
\title{Get descriptive Z values function}

\description{
Retrieves descriptive Z values from ExaGeoStat data based on type.
}

\usage{
get_Z_measurement_vector(data,type)
}

\arguments{
\item{data}{A list of ExaGeoStatData that contains the locations.}
\item{type}{A string specifies the type of descriptor value to retrieve (e.g., "Chameleon", "HiCMA").}
}

\value{
A numeric vector of descriptive Z values.
}

\examples{
ncores <- 2
ngpus <- 0
computation <- "exact"
p <- 1
q <- 1
hardware <- new(Hardware, computation, ncores, ngpus, p, q)

dimension = "3D"
problem_size <- 4
empty_data <- new(Data, problem_size, dimension)

dts <- 2
kernel <- "univariate_matern_stationary"
initial_theta <- c(1,0.1,0.5)
exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta,
problem_size=problem_size, dts=dts, dimension=dimension)
Z <- get_Z_measurement_vector(data=exageostat_data, type="chameleon")
}
