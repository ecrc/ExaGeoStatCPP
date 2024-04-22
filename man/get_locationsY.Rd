
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file get_locationsY.Rd
% @brief roxygen2 documentation for the R Interface of get_locationsY function for Y coordinates.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{get_locationsY}
\alias{get_locationsY}
\title{Get Locations Y Function}

\description{
This function retrieves the Y locations from the provided data.
}

\usage{
get_locationsY(data)
}

\arguments{
\item{data}{A list of ExaGeoStatData that contains the locations.}
}

\value{
A numeric vector of Y locations.
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
exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta,
problem_size=problem_size, dts=dts, dimension=dimension)
y <- get_locationsY(data=exageostat_data)
}
