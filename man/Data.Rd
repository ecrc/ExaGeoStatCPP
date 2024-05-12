
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file Data.Rd
% @brief roxygen2 documentation for the R Interface of ExaGeoStatData class
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{ExaGeoStatData}
\alias{ExaGeoStatData}
\title{ExaGeoStatData Class}

\description{
The \code{ExaGeoStatData} class is designed to facilitate the handling and manipulation of geospatial statistics data within the ExaGeoStat framework.
It provides a structured way to store and manage data points based on their dimensions.
Instances of this class can hold a specified number of location points, and support different data dimensions including two-dimensional (2D), three-dimensional (3D), and spatiotemporal (ST) configurations.
}

\section{Constructor}{
  \code{\link{ExaGeoStatData}} Creates a new instance of the
  \code{ExaGeoStatData} class by calling:
  \preformatted{
    new(ExaGeoStatData, problem_size, dimension)
  }
  \describe{
    \item{\code{size}}{An integer represents the size of the locations data.}
    \item{\code{dimension}}{A string represents the dimensions of the data.
    Available dimensions are "2D", "3D", and "ST".}
  }
}

\value{
An object of class \code{ExaGeoStatData} represents a data component with the specified size and dimension.
}

\examples{
problem_size <- 4
dimension = "3D"
empty_data <- new(Data, problem_size, dimension)
}

\keyword{S4 class}
