
% Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file Data.Rd
% @brief roxygen2 documentation for the R wrapper of ExaGeoStatData class
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{ExaGeoStatData}
\alias{ExaGeoStatData}
\title{ExaGeoStatData Class}

\description{
The ExaGeoStatData class represents a data component in the ExaGeoStat system, that manages geo-statistical data with functions
for location and descriptor manipulation. It is initialized with the size and dimension of the data.
}

\section{Constructor}{
 \code{\link{ExaGeoStatData}} Creates a new instance of the \code{ExaGeoStatData} class.
 \code{ExaGeoStatData(size, dimension)}
 \describe{
    \item{\code{size}}{An integer representing the size of the locations data.}
    \item{\code{dimension}}{A string representing the dimensions of the data. - available dimension ("2D", "3D", "ST")}
    }
}

\value{
An object of class \code{ExaGeoStatData} representing a data component with the specified size and dimension.
}

\examples{
problem_size <- 4
dimension = "3D"
empty_data <- new(Data, problem_size, dimension)
}

\keyword{S4 class}
