
% Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file Hardware.Rd
% @brief roxygen2 documentation for the R wrapper of ExaGeoStatHardware class.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{ExaGeoStatHardware}
\alias{ExaGeoStatHardware}
\title{ExaGeoStatHardware Class}

\description{
The ExaGeoStatHardware class represents a hardware component in the ExaGeoStat system.
It is initialized with computation mode, and two integers representing number of CPU cores and number of GPU cores.
}

\section{Constructor}{
 \code{\link{ExaGeoStatHardware}} Creates a new instance of the \code{ExaGeoStatHardware} class.
 \code{ExaGeoStatHardware(computation, num_of_cpus, num_of_gpus)}
 \describe{
 \item{\code{computation}}{A string representing the computation mode.}
 \item{\code{num_of_cpus}}{An integer representing number of CPU cores.}
 \item{\code{num_of_gpus}}{An integer representing number of GPU cores.}
}
}

\section{Methods}{
\subsection{finalize_hardware}{
\code{finalize_hardware()} Manually finalizes the hardware by resetting the context.
}
}

\value{
An object of class \code{ExaGeoStatHardware} representing a hardware component with the specified component and number of CPU cores and GPU cores.
}

\examples{
ncores <- 2
ngpus <- 0
computation <- "exact"
hardware <- new(Hardware, computation, ncores, ngpus)

hardware$finalize_hardware()
}

\keyword{S4 class}
