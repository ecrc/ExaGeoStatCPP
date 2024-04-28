
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file Hardware.Rd
% @brief roxygen2 documentation for the R Interface of ExaGeoStatHardware class.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{ExaGeoStatHardware}
\alias{ExaGeoStatHardware}
\title{ExaGeoStatHardware Class}

\description{
The ExaGeoStatHardware class represents a hardware component in the ExaGeoStat system.
It is initialized with computation mode, and two integers represents number of CPU cores and number of GPU cores.
}

\section{Constructor}{
 \code{\link{ExaGeoStatHardware}} Creates a new instance of the
 \code{ExaGeoStatHardware} class by calling:
 \preformatted{
    new(Hardware, computation, ncores, ngpus, p, q)
 }
 \describe{
     \item{\code{computation}}{A string specifies the computation method, either "exact" or "dst" or "tlr".}
     \item{\code{ncores}}{An integer represents number of CPU cores.}
     \item{\code{ngpus}}{An integer represents number of GPU cores.}
     \item{\code{p}}{An integer represents P grid dimension.}
     \item{\code{q}}{An integer represents Q grid dimension.}
}
}

\section{Methods}{
\subsection{finalize_hardware}{
\code{finalize_hardware()} manually finalizes the hardware by resetting the context.
}
}

\value{
An object of class \code{ExaGeoStatHardware} represents a hardware component with the specified component and number of CPU cores and GPU cores.
}

\examples{
ncores <- 2
ngpus <- 0
p <- 2
q <- 2
computation <- "exact"
hardware <- new(Hardware, computation, ncores, ngpus, p, q)
hardware$finalize_hardware()
}

\keyword{S4 class}
