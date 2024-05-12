
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file fisher.Rd
% @brief roxygen2 documentation for the R Interface of R_ExaGeoStatFisher function.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{fisher}
\alias{fisher }
\title{Fisher function}

\usage{
fisher(kernel, distance_matrix = "euclidean", estimated_theta, dts,
lts = 0, dimension = "2D", train_data, test_data)
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
\item{distance_matrix}{ A string specifies the distance metric, either "euclidean" or "great_circle". Default is "euclidean".}
\item{estimated_theta}{A list of estimated theta parameters.}
\item{dts}{A numeric value represents the time step size.}
\item{lts}{A numeric value represents the length step size. Default is 0.}
\item{dimension}{A string specifies the data dimension, either "2D" or "3D". Default is "2D".}
\item{train_data}{ A numeric vector contains the locations and z measurements for training.}
\item{test_data}{ A numeric vector contains the locations for testing.}
}

\value{A vector contains the Fisher information matrix elements.}

\description{This function computes the Fisher information matrix for a given dataset and theta vector,
using a specified kernel and distance metric. It also allows for the inclusion of missing values and the specification of data dimensions.}

\examples{
dimension = "2D"
ncores <- 1
ngpus <- 0
dts <- 2
kernel <- "univariate_matern_stationary"
estimated_theta <- c(1,0.1,0.5)
computation <- "exact"
p <- 1
q <- 1
hardware <- new(Hardware, computation, ncores, ngpus, p, q)

z_value <- c(-1.272336140360187606, -2.590699695867695773,
             0.512142584178685967, -0.163880452049749520)
locations_x <- c(0.092042420080872822,  0.193041886015106440,
             0.330556191348134576, 0.181612878614480805)
locations_y <- c(0.928648813611047563, 0.103883421072709245,
             0.135790035858701447,  0.434683756771190977)

test_x <- c(0.347951, 0.62768)
test_y <- c(0.806332, 0.105196)

fisher_matrix <- fisher(train_data=list(locations_x, locations_y, z_value),
test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta)
}
