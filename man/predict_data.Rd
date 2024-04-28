
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file predict_data.Rd
% @brief roxygen2 documentation for the R Interface of R_ExaGeoStatPredictData function.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{predict_data}
\alias{predict_data}
\title{Predict Data function}

\usage{
predict_data(kernel, distance_matrix = "euclidean", estimated_theta,
dts, lts = 0, dimension = "2D", train_data, test_data)
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
 \item{distance_matrix}{ A string specifies the distance metric, either "euclidean" or "great_circle". Default is "euclidean"}
 \item{estimated_theta}{A list of estimated theta parameters}
 \item{dts}{A numeric value represents the time step size}
 \item{lts}{A numeric value represents the length step size. Default is 0}
 \item{dimension}{A string specifies the data dimension, either "2D" or "3D". Default is "2D"}
 \item{train_data}{ A numeric vector contains the locations and z measurements for training}
 \item{test_data}{ A numeric vector contains the locations for testing.}
}

\value{
 A vector of predicted z values
}

\description{
This function predicts data based on the provided kernel, distance matrix, estimated theta, and other parameters.
}

\examples{
ncores <- 2
ngpus <- 0
problem_size <- 4
dts <- 2
computation <- "exact"
hardware <- new(Hardware, computation, ncores, ngpus, 1, 1)
kernel <- "univariate_matern_stationary"
estimated_theta <- c(1,0.1,0.5)

z_value <- c( -1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
              -0.163880452049749520)
locations_x <- c(0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
              0.370473792629892440)
locations_y <- c(0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
              0.400778210116731537)
test_x <- c(0.347951, 0.62768)
test_y <- c(0.806332, 0.105196)

predict_data(train_data=list(locations_x, locations_y, z_value),
test_data=list(test_x, test_y), kernel=kernel, dts=dts,
estimated_theta=estimated_theta)
}
