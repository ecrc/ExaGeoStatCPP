
% Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file model_data.Rd
% @brief roxygen2 documentation for the R Interface of R_ExaGeoStatModelData function.
% @version 1.1.0
% @author Mahmoud ElKarargy
% @date 2024-03-17

\name{model_data}
\alias{model_data}
\title{Model Data function}

\usage{
model_data(computation = "exact", kernel, distance_matrix = "euclidean", lb,
ub, tol = 4, mle_itr, dts, lts = 0, dimension = "2D", band = 0, max_rank = 500,
data = NULL, matrix = NULL, x = NULL, y = NULL, z = NULL)
}

\arguments{
\item{computation}{A string specifies the computation method, either "exact" or "dst" or "tlr". Default is "exact".}
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
\item{distance_matrix}{A string specifies the distance metric, either "euclidean" or "great_circle". Default is "euclidean".}
\item{lb}{A numeric value represents the lower bound for the computation.}
\item{ub}{A numeric value represents the upper bound for the computation.}
\item{tol}{A numeric value specifies the tolerance for the computation. Default is 4.}
\item{mle_itr}{A numeric value specifies the maximum number of iterations for the computation.}
\item{dts}{A numeric value represents the time step size.}
\item{lts}{A numeric value represents the length step size. Default is 0.}
\item{dimension}{A string specifies the data dimension, either "2D" or "3D". Default is "2D".}
\item{band}{A numeric value Bandwidth for band matrices, applicable in certain computational kernels, Default is 0.}
\item{max_rank}{A numeric value specifies the Maximum rank for low-rank approximations, Default is 500.}
\item{data}{A list of data vectors. Default is `R_NilValue`.}
\item{matrix}{A matrix object. Default is `R_NilValue`.}
\item{x}{A numeric vector. Default is `R_NilValue`.}
\item{y}{A numeric vector. Default is `R_NilValue`.}
\item{z}{A numeric vector. Default is `R_NilValue`.}
}

\value{
A vector contains the starting theta.
}

\description{
This function models data based on the provided computation method, kernel, distance matrix, and other parameters.
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
lower_bound <- c(0.1,0.1,0.1)
upper_bound <- c(5,5,5)

z_value <- c( -1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
              -0.163880452049749520)
locations_x <- c(0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
              0.370473792629892440)
locations_y <- c(0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
              0.400778210116731537)

theta <- model_data(kernel=kernel, lb=lower_bound, ub=upper_bound,
mle_itr=10, dts=dts, matrix=z_value, x=locations_x, y=locations_y)
}
