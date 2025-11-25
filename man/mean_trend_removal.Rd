
% Copyright (c) 2017-2025 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file mean_trend_removal.Rd
% @brief roxygen2 documentation for the R Interface of R_ExaGeoStatMeanTrendRemoval function
% @version 2.0.0
% @author Mahmoud ElKarargy
% @date 2025-11-25

\name{mean_trend_removal}
\alias{mean_trend_removal}
\title{Mean Trend Removal for Climate Data}

\usage{
mean_trend_removal(lon, startyear, endyear, dts, data_path, forcing_data_path,
results_path, starting_theta, lb, ub, tolerance = 7, max_mle_iterations = 30,
cores = 1, gpus = 0, p = 1, q = 1, log_path = "")
}

\arguments{
  \item{lon}{A numeric value specifying the number of longitude points in the grid.}
  \item{startyear}{A numeric value specifying the start year for the analysis period.}
  \item{endyear}{A numeric value specifying the end year for the analysis period.}
  \item{dts}{A numeric value specifying the tile size for dense matrix computations.}
  \item{data_path}{A string specifying the path to the ERA5/climate data directory containing NetCDF files.}
  \item{forcing_data_path}{A string specifying the path to the forcing data CSV file.}
  \item{results_path}{A string specifying the path where results will be saved (z_*.csv and params.csv).}
  \item{starting_theta}{A numeric vector of initial parameter values for optimization.}
  \item{lb}{A numeric vector of lower bounds for optimization parameters.}
  \item{ub}{A numeric vector of upper bounds for optimization parameters.}
  \item{tolerance}{A numeric value specifying tolerance level for the optimization algorithm (10^(-tolerance)). Default is 7.}
  \item{max_mle_iterations}{A numeric value specifying maximum number of MLE iterations. Default is 30.}
  \item{cores}{A numeric value specifying the number of CPU cores to use. Default is 1.}
  \item{gpus}{A numeric value specifying the number of GPUs to use. Default is 0.}
  \item{p}{A numeric value specifying the P dimension of the MPI process grid. Default is 1.}
  \item{q}{A numeric value specifying the Q dimension of the MPI process grid. Default is 1.}
  \item{log_path}{A string specifying the path for logging. Optional.}
}

\value{
A list containing:
  \item{results_path}{Path where the output files (z_*.csv and params.csv) are saved.}
  \item{num_time_points}{Number of time points processed (number of z_*.csv files created).}
  \item{num_locations}{Number of spatial locations processed.}
}

\description{
This function performs Mean Trend Removal on climate time series data, removing mean trends and
producing normalized residuals for further analysis. The function works with both StarPU and PaRSEC
runtimes and is typically the first step in the Global Climate Emulator pipeline.

The function reads climate data from NetCDF files, applies mean trend removal using the TrendModel
kernel, and outputs normalized residuals (z_*.csv files) and optimized parameters (params.csv) that
can be used as input for the Climate Emulator.
}

\details{
Mean Trend Removal is a data preprocessing step that:
\itemize{
  \item Reads hourly temperature data from ERA5 NetCDF files
  \item Reads forcing data (e.g., solar radiation patterns)
  \item Builds a design matrix using harmonic functions
  \item Optimizes parameters using Maximum Likelihood Estimation (MLE)
  \item Computes normalized residuals by removing the estimated mean trend
  \item Outputs results for each time point and location
}

The function requires MPI for parallel processing. Use \code{mpirun} or equivalent to launch:
\code{mpirun -n 2 Rscript your_script.R}

The function only supports the "TrendModel" kernel, which is specifically designed for climate
mean trend analysis.
}

\section{Output Files}{
The function creates the following outputs in \code{results_path}:
\itemize{
  \item \code{z_YYYY.csv}: Normalized residuals for each time point (one file per time point)
  \item \code{params.csv}: Optimized parameters for each location
}
}

\section{Runtime Requirements}{
\itemize{
  \item Requires MPI (2 or more processes recommended)
  \item Works with both StarPU and PaRSEC runtimes
  \item NetCDF library required for reading ERA5 data
  \item NLOPT library required for optimization
}
}

\examples{
\dontrun{
# Example: Process 3 years of data for 144 longitude points
# Note: Requires MPI - run with: mpirun -n 2 Rscript script.R

library(ExaGeoStatCPP)

result <- mean_trend_removal(
  lon = 144,
  startyear = 2020,
  endyear = 2022,
  dts = 720,
  data_path = "/path/to/ERA_data/",
  forcing_data_path = "/path/to/forcing_new.csv",
  results_path = "/path/to/results/",
  starting_theta = c(0.9),
  lb = c(0.001),
  ub = c(0.95),
  tolerance = 7,
  max_mle_iterations = 30,
  cores = 4,
  gpus = 0,
  p = 1,
  q = 2
)

cat("Results saved to:", result$results_path, "\n")
cat("Time points processed:", result$num_time_points, "\n")
cat("Locations processed:", result$num_locations, "\n")
}
}

\seealso{
\code{\link{climate_emulator}} for the next step in the climate modeling pipeline.
}

