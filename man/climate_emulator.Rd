
% Copyright (c) 2017-2025 King Abdullah University of Science and Technology,
% All rights reserved.
% ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

% @file climate_emulator.Rd
% @brief roxygen2 documentation for the R Interface of R_ExaGeoStatClimateEmulator function
% @version 2.0.0
% @author Mahmoud ElKarargy
% @date 2025-11-25

\name{climate_emulator}
\alias{climate_emulator}
\title{Climate Emulator for Statistical Modeling}

\usage{
climate_emulator(N, dts, timeslot, objects_number, data_path, cores = 1, 
gpus = 0, add_diagonal = 10.0, accuracy = 0, band_dense_dp = 1000, 
hnb = 300, verbose = "detailed", log_path = "")
}

\arguments{
  \item{N}{A numeric value specifying the spatial problem size (typically dts²).}
  \item{dts}{A numeric value specifying the tile size for matrix computations.}
  \item{timeslot}{A numeric value specifying the number of time points to process (z_*.csv files).}
  \item{objects_number}{A numeric value specifying the number of objects for parallel processing.}
  \item{data_path}{A string specifying the directory containing Mean Trend Removal outputs (z_*.csv and params.csv).}
  \item{cores}{A numeric value specifying the number of CPU cores to use. Default is 1.}
  \item{gpus}{A numeric value specifying the number of GPUs to use. Default is 0.}
  \item{add_diagonal}{A numeric value to add to the diagonal for numerical stability. Default is 10.0.}
  \item{accuracy}{A numeric value specifying accuracy for low-rank approximation (10^(-accuracy)). Default is 0.}
  \item{band_dense_dp}{A numeric value specifying the band dense parameter. Default is 1000.}
  \item{hnb}{A numeric value specifying the hierarchical matrix parameter. Default is 300.}
  \item{verbose}{A string specifying verbosity level: "standard" or "detailed". Default is "detailed".}
  \item{log_path}{A string specifying the path for logging. Optional.}
}

\value{
A list containing:
  \item{completed}{Logical value (TRUE) indicating successful completion of climate emulator operations.}
}

\description{
This function runs the Climate Emulator for statistical modeling of climate data. It builds statistical
models from Mean Trend Removal residuals to emulate climate behavior using advanced matrix decomposition
techniques.

IMPORTANT: This function requires PaRSEC runtime and will not work with StarPU. Configure your project
with \code{--use-parsec} flag during build.
}

\details{
The Climate Emulator:
\itemize{
  \item Reads normalized residuals (z_*.csv files) from Mean Trend Removal
  \item Builds spatial covariance matrices for each time point
  \item Applies Tile Low-Rank (TLR) compression using HiCMA
  \item Performs matrix operations (SYRK, Cholesky decomposition)
  \item Computes statistical metrics and diagnostics
  \item Optionally performs inverse spherical harmonics transform
}

The emulator uses hierarchical matrix representations for efficient computation on large-scale
climate data. It supports both CPU and GPU acceleration.
}

\section{Prerequisites}{
Before running Climate Emulator, you must:
\itemize{
  \item Run Mean Trend Removal to generate z_*.csv and params.csv files
  \item Ensure auxiliary files exist: \code{<dts>_Et1.csv}, \code{<dts>_Et2.csv}, \code{<dts>_Ep.csv}
  \item Build ExaGeoStat with PaRSEC runtime: \code{./configure -e --climate-emulator --use-parsec}
  \item Install ExaGeoStat R package: \code{R CMD INSTALL . --configure-args="-r --use-parsec"}
}
}

\section{Input Data Format}{
The function expects the following files in \code{data_path}:
\itemize{
  \item \code{z_0000.csv, z_0001.csv, ...}: Residual data (one file per time point)
  \item \code{params.csv}: Optimized parameters from Mean Trend Removal
  \item \code{<dts>_Et1.csv}: Auxiliary matrix Et1
  \item \code{<dts>_Et2.csv}: Auxiliary matrix Et2
  \item \code{<dts>_Ep.csv}: Auxiliary matrix Ep
}
}

\section{Problem Size}{
The relationship between parameters:
\itemize{
  \item N = dts² (e.g., dts=72 → N=5184)
  \item timeslot ≤ number of z_*.csv files available
  \item objects_number typically equals dts for optimal parallelization
}
}

\section{Performance Tuning}{
Key parameters for performance optimization:
\itemize{
  \item \code{accuracy}: Lower values (higher accuracy) increase computation time
  \item \code{band_dense_dp}: Controls band-dense switching in TLR
  \item \code{hnb}: Hierarchical block size affects memory usage
  \item \code{add_diagonal}: Increases numerical stability but may affect accuracy
}
}

\examples{
\dontrun{
# Example 1: Basic Climate Emulator run
# Requires: Mean Trend Removal outputs + auxiliary files + PaRSEC runtime

library(ExaGeoStatCPP)

result <- climate_emulator(
  N = 5184,              # 72 x 72 grid
  dts = 72,
  timeslot = 4000,       # Process 4000 time points
  objects_number = 72,
  data_path = "/path/to/results/",
  cores = 8,
  gpus = 0,
  add_diagonal = 10,
  accuracy = 0,
  band_dense_dp = 1000,
  hnb = 300,
  verbose = "detailed"
)

if (result$completed) {
  cat("Climate Emulator completed successfully!\n")
}

# Example 2: High-accuracy run with GPU acceleration
result_gpu <- climate_emulator(
  N = 10368,             # 72 x 144 grid
  dts = 72,
  timeslot = 8000,
  objects_number = 72,
  data_path = "/path/to/results/",
  cores = 16,
  gpus = 2,
  add_diagonal = 5,
  accuracy = 3,          # 10^(-3) accuracy
  band_dense_dp = 500,
  hnb = 200,
  verbose = "detailed",
  log_path = "/path/to/emulator.log"
)

# Example 3: Full pipeline (Mean Trend Removal → Climate Emulator)
# Step 1: Mean Trend Removal
mtr_result <- mean_trend_removal(
  lon = 144,
  startyear = 2020,
  endyear = 2020,
  dts = 720,
  data_path = "/data/ERA_data/",
  forcing_data_path = "/data/forcing_new.csv",
  results_path = "/results/",
  starting_theta = c(0.9),
  lb = c(0.001),
  ub = c(0.95),
  cores = 4
)

# Step 2: Climate Emulator (after MTR completes)
ce_result <- climate_emulator(
  N = 5184,
  dts = 72,
  timeslot = 4000,
  objects_number = 72,
  data_path = mtr_result$results_path,  # Use MTR output
  cores = 8
)
}
}

\seealso{
\code{\link{mean_trend_removal}} for the preprocessing step that generates input data.
}

\section{Troubleshooting}{
Common issues and solutions:
\itemize{
  \item "Climate Emulator is NOT supported with StarPU": Rebuild with \code{--use-parsec}
  \item "z_XXXX.csv not found": Ensure Mean Trend Removal completed and timeslot ≤ available files
  \item "Numerical instability": Increase \code{add_diagonal} parameter
  \item Out of memory: Reduce \code{hnb} or increase available RAM
}
}

