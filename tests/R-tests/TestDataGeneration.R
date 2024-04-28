
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST)#

# @file TestDataGeneration.R
# @brief
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-09

library("ExaGeoStatCPP")

paste("---------------------------------------------------------------------------------------------")
paste("ExaGeoStat with Data Generation only - saving data with default path")

# Variables
dimension = "3D"
ncores <- 4
ngpus <- 0
problem_size <- 16
dts <- 8
lts <- 0
computation <- "exact"
kernel <- "univariate_matern_stationary"
initial_theta <- c(1,0.1,0.5)

# You need to provide a log_path to save the data.
log_path <- getwd()
# data path is where to read data from
data_path <- ""
# observations file path is where to read observation file
observations_file <- ""
# recovery file path is where to read recovery file
recovery_file <- ""
p <- 1
q <- 1

hardware <- new(Hardware, computation, ncores, ngpus, p, q)
exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta, problem_size=problem_size, dts=dts, dimension=dimension, log_path=log_path)

# Print the data..
paste("** Locations x")
x <- get_locationsX(data=exageostat_data)
print(x)
paste("** Locations y")
y <- get_locationsY(data=exageostat_data)
print(y)
paste("** Locations z")
z <- get_locationsZ(data=exageostat_data)
print(z)
paste("** DescZ measurements values")
Z <- get_Z_measurement_vector(data=exageostat_data, type="chameleon")
print(Z)

paste("---------------------------------------------------------------------------------------------")
paste("ExaGeoStat with Data Generation - Reading Data")

# data path is where to read data from
data_path <- "./synthetic_ds/SYN_16_1"

exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta, problem_size=problem_size, dts=dts, dimension=dimension, data_path=data_path)

paste("** Locations x - after reading")
x <- get_locationsX(data=exageostat_data)
print(x)
paste("** Locations y - after reading")
y <- get_locationsY(data=exageostat_data)
print(y)
paste("** Locations z - after reading")
z <- get_locationsZ(data=exageostat_data)
print(z)
paste("** DescZ measurements values - after reading")
Z <- get_Z_measurement_vector(data=exageostat_data, type="chameleon")
print(Z)