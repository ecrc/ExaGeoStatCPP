
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
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

dimension = "3D"
ncores <- 4
ngpus <- 0
problem_size <- 16
dts <- 8
lts <- 0
computation <- "exact"

save_data <- TRUE
# if you didn't change the log_path it will be saved in the default path of project_path/synthetic_ds
log_path <- ""
# data path is where to read data from
data_path <- ""
# observations file path is where to read observation file
observations_file <- ""
# recovery file path is where to read recovery file
recovery_file <- ""

hardware <- new(Hardware, computation, ncores, ngpus)
config <- configurations_init(n=problem_size, cores_gpus=c(ncores, ngpus), kernel="univariate_matern_stationary", computation=computation, tile_size=c(dts,lts), iTheta=c(1,0.1,0.5), lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=5, prediction=c(5,1,1,1,1), paths=c(log_path, data_path, observations_file, recovery_file), save_data=save_data, dimension=dimension)
data_source <- new(Data, problem_size, dimension)

exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)
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

save_data <- FALSE
log_path <- ""
# data path is where to read data from
data_path <- "./synthetic_ds/SYN_16_1"
# observations file path is where to read observation file
observations_file <- ""
# recovery file path is where to read recovery file
recovery_file <- ""

config <- configurations_init(n=problem_size, cores_gpus=c(ncores, ngpus), kernel="univariate_matern_stationary", computation=computation, tile_size=c(dts,lts), iTheta=c(1,0.1,0.5), lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=5, prediction=c(5,1,1,1,1), paths=c(log_path, data_path, observations_file, recovery_file), save_data=save_data, dimension=dimension)
data_source <- new(Data, problem_size, dimension)
exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)

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