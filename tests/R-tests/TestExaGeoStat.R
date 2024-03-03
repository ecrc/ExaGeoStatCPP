
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST)#

# @file TestExaGeoStat.R
# @brief
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-09
#

paste("---------------------------------------------------------------")
paste("ExaGeoStat with all Modules - Dense")
library("ExaGeoStatCPP")

ncores <- 1
ngpus <- 0
problem_size <- 16
dts <- 8
lts <- 0
computation <- "exact"
dimension = "2D"

hardware <- new(Hardware, computation, ncores, ngpus)
config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(dts,lts), iTheta=c(1,0.1,0.5), lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=5, prediction=c(6,1,1,1,1), dimension=dimension)
data_source <- new(Data, problem_size, dimension)

exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)
model_data(hardware=hardware, config=config, data=exageostat_data)
predict_data(hardware=hardware, config=config, data=exageostat_data)

paste("---------------------------------------------------------------------------------------------")
paste("ExaGeoStat with Data Generation - saving data with default path")

dimension = "3D"
save_data <- TRUE
# if you didn't change the log_path it will be saved in the default path of project_path/synthetic_ds
log_path <- ""
# data path is where to read data from
data_path <- ""
# observations file path is where to read observation file
observations_file <- ""
# recovery file path is where to read recovery file
recovery_file <- ""
config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(dts,lts), iTheta=c(1,0.1,0.5), lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=5, prediction=c(5,1,1,1,1), paths=c(log_path, data_path, observations_file, recovery_file), save_data=save_data, dimension=dimension)
data_source <- new(Data, problem_size, dimension)
exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)
## Print the data..
#paste("Locations x")
#x <- get_locationsX(data=exageostat_data)
#print(x)
#paste("Locations y")
#y <- get_locationsY(data=exageostat_data)
#print(y)
paste("Locations z")
z <- get_locationsZ(data=exageostat_data)
print(z)
#paste("descZ measurements values")
#Z <- get_Z_measurement_vector(data=exageostat_data, type="chameleon")
#print(Z)


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

config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(dts,lts), iTheta=c(1,0.1,0.5), lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=5, prediction=c(5,1,1,1,1), paths=c(log_path, data_path, observations_file, recovery_file), save_data=save_data, dimension=dimension)
data_source <- new(Data, problem_size, dimension)
exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)


paste("---------------------------------------------------------------------------------------------")

paste("ExaGeoStat with data Modeling only")

config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(8,lts), iTheta=c(1,0.1,0.5), max_rank = 500, lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=10)
exageostat_data <- new(Data, problem_size, "2D")
numbers <- c( -1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
             -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
             0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
             -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
             0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
             0.290822066007430102)
z_value = matrix(numbers, nrow = 1, ncol = 16, byrow = TRUE)

numbers <- c(0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
              0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
              0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
              0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
              0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
              0.877592126344701295)
locations_x = matrix(numbers, nrow = 1, ncol = 16, byrow = TRUE)

numbers <- c(0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
            0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
            0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
            0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
            0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
            0.942824444953078489)
locations_y = matrix(numbers, nrow = 1, ncol = 16, byrow = TRUE)
model_data(hardware=hardware, config=config, data=exageostat_data, matrix=z_value, x=locations_x, y=locations_y)

paste("---------------------------------------------------------------")
paste("ExaGeoStat with data Prediction only - all prediction functions")
#config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(8,lts), eTheta=c(0.9, 0.09, 0.4), iTheta=c(1,0.1,0.5), max_rank = 500, lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), prediction=c(5,1,1,1,1))
#predict_data(hardware=hardware, config=config, data=exageostat_data, matrix=z_value, x=locations_x, y=locations_y)
