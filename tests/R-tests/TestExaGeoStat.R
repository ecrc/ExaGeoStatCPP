
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

hardware <- new(Hardware, computation, ncores, ngpus)
config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(dts,lts), iTheta=c(1,0.1,0.5), lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=10, prediction=c(5,1,1,1,1))
data_source <- new(Data, problem_size, "2D")

exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)
model_data(hardware=hardware, config=config, data=exageostat_data)
predict_data(hardware=hardware, config=config, data=exageostat_data)

paste("---------------------------------------------------------------")
hardware$finalize_hardware()
paste("ExaGeoStat with data generation and Modeling only")

problem_size <- 10
ncores <- 5

hardware <- new(Hardware, computation, ncores, ngpus)
config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(8,lts), iTheta=c(1,0.1,0.5), max_rank = 500, lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), mle_itr=5)
data_source <- new(Data, problem_size, "2D")
exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)
model_data(hardware=hardware, config=config, data=exageostat_data)

paste("---------------------------------------------------------------")
paste("ExaGeoStat with data generation and Prediction using MLOE_MMOM only")

config <- configurations_init(n=problem_size, kernel="univariate_matern_stationary", computation=computation, tile_size=c(8,lts), eTheta=c(1,0.1,0.3), iTheta=c(1,0.1,0.5), max_rank = 500, lb_ub=list(c(0.1,0.1,0.1),c(5,5,5)), prediction=c(5,0,0,0,1))
data_source <- new(Data, problem_size, "2D")
exageostat_data <- simulate_data(hardware=hardware, config=config, data=data_source)
predict_data(hardware=hardware, config=config, data=exageostat_data)


