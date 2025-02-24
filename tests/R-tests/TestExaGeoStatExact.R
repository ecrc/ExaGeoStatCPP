
# Copyright (c) 2017-2025 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST)#

# @file TestExaGeoStatAPI.R
# @brief
# @version 1.1.0
# @author Sameh Abdulah
# @date 2025-02-24

library("ExaGeoStatCPP")

paste("---------------------------------------------------------------")
paste("ExaGeoStat with all Modules - Dense")

ncores <- 30
ngpus <- 0
problem_size <- 1600
dts <- 320
lts <- 0
computation <- "exact"
dimension = "2D"
kernel <- "univariate_matern_stationary"
initial_theta <- c(1,0.1,0.5)
lower_bound <- c(0.1,0.1,0.1)
upper_bound <- c(5,5,5)
p <- 1
q <- 1
opt_itrs <- 100

hardware <- new(Hardware, computation, ncores, ngpus, p, q)

exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta, problem_size=problem_size, dts=dts, dimension=dimension)
estimated_theta <- model_data(matrix=exageostat_data$m, x=exageostat_data$x, y=exageostat_data$y, kernel=kernel, dts=dts, dimension=dimension,lb=lower_bound, ub=upper_bound, mle_itr=opt_itrs)

test_x <- c(0.2, 0.330)
test_y <- c(0.104, 0.14)
predict_data(train_data=list(x=exageostat_data$x, y=exageostat_data$y, exageostat_data$m), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta)

