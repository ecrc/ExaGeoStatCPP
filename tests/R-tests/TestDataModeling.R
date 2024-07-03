
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
paste("ExaGeoStat with data Modeling only - dst")

dimension = "2D"
ncores <- 2
ngpus <- 0
problem_size <- 16
dts <- 8
lts <- 0
computation <- "diag_approx"
kernel <- "univariate_matern_stationary"
initial_theta <- c(1,0.1,0.5)
lower_bound <- c(0.1,0.1,0.1)
upper_bound <- c(5,5,5)
p <- 1
q <- 1

hardware <- new(Hardware, computation, ncores, ngpus, p, q)

z_value <- c( -1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
             -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
             0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
             -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
             0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
             0.290822066007430102)

locations_x <- c(0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
              0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
              0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
              0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
              0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
              0.877592126344701295)

locations_y <- c(0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
            0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
            0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
            0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
            0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
            0.942824444953078489)

theta <- model_data(matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, dimension=dimension,lb=lower_bound, ub=upper_bound, mle_itr=10, computation=computation, band=1)
hardware$finalize_hardware()

paste("---------------------------------------------------------------------------------------------")
paste("ExaGeoStat with data Modeling only - tlr")

lts <- 8
computation <- "tlr"

hardware <- new(Hardware, computation, ncores, ngpus, 1, 1)
theta <- model_data(matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, lts=lts, dimension=dimension,lb=lower_bound, ub=upper_bound, mle_itr=10, computation=computation, band=1)
hardware$finalize_hardware()

paste("---------------------------------------------------------------------------------------------")
paste("ExaGeoStat with data Modeling only - exact")

lts <- 0
computation <- "exact"

hardware <- new(Hardware, computation, ncores, ngpus, 1, 1)
theta <- model_data(matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, dimension=dimension,lb=lower_bound, ub=upper_bound, mle_itr=10, computation=computation)
