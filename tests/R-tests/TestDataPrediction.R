
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST)#

# @file TestDataPrediction.R
# @brief
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-09

library("ExaGeoStatCPP")

paste("---------------------------------------------------------------")
paste("ExaGeoStat with data Prediction only - mspe")

dimension = "2D"
ncores <- 2
ngpus <- 0
problem_size <- 16
dts <- 8
lts <- 0
computation <- "exact"
kernel <- "univariate_matern_stationary"
estimated_theta <- c(1,0.1,0.5)

hardware <- new(Hardware, computation, ncores, ngpus)

z_value <- c(-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
              -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
              0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
              -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
              0.971213790000161170, 0.538973474182433021)

locations_x <- c(0.092042420080872822,  0.193041886015106440,  0.330556191348134576,
                  0.181612878614480805,  0.370473792629892440, 0.652140077821011688,
                  0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                  0.465445944914930965,  0.528267338063630132,  0.974792095826657490,
                  0.552452887769893985, 0.877592126344701295)

locations_y <- c(0.928648813611047563, 0.103883421072709245,  0.135790035858701447,  0.434683756771190977,
                 0.400778210116731537,  0.168459601739528508, 0.105195696955825133,
                 0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                 0.627679865720607300,  0.958236057068741931,
                 0.573571374074921758, 0.568657969024185528)

test_x <- c(0.347951, 0.62768)
test_y <- c(0.806332, 0.105196)
predict_data(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta)

paste("---------------------------------------------------------------")
paste("ExaGeoStat with data Prediction only - idw")

test_measurements = c(-1.05428, -1.47441)
idw_error = idw(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta, test_measurements=test_measurements)
paste("idw error values:")
paste(idw_error)

paste("---------------------------------------------------------------")
paste("ExaGeoStat with data Prediction only - fisher")

fisher_matrix <- fisher(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta)
paste("fisher matrix values:")
paste(fisher_matrix)

paste("---------------------------------------------------------------")
paste("ExaGeoStat with data Prediction only - MLOE-MMOM")
true_theta <- c(1.1,0.2,0.5)

result_mloe_mmom = mloe_mmom(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta, true_theta=true_theta)
cat(sprintf("MLOE = %.6f", result_mloe_mmom[1]), "\n")
cat(sprintf("MMOM = %.6f", result_mloe_mmom[2]), "\n")
