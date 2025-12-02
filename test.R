.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1")
library(ExaGeoStatCPP)

ncores <- 4
ngpus <- 0
problem_size <- 1600
dts <- 320
lts <- 0
computation <- "exact"
dimension <- "2D"
kernel <- "UnivariateMaternNuggetsStationary"
initial_theta <- c(1,0.1,0.5,0.1)
lower_bound <- c(0.01,0.2,0.01,0.01)
upper_bound <- c(5,5,5,5)
p <- 1
q <- 1
opt_itrs <- 50
acc <- 1e-9

hardware <- new(Hardware, computation, ncores, ngpus, p, q)

exageostat_data <- simulate_data(
  kernel = kernel,
  initial_theta = initial_theta,
  problem_size = problem_size,
  dts = dts,
  dimension = dimension
)

estimated_theta <- model_data(
  matrix=exageostat_data$m,
  x=exageostat_data$x,
  y=exageostat_data$y,
  kernel=kernel, dts=dts,
  dimension=dimension,
  lb=lower_bound,
  ub=upper_bound,
  mle_itr=opt_itrs)

test_x <- c(0.2, 0.330)
test_y <- c(0.104, 0.14)
test_z <- c(-0.10838, -0.10838)

predict_data(
  kernel=kernel,
  estimated_theta=estimated_theta,
  dts=dts,
  train_data=list(exageostat_data$x, exageostat_data$y, exageostat_data$m),
  test_data=list(test_x, test_y),
  test_measurements=test_z
  )
