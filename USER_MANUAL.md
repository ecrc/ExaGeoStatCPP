ExaGeoStat User Manual
================
## Configurations
* Run help of config.sh to know the needed arguments to run with your specific options.

        ./config.sh -h
* To Enable support of HiCMA, Use the following argument.

        ./config.sh -H
* To Enable support of Chameleon, Use the following argument.

        ./config.sh -C
* To enable examples add **```-e```** 

* To enable tests add **```-t```** 

* To enable CUDA add **```-c```**

* To enable MPI add **```-m```**

* To enable Verbose add **```-v```**

* To enable debug mode add **```-d```**

* To change the installation path of the dependencies use **```-i <installation/path>```**

Please note, currently we are supporting using only one of the two options of HiCMA or Chameleon.

## Building

* Run help of clean_build.sh to know the additional arguments options.

        ./clean_build.sh -h
* Run clean_build.sh to clean, Build and Install project.

        ./clean_build.sh
* To enable verbose printing Run the following command.

        ./clean_build.sh -v
* To enable building with specific number of threads Run the following command.

        ./clean_build.sh -j <thread_number>

Supported Covariance Functions/ Kernels:
======================
1. univariate_matern_stationary
2. univariate_exp_non_gaussian
3. univariate_matern_dbeta
4. univariate_matern_ddbeta_beta
5. univariate_matern_ddbeta_nu
6. univariate_matern_ddnu_nu
7. univariate_matern_ddsigma_square
8. univariate_matern_ddsigma_square_beta
9. univariate_matern_ddsigma_square_nu
10. univariate_matern_dnu
11. univariate_matern_dsigma_square
12. univariate_matern_non_gaussian
13. univariate_matern_non_sta
14. univariate_matern_non_stationary
15. univariate_matern_nuggets_stationary
16. univariate_spacetime_matern_stationary
17. bivariate_matern_flexible
18. bivariate_matern_parsimonious
19. bivariate_spacetime_matern_stationary
20. Trivariate_matern_parsimonious


## Arguments

* {Mandatory} To set the problem size (N)

        --N=<problem_size>
* {Mandatory} To set the kernel

        --kernel=<supported_kernel>
* {Mandatory} To set the dense tile size in case of using Chameleon

        --dts=<value>
* {Mandatory} To set the low tile size in case of using HiCMA

        --lts=<value>
* {Optional} To set the dimension, Default is 2D

        --dimension=<2D/3D/ST>
* {Optional} To set the p grid

        --p_grid=<value>
* {Optional} To set the q grid

        --q_grid=<value>
* {Optional} To set the time slot, Default is 1

        --time_slot=<value>
* {Optional} To set the computation, Default is dense

        --computation=<dense/tlr/dst>
* {Optional} To set the precision, Default is double

        --precision=<single/double/mix>
* {Optional} To set the number of cores, Default is 1

        --cores=<value>
* {Optional} To set the number of GPUs, Default is 0

        --gpus=<value>
* {Optional} To set the number of unknown observation to be predicted, Default is 0

        --Zmiss=<value>
* {Optional} To set the path of observation file

        --observations_file=<path/to/file>
* {Optional} To set the max rank, Default is 1

        --max_rank=<value>
* {Optional} To set the lower bounds of optimization

        --olb=<value:value:....:value>
* {Optional} To set the upper bounds of optimization

        --oub=<value:value:....:value>
* {Optional} To set the initial theta

        --itheta=<value:value:....:value>
* {Optional} To set the target theta

        --ttheta=<value:value:....:value>
* {Optional} To set the seed value, Default is 0

        --seed=<value>
* {Optional} To set the run mode value, Default is standard

        --run_mode=<verbose/standard>
* {Optional} To set the path of log files to be written, Default is ./exageostat-cpp/synthetic_ds/

        --log_path=<path/to/file>
* {Optional} To enable generating of synthetic data, Default is ON

        --synthetic_data
* {Optional} To enable Out of core technology, Default is OFF

        --OOC 
* {Optional} To enable approximation mode, Default is ON

        --approximation_mode
* {Optional} To enable writing log files, Default is OFF

        --log