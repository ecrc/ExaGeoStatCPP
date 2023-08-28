ExaGeoStat User Manual
================

# Content

1. Configurations of the software.
2. Building ExaGeoStat.
3. Supported Covariance kernels.
4. Arguments.
5. List of Descriptors.
6. Supported operations.
7. Contributing

## Configurations

* Run help of config.sh to know the needed arguments to run with your specific options.
```commandline
./config.sh -h
```
* To Enable support of Chameleon, Use the following argument.
```commandline
./config.sh -C
```
* To Enable support of HiCMA, Use the following argument.
```commandline
./config.sh -H
```
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
```commandline
./clean_build.sh -h
```
* Run clean_build.sh to clean, Build and Install project.
```commandline
./clean_build.sh
```
* To enable verbose printing Run the following command.
```commandline
./clean_build.sh -v
```
* To enable building with specific number of threads Run the following command.
```commandline
./clean_build.sh -j <thread_number>
```

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

## List of Descriptors
### Covariance Matrix Descriptors
1. DESCRIPTOR_C: Covariance matrix C descriptor
2. DESCRIPTOR_C11: Covariance Matrix C11 descriptor.
3. DESCRIPTOR_C12: Covariance Matrix C21 descriptor.
4. DESCRIPTOR_C21: Covariance Matrix C12 descriptor.
5. DESCRIPTOR_C22: Covariance Matrix C22 descriptor.
### Measurements Matrix Descriptors
1. DESCRIPTOR_Z: Measurements Z descriptor.
2. DESCRIPTOR_Z_COPY: A copy of Measurements Z descriptor.
3. DESCRIPTOR_Z_OBSERVATIONS: Observed Measurements Z descriptor.
4. DESCRIPTOR_Z_Actual: Actual Measurements Z descriptor.
### Measurements Sub-Matrix Descriptors
1. DESCRIPTOR_Z_1: Measurements Z1 sub-matrix descriptor.
2. DESCRIPTOR_Z_2: Measurements Z2 sub-matrix descriptor.
### Dot Product Descriptors
1. DESCRIPTOR_PRODUCT: Dot product descriptor.
2. DESCRIPTOR_PRODUCT_1: Dot product descriptor.
3. DESCRIPTOR_PRODUCT_2: Dot product descriptor.
### Determinant Descriptors
1. DESCRIPTOR_DETERMINANT: Determinant descriptor.
### Mean Square Error Descriptors
1. DESCRIPTOR_MSE: Mean Square Error descriptor.
### HiCMA Descriptors
1. DESCRIPTOR_CRK: HiCMA descCrk descriptor.
2. DESCRIPTOR_C12RK: HiCMA descCrk descriptor.
3. DESCRIPTOR_C22RK: HiCMA descCrk descriptor.
4. DESCRIPTOR_CD: HiCMA descCD descriptor.
5. DESCRIPTOR_C12D: HiCMA descCD descriptor.
6. DESCRIPTOR_C22D: HiCMA descCD descriptor.
7. DESCRIPTOR_CUV : HiCMA descCUV descriptor
8. DESCRIPTOR_C12UV : HiCMA descCUV descriptor
9. DESCRIPTOR_C22UV : HiCMA descCUV descriptor



## Supported operations

### Provide your arguments

To use any operations, you firstly need to provide your arguments to the operation.
This is done by the Configurations module. You have two ways to set your arguments.

1. Provide your arguments with the command line.

    // Create a new configurations object.
    Configurations configurations;
    // Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);

2. Set your arguments manually in the code.

    Configurations synthetic_data_configurations;
    synthetic_data_configurations.SetProblemSize(10);
    synthetic_data_configurations.SetKernelName("BivariateSpacetimeMaternStationary");
    synthetic_data_configurations.SetPrecision(exageostat::common::double);

### Initialise the Hardware

To use any operations you need to initialise the hardware with your selection of number of cores and gpus.
```c++
auto hardware = ExaGeoStatHardware(computation, number of cores, number of gpus);
```

### Data Generation

ExaGeoStat can be used with different types of data, including:

Synthetic data: Synthetic data is generated by the software according to the user arguments.

Real data: ExaGeoStat can also be used with real data, such as data from satellite imagery or weather sensors. Real data
can be used to train the software to better predict the values of new data.

##### Synthetic data generation

After providing your arguments with the Configurations module, you need to do the following two steps.
```c++
// Create a new ExaGeoStat data that holds the locations data and descriptors data.
ExaGeoStatData<double> data(configurations.GetProblemSize(), configurations.GetDimension(), hardware);
// Generate data by passing your arguments throuh the configurations, your hardware and your container of the data which will be filled with the new generated data.
ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, configurations, data);
```

##### Real data generation

- Not supported for now.

### Data Modeling.

To use data modeling you have to do this operation.
```c++
// you have to pass your arguments throuh the configurations, your hardware and your data.
ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data);
```

### Data prediction

- Not supported for now.
