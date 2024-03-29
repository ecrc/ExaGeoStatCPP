ExaGeoStatCPP User Manual
=========================

# Content

1. [Configurations of the software](#configurations)
2. [Building ExaGeoStatCPP](#building)
3. [Supported Covariance kernels](#supported-covariance-functions-kernels-)
4. [Arguments](#arguments)
5. [List of Descriptors](#list-of-descriptors)
6. [Supported operations](#supported-operations)
7. [Contributing](#contributing)

## Configurations

* Run the help of `config.sh` to know the needed arguments for your specific options.

```commandline
./config.sh -h
```

* To Enable support of HiCMA, add `-H` disabled by default.
* To enable examples, add `-e` enabled by default.
* To enable tests, add `-t` disabled by default.
* To enable heavy tests, add `-T` disabled by default.
* To enable CUDA, add `-c` disabled by default.
* To enable MPI, add `-m` disabled by default.
* To enable verbose output, add `-v` disabled by default.
* To change the installation path of the dependencies, use `-i <installation/path>` project_path/installdir/_deps/ by default on Unix systems.
* To enable manually passing mkl as BLA vendor, add `--use-mkl` MKL by default.
* To enable packaging system for distribution, add `-p` disabled by default.

## Building

* Run the help of `clean_build.sh` to know additional argument options.

```commandline
./clean_build.sh -h
```
* Run clean_build.sh to build the project.
```commandline
./clean_build.sh
```
* To build and install the project, Run the following command.
```commandline
./clean_build.sh -i
```
* To build and enable verbose printing, Run the following command.
```commandline
./clean_build.sh -v
```
* To enable building with a specific number of threads, run the following command.
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
13. univariate_matern_nuggets_stationary
14. univariate_spacetime_matern_stationary
15. bivariate_matern_flexible
16. bivariate_matern_parsimonious
17. bivariate_spacetime_matern_stationary
18. trivariate_matern_parsimonious

* To add your kernel, please refer to [Contribution Guidelines](CONTRIBUTING.md)
## Arguments

* {Mandatory} To set the problem size (N)

        --N=<problem_size>
* {Mandatory} To set the kernel

        --kernel=<supported_kernel>
* {Mandatory} To set the dense tile size in the case of Chameleon

        --dts=<value>
* {Mandatory} To set the low tile size in case of HiCMA

        --lts=<value>
* {Optional} To set the dimension, the default is 2D

        --dimension=<2D/3D/ST>
* {Optional} To set the p grid, the default is 1

        --p=<value>
* {Optional} To set the q grid, the default is 1

        --q=<value>
* {Optional} To set the time slot, the default is 1

        --time_slot=<value>
* {Optional} To set the computation, the default is dense

        --computation=<dense/tlr/dst>
* {Optional} To set the precision, the default is double

        --precision=<single/double>
* {Optional} To set the number of cores, the default is 1

        --cores=<value>
* {Optional} To set the number of GPUs, the default is 0

        --gpus=<value>
* {Optional} To set the number of unknown observations to be predicted, the default is 0

        --Zmiss=<value>
* {Optional} To set the path of the observation file

        --observations_file=<path/to/file>
* {Optional} To set the max rank, the default is 1

        --max_rank=<value>
* {Optional} To set the lower bounds of optimization

        --olb=<value:value:....:value>
* {Optional} To set the upper bounds of optimization

        --oub=<value:value:....:value>
* {Optional} To set the initial theta

        --itheta=<value:value:....:value>
* {Optional} To set the target theta

        --ttheta=<value:value:....:value>
* {Optional} To set the estimated theta

        --etheta=<value:value:....:value>
* {Optional} To set the seed value, the default is 0

        --seed=<value>
* {Optional} To set the verbose value, the default is standard

        --verbose=<quiet/standard/detailed>
* {Optional} To set the path of log files to be written, the default is ./exageostat-cpp/synthetic_ds/

        --log_path=<path/to/file>
* {Optional} To enable reading a CSV file containing real data, if not entered the default is the generation of synthetic data

        --data_path=<path/to/file>
* {Optional} To enable out-of-core (OOC), the default is OFF

        --OOC 
* {Optional} To enable approximation mode, the default is ON

        --approximation_mode
* {Optional} To enable writing log files, the default is OFF

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
5. DESCRIPTOR_Z_MISS: Missing Measurements Z descriptor.
### Measurements Sub-Matrix Descriptors
1. DESCRIPTOR_Z_1: Measurements Z1 sub-matrix descriptor.
2. DESCRIPTOR_Z_2: Measurements Z2 sub-matrix descriptor.
3. DESCRIPTOR_Z_3: Measurements Z3 sub-matrix descriptor.
### Dot Product Descriptors
1. DESCRIPTOR_PRODUCT: Dot product descriptor.
2. DESCRIPTOR_PRODUCT_1: Dot product descriptor.
3. DESCRIPTOR_PRODUCT_2: Dot product descriptor.
4. DESCRIPTOR_PRODUCT_3: Dot product descriptor.
### Determinant Descriptors
1. DESCRIPTOR_DETERMINANT: Determinant descriptor.
### Mean Square Prediction Error Descriptors
1. DESCRIPTOR_MSPE: Mean Square Error descriptor.
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



## Supported Operations

### Provide Arguments

To use any operations, you must initially supply the necessary arguments 
to the operation via the Configurations module. There are two methods available
for setting your arguments:

1. Provide your arguments with the command line.
```c++
    // Create a new configuration object.
    Configurations configurations;
    // Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
```
2. Set your arguments manually in the code.
```c++
    Configurations synthetic_data_configurations;
    synthetic_data_configurations.SetProblemSize(10);
    synthetic_data_configurations.SetKernelName("BivariateSpacetimeMaternStationary");
    synthetic_data_configurations.SetPrecision(exageostat::common::double);
```
### Initialize and Finalize the Hardware

To use any operations, you must initialize the hardware by selecting the number of CPUs and/or GPUs.
```c++
// Initialize an instance of the hardware
auto hardware = ExaGeoStatHardware(computation, number of cores, number of gpus);

// Other code goes here

// Finalize the hardware instance.
hardware.FinalizeHardware()
```
The subsequent arguments are as follows:

 - `computation`: Specifies the computation mode for the solver.
 - `number of cores`: Indicates the number of CPU cores to be used for the solver.
 - `number of gpus`: Specifies the number of GPUs to be used for the solver.


##### *ExaGeoStatR wrapper*
```R
hardware <- new(Hardware, computation, number of cores, number of gpus);
hardware$finalize_hardware()
```
First arguement represents the name of the R class that wrapps its correponding C++ class, the rest of arguments are the same as the C++ version

### Types of Data

`ExaGeoStatCPP` can be used with 2 types of data:

- Synthetic data i.e. generated by the software according to the user arguments.
- Real data e.g. data from satellite imagery or weather sensors. Real data can be used to train the software to predict the values of new data better.

#### Using Synthetic Data

Here we generate the data to be used by providing the needed arguments with the Configurations module, and then using the following code:
```c++
// Create a new ExaGeoStat data that holds the locations and descriptors data.
std::unique_ptr<ExaGeoStatData<double>> data;

// Generate data by passing your arguments through the configurations, hardware, 
//and container of the data, which will be filled with the newly generated data.
ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);
```

#### Using Real Data

Here we use already existing data by providing the path to it:

- The Data Path must be passed to Configuration
```
--log_path=<path/to/file>
```

And then using the following code:

```c++
// Create a new ExaGeoStat data that holds the locations data and descriptors data.
std::unique_ptr<ExaGeoStatData<double>> data;

// Generate data by passing your arguments through the configurations, your hardware and your container of the data which will be filled with the new generated data.
ExaGeoStat<double>::ExaGeoStatLoadData(hardware, configurations, data);
```
The subsequent arguments are as follows:

 - `kernel`: Specifies the kernel function to be used in the simulation.
 - `initial_theta`: Sets the initial values for the parameters of the simulation.
 - `problem_size`: Defines the size of the problem or dataset to be simulated.
 - `dts`: Specifies the time steps or intervals for the simulation.
 - `dimension`: Indicates the dimensionality of the data to be simulated.

### Location Getters

ExaGeoStat supports locations data of dimension upto 3D, and therefore we have getters for X, Y and Z coordinates.

#### X-coordinate Getter
```c++
double *locations_x = exageostat_data->GetLocations()->GetLocationX();
```
##### *ExaGeoStatR wrapper*
```R 
locations_x <- get_locationsX(data=exageostat_data)
```
The subsequent arguments are as follows:

- `exagostat_data`: pointer to ExaGeoStatData object containing the spatial data..

#### Y-coordinate Getter
```c++
double *locations_y = exageostat_data->GetLocations()->GetLocationY();
```
##### *ExaGeoStatR wrapper*
```R 
locations_y <- get_locationsY(data=exageostat_data)
```
The subsequent arguments are as follows:

- `exagostat_data`: pointer to ExaGeoStatData object containing the spatial data..

#### Z-coordinate Getter
```c++
double *locations_z = exageostat_data->GetLocations()->GetLocationZ();
```
##### *ExaGeoStatR wrapper*
```R 
locations_x <- get_locationsZ(data=exageostat_data)
```
The subsequent arguments are as follows:

- `exagostat_data`: pointer to ExaGeoStatData object containing the spatial data.

#### Descriptive Z Values Getter
This function is used to retrieve descriptive Z values from ExaGeoStat data based on type of descriptor.
```c++
//in this example, we use a chameleon descriptor. Similar code can be used for hicma descriptor.
DescriptorType descriptor_type = CHAMELEON_DESCRIPTOR;
void *descriptor = exageostat_data->GetDescriptorData()->GetDescriptor(descriptor_type, DESCRIPTOR_Z).chameleon_desc;
double *desc_Z_values = exagostat_data->GetDescriptorData()->GetDescriptorMatrix(descriptor_type, descriptor);
```
The used variables are as follows:

- `descriptor_type`: enum denoting the descriptor type,e.g. CHAMELEON_DESCRIPTOR,HICMA_DESCRIPTOR.
- `exagostat_data`: pointer to ExaGeoStatData object containing the spatial data.
- `desc_Z_values`: pointer to descriptor matrix.

##### *ExaGeoStatR wrapper*
```R 
desc_Z_values <- get_Z_measurement_vector(data=exageostat_data, type="chameleon")
```
The subsequent arguments are as follows:

- `exagostat_data`: pointer to ExaGeoStatData object containing the spatial data.
- `type`: string specifying the type of descriptor value to retrieve.

### Data Modeling

To use data modeling, you have to do this operation.
```c++
//You have to pass your arguments through the configurations, your hardware, and your data.
ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data, z_matrix);
```


##### *ExaGeoStatR wrapper*
```R
theta <- model_data(data=empty_data, matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, dimension=dimension, lb=lower_bound, ub=upper_bound, mle_itr=10, computation=computation)
```
This function models data based on specified parameters:

 - `data`: The dataset to be used for modeling.
 - `matrix`: The matrix of values to be used in the model.
 - `x`: The x-coordinates of the data points.
 - `y`: The y-coordinates of the data points.
 - `kernel`: The kernel function to be applied in the model.
 - `dts`: The time steps or intervals for the model.
 - `dimension`: The dimensionality of the data.
 - `lb`: The lower bound for the model parameters.
 - `ub`: The upper bound for the model parameters.
 - `mle_itr`: The number of iterations for maximum likelihood estimation.
 - `computation`: The computation mode for the model.


### Data Prediction
```c++
//You have to pass your arguments through the configurations, your hardware, and your data.
ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_matrix);
```


##### *ExaGeoStatR wrapper*
```R
predict_data(data=empty_data, matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, estimated_theta=estimated_theta, z_miss=5)
```
This function predicts data based on specified parameters:

 - `data`: The dataset to be used for prediction.
 - `matrix`: The matrix of values to be used in the prediction.
 - `x`: The x-coordinates of the data points.
 - `y`: The y-coordinates of the data points.
 - `kernel`: The kernel function to be applied in the prediction.
 - `dts`: The time steps or intervals for the prediction.
 - `estimated_theta`: The estimated parameters from the model.
 - `z_miss`: The number of missing values to be handled in the prediction.

### Fisher Function
1. Pass the fisher arguments to the Configurations.
```
--fisher
```

2. Call the Data Prediction function.
```c++
// you have to pass your arguments through the configurations, your hardware and your data.
ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_matrix);
```


##### *ExaGeoStatR wrapper*
```R
fisher(data=empty_data, matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, estimated_theta=estimated_theta, z_miss=5)
```
This function predicts data based on specified parameters:

 - `data`: The dataset to be used for the Fisher method.
 - `matrix`: The matrix of values to be used in the Fisher method.
 - `x`: The x-coordinates of the data points.
 - `y`: The y-coordinates of the data points.
 - `kernel`: The kernel function to be applied in the Fisher method.
 - `dts`: The time steps or intervals for the Fisher method.
 - `estimated_theta`: The initial estimated parameters for the Fisher method.
 - `z_miss`: The number of missing values to be handled in the Fisher method.

### MLOE-MMOM Function
1. Pass the MLOE-MMOM arguments to the Configurations.

```
--mloe-mmom
```

2. Call the Data Prediction function.
```c++
// you have to pass your arguments through the configurations, your hardware and your data.
ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_matrix);
```


##### *ExaGeoStatR wrapper*
```R
mloe_mmom(data=empty_data, matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, estimated_theta=estimated_theta, true_theta=estimated_theta, z_miss=5)
```
This function predicts data based on specified parameters:

 - `data`: The dataset to be used for the MLE.
 - `matrix`: The matrix of values to be used in the MLE.
 - `x`: The x-coordinates of the data points.
 - `y`: The y-coordinates of the data points.
 - `kernel`: The kernel function to be applied in the MLE.
 - `dts`: The time steps or intervals for the MLE.
 - `estimated_theta`: The initial estimated parameters for the MLE.
 - `true_theta`: The true parameters for the MLE, used for comparison or validation.
 - `z_miss`: The number of missing values to be handled in the MLE.

### IDW Function
1. Pass the IDW arguments to the Configurations.

```
--idw
```

2. Call the Data Prediction function.
```c++
// you have to pass your arguments through the configurations, your hardware and your data.
ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_matrix);
```


##### *ExaGeoStatR wrapper*
```R
idw(data=empty_data, matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, estimated_theta=estimated_theta, z_miss=5)
```
This function predicts data based on specified parameters:

 - `data`: The dataset to be used for the IDW.
 - `matrix`: The matrix of values to be used in the IDW.
 - `x`: The x-coordinates of the data points.
 - `y`: The y-coordinates of the data points.
 - `kernel`: The kernel function to be applied in the IDW.
 - `dts`: The time steps or intervals for the IDW.
 - `estimated_theta`: The initial estimated parameters for the IDW.
 - `z_miss`: The number of missing values to be handled in the IDW.

## Contributing
[Contribution Guidelines](CONTRIBUTING.md)
 - 