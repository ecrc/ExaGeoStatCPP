# ExaGeoStatCPP User Manual

## Content

1. [ExaGeoStatCPP v1.1.0](#ExaGeoStatCPP)
2. [Configurations of the software](#configurations)
3. [Building ExaGeoStatCPP](#building)
4. [Arguments](#arguments)
5. [List of Descriptors](#list-of-descriptors)
6. [Supported operations](#supported-operations)
7. [Manuals](#Manuals)
8. [Contributing](#contributing)

## ExaGeoStatCPP
> Current Version of ExaGeoStatCPP: 1.1.0
### Supported Operations:
1. (Data Generation): Generating large geospatial synthetic datasets using  dense, Diagonal Super-Tile (DST) and Tile Low-Rank (TLR) approximation techniques.
2. (Data Modeling): Modeling large geospatial datasets on dense, DST and TLR approximation techniques through the Maximum likelihood Estimation (MLE) operation.
3. (Data Prediction): Predicting missing measurements on given locations using dense, DST, and TLR approximation techniques.
4. (MLOE/MMOM): Computing the Mean Loss of Efficiency (MLOE), Mean Misspecification of the Mean Square Error (MMOM), and Root mean square MOM (RMOM) to describe the prediction performance over the whole observation region.
5. (Fisher Information Matrix (FIM)): Quantifying the information content that a variable x carries about a parameter $\theta$ within a Gaussian distribution.

### Supported Covariance Functions:
1. Univariate Matérn (Gaussian/Stationary)
2. Univariate Matérn with Nugget (Gaussian/Stationary)
3. Flexible Bivariate Matérn (Gaussian/Stationary)
4. Parsimonious Bivariate Matérn (Gaussian/Stationary)
5. Parsimonious trivariate Matérn (Gaussian/Stationary)
6. Univariate Space/Time Matérn (Gaussian/Stationary)
7. Bivariate Space/Time Matérn (Gaussian/Stationary)
8. Tukey g-and-h Univariate Matérn (non-Gaussian/Stationary)
9. Tukey g-and-h Univariate Power Exponential (non-Gaussian/Stationary)
> To add your kernel, please refer to [Contribution Guidelines](CONTRIBUTING.md)

### Programming models:
1. MPI
2. Task-based programming models

### External libraries:
1. NLOPT [https://nlopt.readthedocs.io/en/latest/](https://nlopt.readthedocs.io/en/latest/)
2. GSL [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)
3. HWLOC [https://www.open-mpi.org/projects/hwloc/](https://www.open-mpi.org/projects/hwloc/)
4. StarPU dynamic runtime system  [https://starpu.gitlabpages.inria.fr/](https://starpu.gitlabpages.inria.fr/)
5. HCORE [https://github.com/ecrc/hcore](https://github.com/ecrc/hcore)
6. HiCMA [https://github.com/ecrc/hicma](https://github.com/ecrc/hicma)
7. Stars-H [https://github.com/ecrc/stars-h](https://github.com/ecrc/stars-h)
8. Chameleon [https://gitlab.inria.fr/solverstack/chameleon](https://gitlab.inria.fr/solverstack/chameleon)

### Project Hierarchy

* **```cmake```** A directory contains essential CMake modules that facilitate the importation and location of required dependencies.
* **```docs```** A directory contains all the necessary documents.
* **```examples```** A directory contains a comprehensive collection of demo code that illustrates the framework's application and demonstrates its features and capabilities.
* **```inst```** A directory contains all the system's header files, mirroring the structure of the src directory.
* **```man```** A directory contains all the R functions documentation.
* **```scripts```** A directory contains benchmarking scripts.
* **```src```** A directory contains all the source files.
* **```tests```** A directory contains all the test files and follows the same structure as the src folder.
* **```clean_build.sh```** A script is designed to compile the software tests once all dependencies are installed, and it is set to build everything by default.
* **```CMakeLists.txt```** The top-level CMake file to configure the build system.
* **```configure```** A Script used to generate the building system inside a 'bin' directory.


## Configurations

* Run the help of `configure` to know the needed arguments for your specific options.
  ``` bash
  ./configure -h
  ```

* To enable R interface, add `-r` disabled by default.
* To enable support of HiCMA, add `-H` disabled by default.
* To enable examples, add `-e` enabled by default.
* To enable tests, add `-t` disabled by default.
* To enable heavy tests, add `-T` disabled by default.
* To enable CUDA, add `-c` disabled by default.
* To enable MPI, add `-m` disabled by default.
* To enable verbose output, add `-v` disabled by default.
* To change the installation path of the dependencies, use `-i <installation/path>` the default is project_path/installdir/_deps/ on Unix systems.
* To enable packaging system for distribution, add `-p` disabled by default.
* To enable showing code warnings, add `-w` disabled by default.
* To manually set mkl as blas vendor, add `--use-mkl`. MKL is required as blas vendor and it's automatically detected but in some environments it need to be manually set.
* To enable PaRSEC as a runtime system, add `--use=parsec`, StarPU by default.

## Building

* Run the help of `clean_build.sh` to know additional argument options.

  ```bash
  ./clean_build.sh -h
  ```
* Run clean_build.sh to build the project.
  ```bash
  ./clean_build.sh
  ```
* To enable the installation of the project, add `-i` disabled by default.
* To enable verbose printing, add `-v` disabled by default.
* To enable building with a specific number of threads, add `-j <thread_number>` running with maximum number of threads by default.


## Arguments
These are the arguments that you can specify when running any C++ example.
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
### Initialize and Finalize the Hardware class

To use any operations, you must initialize the hardware by selecting the number of CPUs and/or GPUs.
```c++
// Initialize an instance of the hardware
auto hardware = ExaGeoStatHardware(computation, number of cores, number of gpus, p, q);

// Other code goes here

// Finalize the hardware instance.
hardware.FinalizeHardware()
```
The subsequent arguments are as follows:

 - `computation`: Specifies the computation mode for the solver.
 - `number of cores`: Indicates the number of CPU cores to be used for the solver.
 - `number of gpus`: Specifies the number of GPUs to be used for the solver.


##### *ExaGeoStat R Interface*
```R
hardware <- new(Hardware, computation, number of cores, number of gpus, p-grid, q-grid);
hardware$finalize_hardware()
```
First arguement represents the name of the R class that wrapps its correponding C++ class, the rest of arguments are the same as the C++ version

### Types of Data

ExaGeoStatCPP can be used with 2 types of data:

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

Here we use existing data by providing the path to it:

- The Data Path must be passed to Configuration
  ```
  data_path <- <path/to/file>
  ```

And then using the following code:

```R
exageostat_data <- simulate_data(kernel=kernel, initial_theta=initial_theta, problem_size=problem_size, dts=dts, dimension=dimension, data_path=data_path)
```

### Location Getters

ExaGeoStat support 2D and 3D spatial locations, and therefore we have getters for X, Y and Z coordinates.

#### X-coordinate Getter
```c++
double *locations_x = exageostat_data->GetLocations()->GetLocationX();
```

#### Y-coordinate Getter
```c++
double *locations_y = exageostat_data->GetLocations()->GetLocationY();
```

#### Z-coordinate Getter
```c++
double *locations_z = exageostat_data->GetLocations()->GetLocationZ();
```

##### *ExaGeoStat R Interface*
```R 
locations <- get_locations(data=exageostat_data)
```
Returns all the locations values. The subsequent arguments are as follows:

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

##### *ExaGeoStat R Interface*
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


##### *ExaGeoStat R Interface*
```R
estimated_theta <- model_data(matrix=z_value, x=locations_x, y=locations_y, kernel=kernel, dts=dts, dimension=dimension,lb=lower_bound, ub=upper_bound, mle_itr=10, computation=computation, band=1)
```
Or
```R
estimated_theta <- model_data(data=exageostat_data, kernel=kernel, dts=dts, dimension=dimension,lb=lower_bound, ub=upper_bound, mle_itr=10)
```

### Data Prediction
```c++
//You have to pass your arguments through the configurations, your hardware, and your data.
ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_matrix);
```


##### *ExaGeoStat R Interface*
```R
predict_data(train_data=list(locations_x, locations_y, locations_z, z_value), test_data=list(test_x, test_y, test_z), kernel=kernel, dts=dts, estimated_theta=estimated_theta)
```

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


##### *ExaGeoStat R Interface*
```R
fisher_matrix <- fisher(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta)
```

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


##### *ExaGeoStat R Interface*
```R
result_mloe_mmom = mloe_mmom(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta, true_theta=true_theta)
```

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


##### *ExaGeoStat R Interface*
```R
idw_error = idw(train_data=list(locations_x, locations_y, z_value), test_data=list(test_x, test_y), kernel=kernel, dts=dts, estimated_theta=estimated_theta, test_measurements=test_measurements)
```

## Manuals
- Find a detailed Manual for R functions in [ExaGeoStatCPP-R-Interface-Manual](docs/ExaGeoStat-R-Interface-Manual.pdf)
- Find a detailed Manual for C++ functions in [ExaGeoStatC-CPP-Manual](docs/ExaGeoStat-CPP-Manual.pdf)
- Doxygen Manual: https://ecrc.github.io/ExaGeoStatCPP

## Contributing
[Contribution Guidelines](CONTRIBUTING.md)
