What is ExaGeoStat?
================

The **Exascale GeoStatistics** project (ExaGeoStat) is a parallel high-performance unified framework for computational
geostatistics on many-core systems. The project aims to optimize the likelihood function for a given spatial data to
provide an efficient way to predict missing observations in the context of climate/weather forecasting applications.
This machine learning framework proposes a unified simulation code structure to target various hardware architectures,
from commodity x86 to GPU accelerator-based shared and distributed-memory systems. ExaGeoStat enables statisticians to
tackle computationally challenging scientific problems at large-scale while abstracting the hardware complexity through
state-of-the-art high-performance linear algebra software libraries.


What is ExaGeoStatCPP?
====================
ExaGeoStatCPP is a C++ API for ExaGeoStat that aims to offer a user-friendly and efficient API for C++ developers, essentially maintaining traditional practices but also embracing contemporary C++ elements like namespaces, templates, and exceptions to enhance functionality.



The **Exascale GeoStatistics** project (ExaGeoStat) is a parallel high-performance unified framework for computational
geostatistics on many-core systems. The project aims to optimize the likelihood function for a given spatial data to
provide an efficient way to predict missing observations in the context of climate/weather forecasting applications.
This machine learning framework proposes a unified simulation code structure to target various hardware architectures,
from commodity x86 to GPU accelerator-based shared and distributed-memory systems. ExaGeoStat enables statisticians to
tackle computationally challenging scientific problems at large-scale while abstracting the hardware complexity through
state-of-the-art high-performance linear algebra software libraries.


Vision of ExaGeoStat/ExaGeoStatCPP
=================

ExaGeoStat/ExaGeoStatCPP is a collaboration between the KAUST Spatial Statistics group and the Extreme Computing Research
Center (ECRC). Its contribution lies not in a new algorithm nor a new dataset,
but in demonstrating the routine use of the larger datasets becoming available to geospatial
statisticians, thanks to the implementation of state-of-the-art statistical algorithms on
High Performance Computing (HPC) hardware.

We have built a standalone software framework (ExaGeoStat/ExaGeoStatCPP) that can run on a variety
of hardware resources, including GPUs and massively distributed systems such as Shaheen-II,
KAUST's Cray XC40 supercomputer, HLRS HPE Apollo (Hawk), ORNL Summit (OLCF-4) supercomputer, and Riken Fugaku supercomputer,
to create a statistical model to predict environmental data (i.e., temperature, flow rates, soil moisture,
wind speed, air pollution, etc.) at spatial locations on which data
is missing, and to exploit large amounts of data to reduce the effect of individual measurement
errors. The best-known methods for such statistical processing have a cost that grows rapidly
in the size of the dataset, namely, in proportion to its cube or third power. Thus, increasing
the size of the dataset by a factor of ten drives up the cost of the computation by a factor of
a thousand while simultaneously driving up the memory requirements by a factor of a hundred.

For instance, according to this cubic growth in complexity, a computation that requires one
minute would require nearly 17 hours on a dataset just ten times larger. This creates a
computational strain on standard statistics software, for which contemporary data sizes
were not anticipated, and even if possible, it puts the computation beyond the interactive
attention span of the analyst. Parallelism (assigning thousands of processors to a single task) and Moore's Law allow
leading-edge computers to handle such "big data"
with ease, but the software bridge must be built. Furthermore, the software interface
must resemble the interactive one with which working statisticians are familiar.

To summarize, the combination of emerging computing capabilities and emerging datasets
promises significant advances in statistical analyses of environmental and many other
phenomena. Such cross-disciplinary advances are natural at KAUST, so this
relatively low-hanging fruit was ours to harvest earliest. Our roadmap now takes ExaGeoStat
a step further on the algorithmic side by integrating tile low-rank matrix approximation.
This low-rank matrix approximation permits the exploitation of the data sparsity of the operator with user-controlled
numerical accuracy. This further expands practical problem sizes for
statisticians with modest computational resources.


Current Version of ExaGeoStatCPP: 1.0.0
======================
Operations:

1. (Data Generation): Generating large geospatial synthetic datasets using  dense, Diagonal Super-Tile (DST) and Tile Low-Rank (TLR) approximation techniques.
2. (Data Modeling): Modeling large geospatial datasets on dense, Diagonal Super-Tile (DST) and Tile Low-Rank (TLR) approximation techniques through the Maximum likelihood Estimation (MLE) operation.
3. (Data Prediction): Predicting missing measurements on given locations using dense, Diagonal Super-Tile (DST), and Tile Low-Rank (TLR) approximation techniques.
4. (MLOE/MMOM): Computing the Mean Loss of Efficiency (MLOE), Mean Misspecification of the Mean Square Error (MMOM), and Root mean square MOM (RMOM) to describe the prediction performance over the whole observation region.
5. (Fisher Information Matrix (FIM)): Quantifying the information content that a variable x carries about a parameter $\theta$ within a Gaussian distribution.

Supported Covariance Functions:
======================

1. Univariate Matérn (Gaussian/Stationary)
2. Univariate Matérn with Nugget (Gaussian/Stationary)
3. Flexible Bivariate Matérn (Gaussian/Stationary)
4. Parsimonious Bivariate Matérn (Gaussian/Stationary)
5. Parsimonious trivariate Matérn (Gaussian/Stationary)
6. Univariate Space/Time Matérn (Gaussian/Stationary)
7. Bivariate Space/Time Matérn (Gaussian/Stationary)
8. Tukey g-and-h Univariate Matérn (non-Gaussian/Stationary)
9. Tukey g-and-h Univariate Power Exponential (non-Gaussian/Stationary)

Programming models:

1. MPI
2. Task-based programming models

External libraries:

1. NLOPT [https://nlopt.readthedocs.io/en/latest/](https://nlopt.readthedocs.io/en/latest/)
2. GSL [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)
3. HWLOC [https://www.open-mpi.org/projects/hwloc/](https://www.open-mpi.org/projects/hwloc/)
4. StarPU dynamic runtime system  [https://starpu.gitlabpages.inria.fr/](https://starpu.gitlabpages.inria.fr/)
5. HCORE [https://github.com/ecrc/hcore](https://github.com/ecrc/hcore)
6. HiCMA [https://github.com/ecrc/hicma](https://github.com/ecrc/hicma)
7. Stars-H [https://github.com/ecrc/stars-h](https://github.com/ecrc/stars-h)
8. Chameleon [https://gitlab.inria.fr/solverstack/chameleon](https://gitlab.inria.fr/solverstack/chameleon)

Project Hierarchy
--------------------

* **```cmake```** A directory contains essential CMake modules that facilitate the importation and location of required dependencies.
* **```docs```** A directory contains all the necessary documents.
* **```examples```** A directory contains a comprehensive collection of demo code that illustrates the framework's application and demonstrates its features and capabilities.
* **```inst```** A directory contains all the system's header files, mirroring the structure of the src directory.
* **```prerequisites```** A directory contains all the necessary prerequisites for the project, as well as default scripts for their installation.
* **```src```** A directory contains all the source files.
* **```tests```** A directory contains all the test files and follows the same structure as the src folder.
* **```clean_build.sh```** A script is designed to compile the software tests once all dependencies are installed, and it is set to build everything by default.
* **```CMakeLists.txt```** The top-level CMake file to configure the build system.
* **```config.sh```** Script used to generate the building system inside a 'bin' directory.

Installation
============

Installation requires at least **CMake of version 3.2**. to build ExaGeoStatCPP,
please follow these instructions:

1. Get from git repository

       git clone git@github.com:ecrc/exageostatcpp

   or

       git clone https://github.com/ecrc/exageostatcpp

2. Go into the ExaGeoStatCPP folder

       cd exageostatcpp

3. Run help of config.sh to know the needed arguments to run with your specific options.

       ./config.sh --h
   or check user manual.

4. Run help of clean_build.sh to know the needed arguments to run with your specific options.

       ./clean_build.sh -h

10. Export the installation paths of the dependencies, e.g.,

        export PKG_CONFIG_PATH=$PWD/installdir/_deps/DEPENDENCY_NAME/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now, you can use the pkg-config executable to collect compiler and linker flags for
EXAGEOSTATCPP.

Using ExaGeoStatCPP
============
Please refer to **```USER_MANUAL```** for detailed instructions.
Please take a look at the end-to-end examples as a reference for using all the operations.

Contribute
=======

[Contribution Guidelines](CONTRIBUTING.md)

References
==========

1. Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "ExaGeoStat: A high performance unified
   software for geostatistics on manycore systems." IEEE Transactions on Parallel and Distributed Systems 29, no. 12 (
   2018): 2771-2784.

2. Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Parallel approximation of the maximum
   likelihood estimation for the prediction of large-scale geostatistics simulations." In 2018 IEEE International Conference
   on Cluster Computing (CLUSTER), pp. 98-108. IEEE, 2018.

3. Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Geostatistical modeling and prediction
   using mixed precision tile Cholesky factorization." In 2019 IEEE 26th international conference on high performance
   computing, data, and analytics (HiPC), pp. 152-162. IEEE, 2019.

4. Mary Lai O. Salvana, Sameh Abdulah, Huang Huang, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "High
   performance multivariate geospatial statistics on manycore systems." IEEE Transactions on Parallel and Distributed
   Systems 32, no. 11 (2021): 2719-2733.

5. Mary Lai O. Salvaña, Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Parallel Space-Time
   Likelihood Optimization for Air Pollution Prediction on Large-Scale Systems." In the Proceedings of the Platform for
   Advanced Scientific Computing Conference (PASC'22). Association for Computing Machinery, New York, NY, USA, Article
   17, 1–11. ACM, 2022.

6. Sameh Abdulah, Qinglei Cao, Yu Pei, George Bosilca, Jack Dongarra, Marc G. Genton, David E. Keyes, Hatem Ltaief, and
   Ying Sun. "Accelerating geostatistical modeling and prediction with mixed-precision computations: A high-productivity
   approach with PaRSEC." IEEE Transactions on Parallel and Distributed Systems 33, no. 4 (2021): 964-976.

7. Sagnik Mondal, Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Parallel Approximations
   of the Tukey g-and-h Likelihoods and Predictions for Non-Gaussian Geostatistics." 2022 IEEE International Parallel
   and Distributed Processing Symposium (IPDPS), Lyon, France, 2022, pp. 379-389. IEEE, 2022.

8. Qinglei Cao, Sameh Abdulah, Rabab Alomairy, Yu Pei, Pratik Nag, George Bosilca, Jack Dongarra et al. "Reshaping
   geostatistical modeling and prediction for extreme-scale environmental applications." In 2022 SC22: International
   Conference for High-Performance Computing, Networking, Storage and Analysis (SC), pp. 13-24. IEEE Computer Society, 2022.
   (ACM GORDON BELL PRIZE Finalist).
         
9. Sagnik Mondal, Sameh Abdulah, Hatem Ltaief, Ying Sun, Marc G. Genton, and David E. Keyes. "Tile low-rank approximations
    of non-Gaussian space and space-time Tukey g-and-h random field likelihoods and predictions on large-scale systems." 
    Journal of Parallel and Distributed Computing 180 (2023): 104715.
    
10. Qinglei Cao, Sameh Abdulah, Hatem Ltaief, Marc G. Genton, David E. Keyes, and George Bosilca. "Reducing Data Motion
    and Energy Consumption of Geospatial Modeling Applications Using Automated Precision Conversion." In 2023 IEEE International Conference
   on Cluster Computing (CLUSTER), IEEE, 2023.
