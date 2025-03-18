
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarpuRuntimeSolver.cpp
 * @brief This file contains the implementation of StarpuRuntimeSolver class.
 * @details StarpuRuntimeSolver is a concrete implementation of the RuntimeSolversMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-11-04
**/

#include <runtime-solver/concrete/StarpuRuntimeSolver.hpp>
#include <data-units/ModelingDataHolders.hpp>
#include <utilities/Logger.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

#include <fstream>
#include <nlopt.hpp>
#include <mpi.h>

extern "C"{
    #include <netcdf.h>
    #include <pnetcdf.h>
}


using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::runtimesolver;
using namespace exageostat::dataunits;
using namespace nlopt;
using namespace std;


int openFileNCmpi(char *filename)
	//!Open the NetCDF file.
	/*!
	 * returns  ncid: netCDF file identifier.
	 * @param[in] data: MLE_data struct with different MLE inputs.
	 * @param[in] filename: netCDF location file.
	 */
{
	VERBOSE("Open the given NetCDF file... ");
    printf("%s\n",filename);
	int ncid, retval;

	if (retval = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid))
		//	    if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
		printf("Error: %s\n", nc_strerror(retval));
	VERBOSE("Done");
	return ncid;

}

void closeFileNCmpi(int ncid)
	//!Close the NetCDF file.
	/*!
	 * @param[in] ncid: netCDF file identifier.
	 * @param[in] data: MLE_data struct with different MLE inputs.
	 */
{
	int retval;

	VERBOSE("Close the NetCDF file... ");
	if ((retval = ncmpi_close(ncid)))
		//              if ((retval = nc_close(ncid)))
		printf("Error: %s\n", nc_strerror(retval));
	VERBOSE("Done");
}

int isLeapYear(int year) {
    if (year % 400 == 0) {
        return 1;
    } else if (year % 100 == 0) {
        return 0;
    } else if (year % 4 == 0) {
        return 1;
    } else {
        return 0;
    }
}

double * readObsFile(char *obsfile, int n)
//! Read real observation file
/* Returns observations struct.
 *  @param[in] obsfile: observation file  path.
 *  @param[in]  n : number of spatial locations (observations).
 *  @param[out] streamdata: observations vector.
 */
{

    FILE * fp;
    char * line     = NULL;
    size_t len      = 0;
    ssize_t read;
    int count       = 0;

    double *z_vec = (double *) malloc(n * sizeof(double));

    fp = fopen(obsfile, "r");
    if (fp == NULL)
    {
        printf("readObsFile:cannot open observations file: %s\n", obsfile);
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &len, fp)) != -1)
        z_vec[count++]=atof(line);

    fclose(fp);
    free(line);

    return z_vec;
}

template<typename T>
T StarpuRuntimeSolver<T>::ModelingOperations(std::unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfigurations,
                                              T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

    double opt_f;
    int M=10;
    int t=365*24;  //Fixed
    int no_years= 751;
    int num_params = 0;
    size_t N, dts, lts = 0, log;

    auto forcing_file = aConfigurations.GetForcingDataPath();
    auto MAX_ITERATIONS = 30;
    std::ifstream file(forcing_file);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open observations file: " + forcing_file);
    }

//    std::vector<double> observations;
//
//    std::string line;
//    while (std::getline(file, line)) {
//        try {
//            observations.push_back(std::stod(line));
//        } catch (const std::exception& e) {
//            throw std::runtime_error("Error parsing observation value: " + line);
//        }
//    }
    fprintf(stderr, "Number of params to optimize: %d\n", aKernel.GetParametersNumbers());

    double * forcing = (double *) malloc(751 * sizeof(double));
    forcing = readObsFile("/home/mahmoudelkarargy/Documents/KAUST/stage0_files/forcing_new.csv", 751);
    fprintf(stderr, "forcing: %0.16f %e\n", forcing[0], forcing[1]);
    apMeasurementsMatrix = (T*) forcing;
    fprintf(stderr, "data.forcing: %f %f\n", apMeasurementsMatrix[0], apMeasurementsMatrix[1]);
    // Optimizater initialization

    int v= 365;//*no_years;
    N = 365 * 24* (2023-1988);// N =time_len;

	//NLOPT_LN_BOBYQA
    opt optimizing_function(nlopt::LN_BOBYQA, aKernel.GetParametersNumbers());
    // Initialize problem's bound.
    if(aConfigurations.GetLowerBounds().empty()){
        throw std::runtime_error("Please set the lower bounds in order to use this Kernel.");
    }
    if(aConfigurations.GetUpperBounds().empty()){
        throw std::runtime_error("Please set the upper bounds in order to use this Kernel.");
    }
    optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
    optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
    optimizing_function.set_ftol_abs(aConfigurations.GetTolerance());
    optimizing_function.set_maxeval(MAX_ITERATIONS);

	fprintf(stderr, "data.forcing: %f %f\n", apMeasurementsMatrix[0], apMeasurementsMatrix[1]);

    // TODO: ALL OF THESE WILL BE REMOVED TILL WE REVIEW IT WITH SAMEH
    aConfigurations.SetProblemSize(N);
    fprintf(stderr,"********************SUMMARY**********************\n");
    fprintf(stderr,"#Synthetic Dataset\n");
    fprintf(stderr,"Number of Locations: %d\n", aConfigurations.GetProblemSize());
    fprintf(stderr,"#Threads per node: %d\n", aConfigurations.GetCoresNumber());
    fprintf(stderr,"#GPUs: %d\n", aConfigurations.GetGPUsNumbers());
    if(aConfigurations.GetPrecision() == common::DOUBLE)
        fprintf(stderr,"#Double Precision!\n");
    if (aConfigurations.GetComputation() == common::EXACT_DENSE)
        fprintf(stderr,"#Dense Tile Size: %d\n", aConfigurations.GetDenseTileSize());
    fprintf(stderr,"p=%d, q=%d\n", aConfigurations.GetPGrid(), aConfigurations.GetQGrid());
    fprintf(stderr,"***************************************************\n");

    int ncid;
    int retval;
    size_t  t2m_len;
	int  lon_varid, lat_varid, time_varid,t2m_varid;
	MPI_Offset lat_len, lon_len, time_len;
int 	len;
    char numStr[20];
	char path[150] = "/home/mahmoudelkarargy/Documents/KAUST/stage0_files/data_1940.nc";

	ncid 	= openFileNCmpi(path);
	fprintf(stderr, "%s\n", path);

	if ((retval = ncmpi_inq_dimid(ncid,"longitude",&lon_varid)))
		printf("Error: %s\n", nc_strerror(retval));
	if ((retval = ncmpi_inq_dimlen(ncid, lon_varid, &lon_len)))
		printf("Error: %s\n", nc_strerror(retval));
	//**********
	if ((retval = ncmpi_inq_dimid(ncid,"latitude",&lat_varid)))
		printf("Error: %s\n", nc_strerror(retval));
	if ((retval = ncmpi_inq_dimlen(ncid, lat_varid, &lat_len)))
		printf("Error: %s\n", nc_strerror(retval));
	//*********
	if ((retval = ncmpi_inq_dimid(ncid,"time",&time_varid)))
		printf("Error: %s\n", nc_strerror(retval));
	if ((retval = ncmpi_inq_dimlen(ncid, time_varid, &time_len)))
		printf("Error: %s\n", nc_strerror(retval));


	len = lon_len * lat_len * time_len;
	fprintf(stderr, "%d, %d, %d\n", lon_varid, lat_varid, time_varid);
	fprintf(stderr, "%d, %d, %d, %d\n", lon_len, lat_len, time_len, len);
	if ((retval = ncmpi_inq_varid(ncid, "t2m", &t2m_varid)))
		printf("Error: %s\n", nc_strerror(retval));



    size_t no_locs=1440;
	double *t2m_hourly_per_year[no_locs];
	int t2m_hourly_per_year_count[no_locs];
	for(int k=0;k<no_locs;k++)
	{
		t2m_hourly_per_year[k] = (double *) malloc(v * 24* (2023-1988)  *   sizeof(double));
		if(t2m_hourly_per_year[k] == NULL)
			fprintf(stderr, "no enough memory... t2m_hourly_per_year");
		t2m_hourly_per_year_count[k]=0;
	}

    int rank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	// TODO: Why lts & maxrank?
    int min_loc=lts;
    int maxrank = 0;
	int max_loc= maxrank;
	int r=0;
	double sum_obs=0;
	double  scaling_var, offset_var;

//    char numStr[20];
    int* t2m;
    int *t2m_local;
    for(int y=1940;y<1941;y++)
	{
		// Convert integer to string
		string basePath = "/home/mahmoudelkarargy/Documents/KAUST/stage0_files/data_";
        string filePath = "/home/mahmoudelkarargy/Documents/KAUST/stage0_files/data_1940.nc";
        ncid    = openFileNCmpi(path);
		if(isLeapYear(y))
		{
			t2m = (int *) malloc((v*24+24) * no_locs * sizeof(int));
			t2m_local = (int *) malloc((v*24+24) * no_locs * sizeof(int)/nprocs);
			if(t2m == NULL|| t2m_local == NULL)
				fprintf(stderr, "no enough memory... t2m or t2m_local");
			time_len=v*24+24;
		}
		else
		{	t2m = (int *) malloc(v * 24 * no_locs* sizeof(int));
			t2m_local = (int *) malloc((v*24) * no_locs* sizeof(int)/nprocs);
			if(t2m == NULL|| t2m_local == NULL)
				fprintf(stderr, "no enough memory... t2m or t2m_local");
			time_len=v*24;
		}

		double x_read=0;
		START_TIMING(x_read);
        MPI_Offset index[] = {(MPI_Offset)(time_len / nprocs) * rank, (MPI_Offset)lts, 0};
        MPI_Offset count[] = {(MPI_Offset)(time_len / nprocs), 1, (MPI_Offset)no_locs};

		VERBOSE("Start reading......\n");
		if (retval = ncmpi_get_vara_int_all(ncid, t2m_varid, index, count, t2m_local))
			printf("ncmpi_get_vara_int_all: Error: %s -- %s \n", nc_strerror(retval), filePath.c_str());
		if (retval =ncmpi_get_att_double(ncid, t2m_varid, "scale_factor", &scaling_var))
			printf("ncmpi_get_att_double: Error: %s\n", nc_strerror(retval));
		if (retval =ncmpi_get_att_double(ncid, t2m_varid, "add_offset", &offset_var))
			printf("ncmpi_get_att_double: Error: %s\n", nc_strerror(retval));

                VERBOSE("End reading.\n");

		MPI_Allgather(t2m_local, (time_len/nprocs)*no_locs, MPI_INT, t2m, (time_len/nprocs)*no_locs, MPI_INT, MPI_COMM_WORLD);
                VERBOSE("End gathering.\n");

                STOP_TIMING(x_read);
#if defined(CHAMELEON_USE_MPI)
		if(MORSE_My_Mpi_Rank() == 0)
		{
#endif
			fprintf(stderr, "================================ \n");
			fprintf(stderr, "%s\n", filePath.c_str());
			fprintf(stderr, "time to read the NETCDF file: %f secs\n", x_read);
#if defined(CHAMELEON_USE_MPI)
		}
#endif

        int num_locs=0;
		for (int lu=0;lu<no_locs; lu++)
		{
//			fprintf(stderr, ">>> %d %d\n", num_locs, no_locs);
			double sum_temp = 0;
			int sum_tmp=0;
			int r=0;
			if(isLeapYear(y))
			{
				r=0;
//				fprintf(stderr, ">>> %d %d\n", num_locs, no_locs);
				for(int k=lu;k<(v*24+24)*no_locs;k+=no_locs)
				{

					if(r>= 1416 && r<=1439);
					else
					{
							t2m_hourly_per_year[lu][t2m_hourly_per_year_count[lu]++]=((double)t2m[k]*scaling_var) + offset_var - 273.15;// = sum_temp/24.0;
if((t2m_hourly_per_year[lu][t2m_hourly_per_year_count[lu]-1]) >10)
fprintf(stderr, " large value : %lf \n", (t2m_hourly_per_year[lu][t2m_hourly_per_year_count[lu]-1]));
					}
					r++;
				}
			}
			else
			{
				r=0;
				for(int k=lu;k<v*24*no_locs;k+=no_locs)
				{
						t2m_hourly_per_year[lu][t2m_hourly_per_year_count[lu]++] = ((double)t2m[k]*scaling_var) + offset_var - 273.15; //sum_temp/24.0;
if((t2m_hourly_per_year[lu][t2m_hourly_per_year_count[lu]-1]) >10)
fprintf(stderr, " large value : %lf \n", (t2m_hourly_per_year[lu][t2m_hourly_per_year_count[lu]-1]));
					r++;
				}
				num_locs++;
			}
		}
		//*****************************

		closeFileNCmpi(ncid);//	free (path2);
	}
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(aConfigurations.GetComputation());
    linear_algebra_solver->InitiateDescriptors(aConfigurations, *aData->GetDescriptorData(), aKernel.GetVariablesNumber(), nullptr);
	int count = 0;
    auto modeling_data = new mModelingData(aData, aConfigurations, *apMeasurementsMatrix, aKernel);
    optimizing_function.set_max_objective(DataModelingAPI, (void *) modeling_data);

    int l, u;
    for (u=lts;u<lts+1; u++)
		for (l=0;l<no_locs; l++)//no_locs
		{

			fprintf(stderr, "u=%d, l=%d\n", u, l);
			//if(l ==data.hicma_maxrank)
			if(true)
			{
                aData->SetMleIterations(0);
#if defined(CHAMELEON_USE_MPI)
				if(MORSE_My_Mpi_Rank() == 0)
				{
#endif
                    auto *CHAM_descZ = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
                    VERBOSE("Copy Z from vector to decriptor.")
                    runtime::RuntimeFunctions<T>::CopyDescriptorZ(*aData->GetDescriptorData(), CHAM_descZ, (T *) &t2m_hourly_per_year[count][0]);
                    VERBOSE("Done Z copying step.")
	                VERBOSE("************************************************************\n");
//					START_TIMING(data.total_exec_time);
                    optimizing_function.set_max_objective(DataModelingAPI, (void *) modeling_data);
                    aConfigurations.GetStartingTheta()[0] = 0.9;
                    optimizing_function.optimize(aConfigurations.GetStartingTheta(), opt_f);
                    fprintf(stderr, "final theta: %f\n",  aConfigurations.GetStartingTheta()[0]);
					linear_algebra_solver->mean_trend(aData, aConfigurations, aKernel, forcing, aConfigurations.GetStartingTheta().data(), l, u);
                    #if defined(CHAMELEON_USE_MPI)
				}
#endif
			}
			count++;
            }

    fprintf(stderr, "Done");
    return 0;
}

template<typename T>
double StarpuRuntimeSolver<T>::DataModelingAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo) {

    auto config = ((mModelingData<T> *) apInfo)->mpConfiguration;
    auto data = ((mModelingData<T> *) apInfo)->mpData;
    auto measurements = ((mModelingData<T> *) apInfo)->mpMeasurementsMatrix;
    auto kernel = ((mModelingData<T> *) apInfo)->mpKernel;

    // We do Date Modeling with any computation.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(EXACT_DENSE);
    return linear_algebra_solver->ExaGeoStatMLETile(*data, *config, aTheta.data(), measurements, *kernel);
}

//
//-DCMAKE_CXX_FLAGS_RELEASE="-Wl,--no-as-needed -w -fpic"
//-DCMAKE_CXX_FLAGS_DEBUG="-Wl,--no-as-needed -w -fpic"