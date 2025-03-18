
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2023-02-14
**/

#include <fstream>
#include <nlopt.hpp>
#include <mpi.h>

#include <data-transformer/DataTransformer.hpp>
#include <data-loader/concrete/ParsecLoader.hpp>
#include <runtime-solver/RuntimeSolverFactory.hpp>
extern "C"{
    #include <netcdf.h>
    #include <pnetcdf.h>

    #include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

using namespace std;
using namespace nlopt;

using namespace exageostat::dataLoader::parsec;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;
using namespace exageostat::transformers;

template<typename T>
ParsecLoader<T> *ParsecLoader<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new ParsecLoader<T>();
    }
    return mpInstance;
}

template<typename T>
unique_ptr<ExaGeoStatData<T>>
ParsecLoader<T>::LoadData(Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) {

    //create data object
    auto data = make_unique<ExaGeoStatData<T>>(aConfigurations.GetProblemSize(), aConfigurations.GetDimension());
    
    // Stage 0; File Generation
    if(aConfigurations.GetStageZero()){
        StageZeroGeneration(aConfigurations, data, aKernel);
    }
//    // Initiate Descriptors
//    InitFilesDescriptors(aConfigurations, data);
//    // Forward SHT
//    DataTransformer<T>::ForwardSphericalHarmonicsTransform(aConfigurations.GetDenseTileSize(), data);
//    // Forward SHT Reshape
//    DataTransformer<T>::ForwardReshape(aConfigurations, data);
//    // Generate matrix
//    CompressMatrixHelper(aConfigurations, data);

    return data;
}

template<typename T>
int ParsecLoader<T>::OpenNetCDFFile(const string &aFilename) {

    VERBOSE("Open the given " << aFilename << " NetCDF file.");

//    int netcdf_id, return_val;
//    if ((return_val = ncmpi_open(MPI_COMM_WORLD, aFilename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &netcdf_id)))
//        printf("Error: %s\n", nc_strerror(return_val));
//
//    VERBOSE("Done");
//    return netcdf_id;

    int ncid = 0, return_val = 0;

    //    if (return_val = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid))
    if ((return_val = nc_open(aFilename.c_str(), NC_NOWRITE, &ncid)))
	    printf("Error: %s --- %s\n", nc_strerror(return_val), aFilename.c_str());
    VERBOSE("Done \n");
    return ncid;
}

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
	VERBOSE("Done \n");
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
	VERBOSE("Done\n");
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

template<typename T>
void ParsecLoader<T>::StageZeroGeneration(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData, kernels::Kernel<T> &aKernel) {

    auto forcing_file = aConfigurations.GetForcingDataPath();
    auto MAX_ITERATIONS = 30;
    std::ifstream file(forcing_file);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open observations file: " + forcing_file);
    }

    std::vector<double> observations;

    std::string line;
    while (std::getline(file, line)) {
        try {
            observations.push_back(std::stod(line));
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing observation value: " + line);
        }
    }

    fprintf(stderr, "Number of params to optimize: %d\n", aKernel.GetParametersNumbers());

    // Print the first two values for debugging.
    std::cerr << "forcing: " << observations[0] << " " << observations[1] << std::endl;

    // Assigning forcing to a hypothetical data structure.
    std::vector<double> forcing;
    forcing = observations;

    std::cerr << "data.forcing: " << forcing[0] << " " << forcing[1] << std::endl;




    // Optimizater initialization
	//NLOPT_LN_BOBYQA
    opt optimizing_function(nlopt::LN_BOBYQA, aKernel.GetParametersNumbers());
    // Initialize problem's bound.
    if(aConfigurations.GetLowerBounds().empty()){
        throw runtime_error("Please set the lower bounds in order to use this Kernel.");
    }
    if(aConfigurations.GetUpperBounds().empty()){
        throw runtime_error("Please set the upper bounds in order to use this Kernel.");
    }
    optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
    optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
    optimizing_function.set_ftol_abs(aConfigurations.GetTolerance());
    optimizing_function.set_maxeval(MAX_ITERATIONS);

	fprintf(stderr, "data.forcing: %f %f\n", forcing[0], forcing[1]);

    // TODO: ALL OF THESE WILL BE REMOVED TILL WE REVIEW IT WITH SAMEH
    aConfigurations.SetProblemSize(365 * 24* (2023-1988));
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


//    string netcdf_data_path = aConfigurations.GetNetCDFDataPath();
//    int netcdf_id;
//    char path[150] = "/home/mahmoudelkarargy/Documents/KAUST/stage0_files/data_1940.nc";
//    netcdf_id = openFileNCmpi(path);
//    printf("netcdf_id: %d\n", netcdf_id);
//    int return_val, longitude_varid, latitude_varid, time_varid, t2m_varid;
//    size_t longitude_len, latitude_len, time_len;
//
//    if ((return_val = nc_inq_dimid(netcdf_id, "longitude", &longitude_varid)))
//		printf("Error: %s\n", nc_strerror(return_val));
//	if ((return_val = nc_inq_dimlen(netcdf_id, longitude_varid, &longitude_len)))
//		printf("Error: %s\n", nc_strerror(return_val));
//
//	if ((return_val = nc_inq_dimid(netcdf_id, "latitude", &latitude_varid)))
//		printf("Error: %s\n", nc_strerror(return_val));
//	if ((return_val = nc_inq_dimlen(netcdf_id, latitude_varid, &latitude_len)))
//		printf("Error: %s\n", nc_strerror(return_val));
//
//	if ((return_val = nc_inq_dimid(netcdf_id, "time", &time_varid)))
//		printf("Error: %s\n", nc_strerror(return_val));
//	if ((return_val = nc_inq_dimlen(netcdf_id, time_varid, &time_len)))
//		printf("Error: %s\n", nc_strerror(return_val));
//
//	int length = longitude_len * latitude_len * time_len;
//
//	fprintf(stderr, "%d, %d, %d\n", longitude_varid, latitude_varid, time_varid);
//	fprintf(stderr, "%d, %d, %d, %d\n", longitude_len, latitude_len, time_len, length);
//
//	if ((return_val = nc_inq_varid(netcdf_id, "t2m", &t2m_varid)))
//		printf("Error: %s\n", nc_strerror(return_val));

    int ncid;
    int retval;
    size_t  t2m_len;
	int  lon_varid, lat_varid, time_varid,t2m_varid;
	MPI_Offset lat_len, lon_len, time_len;
int 	len;
    char numStr[20];
	char path[150] = "/home/mahmoudelkarargy/Documents/KAUST/stage0_files/data_1940.nc";

	ncid 	= openFileNCmpi(path);
    printf("netcdf_id: %d\n", ncid);
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
    int v= 365;//*no_years;
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
    size_t lts = 0;
    int min_loc=lts;

    printf("lts: %d\n", lts);
    printf("maxrank: %d\n", 0);
	int max_loc= 0;
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
        printf("before reading!! %s\n", filePath.c_str());
//		netcdf_id    = OpenNetCDFFile(netcdf_data_path);
        ncid    = openFileNCmpi(path);
        printf("id: %d\n", ncid);
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
//		size_t index[] = {(time_len/nprocs)*rank, lts, 0};
//		size_t count[] = {time_len/nprocs, 1, no_locs};
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
			//fprintf(stderr, ">>> %d %d\n", num_locs, no_locs);
			double sum_temp = 0;
			int sum_tmp=0;
			//ll=0;
			//to be removed
			// 8784
			int r=0;
			if(isLeapYear(y))
			{
				r=0;
				//fprintf(stderr, ">>> %d %d\n", num_locs, no_locs);
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
//    MORSE_dmle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N,  0, 0);
                            //            fprintf(stderr, "hi2\n");
                            printf("1\n");
	int count = 0;
    auto runtime_solver = runtimesolver::RuntimeSolverFactory<T>::CreateRuntimeSolver();
                            printf("2\n");
//    T result = runtime_solver->ModelingOperations(aData, aConfigurations, nullptr, aKernel);
                            printf("3\n");


//    nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
    printf("how exay dh size lts: %d\n", lts);
    printf("how exay dh size noLocs: %d\n", no_locs);
    for (int u=lts;u<lts+1; u++)
		for (int l=0;l<no_locs; l++)//no_locs
		{

			fprintf(stderr, "u=%d, l=%d\n", u, l);
			//if(l ==data.hicma_maxrank)
			if(true)
			{
//				data.iter_count=0;
                #if defined(CHAMELEON_USE_MPI)
				if(MORSE_My_Mpi_Rank() == 0)
				{
#endif
//                    MORSE_MLE_dzcpy(&data, &t2m_hourly_per_year[count][0]);
                    auto *CHAM_descZ = aData->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).parsec_desc;
                    VERBOSE("Copy Z from vector to decriptor.\n");
                    runtime::RuntimeFunctions<T>::CopyDescriptorZ(*aData->GetDescriptorData(), CHAM_descZ, (T *) &t2m_hourly_per_year[count][0]);
                    VERBOSE("Done Z copying step.\n");
	                VERBOSE("************************************************************\n");
//					START_TIMING(data.total_exec_time);
//					nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
    auto runtime_solver = runtimesolver::RuntimeSolverFactory<T>::CreateRuntimeSolver();
                            printf("2\n");
    T result = runtime_solver->ModelingOperations(aData, aConfigurations, nullptr, aKernel);
                            printf("3\n");
//					starting_theta[0] = 0.9;
//					nlopt_optimize(opt, starting_theta, &opt_f);
//					STOP_TIMING(data.total_exec_time);
					//	starting_theta[0] = 0.272049;
//					fprintf(stderr, "final theta: %f\n", starting_theta[0]);

//					mean_trend(starting_theta, &data, l, u);
                    #if defined(CHAMELEON_USE_MPI)
				}
#endif
			}
			count++;
            }

    fprintf(stderr, "Done\n");
}

template<typename T>
void ParsecLoader<T>::InitFilesDescriptors(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData) {
    
    SYNC_TIME_START();
    int L = aConfigurations.GetDenseTileSize();
    int MB;
    int NB;
    int t = aConfigurations.GetTimeSlot();
    int P = aConfigurations.GetPGrid();
    int nodes = aConfigurations.GetCoresNumber();
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int verbose = Configurations::GetVerbosity() == DETAILED_MODE? 1: 0;
    int gpus = aConfigurations.GetGPUsNumbers();
    int tile_size = aConfigurations.GetDenseTileSize();
    string files_directory_path = aConfigurations.GetDataPath();
    int path_length = files_directory_path.length();
    char filename[path_length + 50];
    char directory_path[path_length];
    sprintf(directory_path, "%s", files_directory_path.c_str());

    MB = L + 1;
    NB = L * 2;
	VERBOSE_PRINT(rank, verbose, ("Reading f_data\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA);
    parsec_matrix_block_cyclic_t *pF_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA).parsec_desc;
    ReadCSVToComplexTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_data_desc, MB, NB, nodes, t, directory_path, rank, verbose, gpus);

    MB = 2*L-1;
    NB = L+1;
	VERBOSE_PRINT(rank, verbose, ("Reading Et1\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1);
    parsec_matrix_block_cyclic_t *pEt1_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Et1.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt1_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Et2\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2);
    parsec_matrix_block_cyclic_t *pEt2_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Et2.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt2_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Ep\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP);
    parsec_matrix_block_cyclic_t *pEp_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Ep.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEp_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = (L*L+L)/2;
    NB = L;
	VERBOSE_PRINT(rank, verbose, ("Reading Slmn\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN);
    parsec_matrix_block_cyclic_t *pSlum_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Slmn.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSlum_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Ie\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE);
    parsec_matrix_block_cyclic_t *pIe_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Ie.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIe_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Io\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO);
    parsec_matrix_block_cyclic_t *pIo_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Io.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIo_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L-1;
    NB = L+1;
	VERBOSE_PRINT(rank, verbose, ("Reading P\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P);
    parsec_matrix_block_cyclic_t *pP_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_P.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pP_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading D\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D);
    parsec_matrix_block_cyclic_t *pD_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_D.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pD_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = L;
	VERBOSE_PRINT(rank, verbose, ("Reading flm\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM);
    parsec_matrix_block_cyclic_t *pFlm_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM).parsec_desc;
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlm_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

	VERBOSE_PRINT(rank, verbose, ("Reading flmERA\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMERA);
    parsec_matrix_block_cyclic_t *pFlmera_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMERA).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_flmERA.csv");
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlmera_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    // Backward
    if(aConfigurations.GetEnableInverse()){

        MB = L+1;
        NB = (L*L+L)/2;
        VERBOSE_PRINT(rank, verbose, ("Reading Zlm\n"));
        aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM);
        parsec_matrix_block_cyclic_t *PZlm_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM).parsec_desc;
        sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Zlm.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), PZlm_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

        MB = 2*L-1;
        NB = 2*L;
        VERBOSE_PRINT(rank, verbose, ("Reading SC\n"));
        aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC);
        parsec_matrix_block_cyclic_t *pSc_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC).parsec_desc;
        sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_SC.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSc_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

        MB = L+1;
        NB = 2*L;
        VERBOSE_PRINT(rank, verbose, ("f_spatial\n"));
        aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL);
        parsec_matrix_block_cyclic_t *pF_spatial_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL).parsec_desc;
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_spatial_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);
    }

    // Init and allocate memory for desc_flmT
    MB = L * L;
    NB = t;
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT);
    parsec_matrix_block_cyclic_t *pFlmt_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT).parsec_desc;
    parsec_matrix_block_cyclic_init(pFlmt_data_desc, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank, L*((MB/nodes%L) ? MB/nodes/L+1 : MB/nodes/L),
                                    NB, MB, NB, 0, 0, MB, NB, nodes, 1, 1, 1, 0, 0);

    pFlmt_data_desc->mat = parsec_data_allocate((size_t)pFlmt_data_desc->super.nb_local_tiles *
                                   (size_t)pFlmt_data_desc->super.bsiz *
                                   (size_t)parsec_datadist_getsizeoftype(pFlmt_data_desc->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&desc_flmT, "desc_flmT");

    // Init and allocate memory for pA_data_desc
    MB = L * L;
    NB = t;
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_A);
    parsec_matrix_block_cyclic_t *pA_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_A).parsec_desc;
    parsec_matrix_block_cyclic_init(pA_data_desc, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank, L, L, MB, NB, 0, 0,
                                    pFlmt_data_desc->super.mb, pFlmt_data_desc->super.nb, P, nodes/P, 1, 1, 0, 0);
    pA_data_desc->mat = parsec_data_allocate((size_t)pA_data_desc->super.nb_local_tiles *
                                   (size_t)pA_data_desc->super.bsiz *
                                   (size_t)parsec_datadist_getsizeoftype(pA_data_desc->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&pA_data_desc, "desc_A");

    if(aConfigurations.GetEnableInverse()){
        int ts_test_M = 2000;
        int ts_test_N = 1;
        // Allocate memory
        double *pFileContent = (double *)malloc(ts_test_M * ts_test_N * sizeof(double));
        sprintf(filename, "%s/%s", files_directory_path.c_str(), "ts_test.csv");
        ReadCSVFileHelper(filename, pFileContent, ts_test_M, ts_test_N);
    }

    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Load Data\n"));
}
template<typename T>
int ParsecLoader<T>::ReadCSVFileHelper(const char* apFilename, double *apFileContent, int aM, int aN) {

    FILE *pFile = fopen(apFilename, "r");
    if (!pFile) {
        printf("File opening failed: %s", apFilename);
        return -1;
    }

    int status = 0;
    for (int i = 0; i < aM; i++) {
        for (int j = 0; j < aN; j++) {
            status = fscanf(pFile, "%lf,", &apFileContent[i * aN + j]);
            if (status != 1) {
                fprintf(stderr, "Error reading file at row %d, column %d\n", i, j);
                fclose(pFile);
                return 1;
            }
        }
    }
    fclose(pFile);

    return 0;
}

template<typename T>
void ParsecLoader<T>::CompressMatrixHelper(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData) {

    int max_rank = aConfigurations.GetMaxRank();
    int n = aConfigurations.GetProblemSize();
    int adaptive_decision = aConfigurations.GetAdaptiveDecision();
    int tol = aConfigurations.GetTolerance();
    int send_full_tile = 0;
    int auto_band = 0;
    int gpus = aConfigurations.GetGPUsNumbers();
    double upper_lower = EXAGEOSTAT_LOWER;
    int L = aConfigurations.GetDenseTileSize();
    int N = aConfigurations.GetProblemSize();
    int NT = (N % L == 0) ? (N/L) : (N/L + 1);

    SYNC_TIME_START();
    MatrixCompress((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), &ExaGeoStatHardware::GetHicmaParams()->norm_global, upper_lower, L, NT, max_rank, n,
                   adaptive_decision, tol, send_full_tile, auto_band, gpus, ExaGeoStatHardware::GetHicmaData(), ExaGeoStatHardware::GetParamsKernel());

    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Matrix genneration Matrix norm: norm_global= %le\n", ExaGeoStatHardware::GetHicmaParams()->norm_global));
}

template<typename T>
void ParsecLoader<T>::ReadData(Configurations &aConfigurations, vector<T> &aMeasurementsMatrix, vector<T> &aXLocations,
                               vector<T> &aYLocations, vector<T> &aZLocations, const int &aP) {
}

template<typename T>
void ParsecLoader<T>::WriteData(const T &aMatrixPointer, const int &aProblemSize, const int &aP, string &aLoggerPath, Locations<T> &aLocations) {

}

template<typename T>
void ParsecLoader<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> ParsecLoader<T> *ParsecLoader<T>::mpInstance = nullptr;
