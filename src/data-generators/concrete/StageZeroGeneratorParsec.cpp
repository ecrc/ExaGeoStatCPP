// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StageZeroGeneratorParsec.cpp
 * @brief Implementation of the StageZeroGeneratorParsec class using PaRSEC/DPLASMA
 * @version 2.0.0
**/

#include <data-generators/concrete/StageZeroGeneratorParsec.hpp>
#include <data-generators/LocationGenerator.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <pnetcdf.h>
#include <nlopt.h>
#include <mpi.h>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <thread>
#include <vector>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define START_TIMING(x) double _start_##x = MPI_Wtime();
#define STOP_TIMING(x) x = MPI_Wtime() - _start_##x;

namespace {
constexpr int kForcingOffset = 238;
}

using namespace exageostat::generators::stagezero;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::results;

template<typename T>
StageZeroGeneratorParsec<T> *StageZeroGeneratorParsec<T>::GetInstance() {
    if (mpInstance == nullptr) {
        // Ensure MPI is initialized before creating the instance
        // int initialized;
        // MPI_Initialized(&initialized);
        // if (!initialized) {
        //     int argc = 0;
        //     char **argv = nullptr;
        //     MPI_Init(&argc, &argv);
        // }
        mpInstance = new StageZeroGeneratorParsec<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
StageZeroGeneratorParsec<T>::CreateData(exageostat::configurations::Configurations &aConfigurations,
                                        exageostat::kernels::Kernel<T> &aKernel) {
    this->Runner(aConfigurations);
    return std::move(this->mData);
}

template<typename T>
void StageZeroGeneratorParsec<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

template<typename T>
void StageZeroGeneratorParsec<T>::Runner(exageostat::configurations::Configurations &aConfigurations) {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    mArgs.mConfigs = &aConfigurations;
    this->ConfigureGenerator();
    this->Allocate();
    this->ReadForcingData();
    this->ReadNetCDFFiles();
    
    // After data gathering, only rank 0 continues with optimization
    // Other ranks can exit since they're no longer needed
    if (rank != 0) {
        fprintf(stderr, "[StageZero PaRSEC] Rank %d: Data gathering complete, exiting (only rank 0 performs optimization)\n", rank);
        this->CleanUp();
        return;
    }
    
    fprintf(stderr, "[StageZero PaRSEC] Rank 0: Starting optimization phase\n");
    this->RunMeanTrend();
    
    this->CleanUp();
}

template<typename T>
void StageZeroGeneratorParsec<T>::ConfigureGenerator() {
    
    // Values for the mean trend removal
    mArgs.mM = 10;
    mArgs.mT = 365*24;
    mArgs.mNoYears = 751;
    
    // Only support trend_model for this workflow
    if (mArgs.mConfigs->GetKernelName() == "TrendModel" || mArgs.mConfigs->GetKernelName() == "trend_model") {
        mArgs.mNumParams = 1;
    } else {
        throw std::invalid_argument("Unsupported kernel for Stage Zero: only 'TrendModel' is supported");
    }
    
    // Derive observation years and N from configuration
    int start_year = mArgs.mConfigs->GetStartYear();
    int end_year   = mArgs.mConfigs->GetEndYear();
    if (end_year < start_year) {
        throw std::runtime_error("EndYear must be >= StartYear");
    }
    int obs_years = (end_year - start_year + 1);
    mArgs.mN = static_cast<size_t>(mArgs.mT) * static_cast<size_t>(obs_years);
    
    // Number of locations
    try { mArgs.mNumLocs = mArgs.mConfigs->GetNumLocs(); }
    catch (...) { mArgs.mNumLocs = mArgs.mConfigs->GetProblemSize(); }
}

template<typename T>
void StageZeroGeneratorParsec<T>::ReadNetCDFFiles() {

    int ncid, retval;
    MPI_Offset lat_len, lon_len, time_len;
    int lon_varid, lat_varid, time_varid, t2m_varid;
    int v = 365;
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int start_year = mArgs.mConfigs->GetStartYear();
    int end_year   = mArgs.mConfigs->GetEndYear();

    auto openFileNCmpi = [&](char *filename) {
        int id, ret;
        if ((ret = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &id))) {
            throw std::runtime_error("Error opening NetCDF file: " + std::string(ncmpi_strerror(ret)));
        }
        return id;
    };

    auto closeFileNCmpi = [&](int id) {
        int ret;
        if ((ret = ncmpi_close(id))) {
            throw std::runtime_error("Error closing NetCDF file: " + std::string(ncmpi_strerror(ret)));
        }
    };

    
    char path[256];
    snprintf(path, sizeof(path), "%sdata_%d.nc", mArgs.mConfigs->GetDataPath().c_str(), start_year);

    ncid = openFileNCmpi(path);

    // Dimension IDs and lengths with error checks
    if ((retval = ncmpi_inq_dimid(ncid, "longitude", &lon_varid)))
        std::cerr << "Error: " << ncmpi_strerror(retval) << std::endl;
    if ((retval = ncmpi_inq_dimlen(ncid, lon_varid, &lon_len)))
        std::cerr << "Error: " << ncmpi_strerror(retval) << std::endl;

    if ((retval = ncmpi_inq_dimid(ncid, "latitude", &lat_varid)))
        std::cerr << "Error: " << ncmpi_strerror(retval) << std::endl;
    if ((retval = ncmpi_inq_dimlen(ncid, lat_varid, &lat_len)))
        std::cerr << "Error: " << ncmpi_strerror(retval) << std::endl;

    if ((retval = ncmpi_inq_dimid(ncid, "time", &time_varid)))
        std::cerr << "Error: " << ncmpi_strerror(retval) << std::endl;

    if ((retval = ncmpi_inq_varid(ncid, "t2m", &t2m_varid)))
        std::cerr << "Error: " << ncmpi_strerror(retval) << std::endl;

    closeFileNCmpi(ncid);

    int min_loc = mArgs.mConfigs->GetLowTileSize(); 
    
    for (int y = start_year; y <= end_year; y++) {
        char path2[256];
        snprintf(path2, sizeof(path2), "%sdata_%d.nc", mArgs.mConfigs->GetDataPath().c_str(), y);

        ncid = openFileNCmpi(path2);

        double scaling_var = 1.0, offset_var = 0.0;
        if ((retval = ncmpi_get_att_double(ncid, t2m_varid, "scale_factor", &scaling_var))) {
            std::cerr << "Warning: missing scale_factor: " << ncmpi_strerror(retval) << std::endl;
        }
        if ((retval = ncmpi_get_att_double(ncid, t2m_varid, "add_offset", &offset_var))) {
            std::cerr << "Warning: missing add_offset: " << ncmpi_strerror(retval) << std::endl;
        }

        int *t2m = nullptr;
        int *t2m_local = nullptr;

        if (IsLeapYear(y)) time_len = v * 24 + 24;
        else time_len = v * 24;

        size_t total_elems = static_cast<size_t>(time_len) * static_cast<size_t>(mArgs.mNumLocs);
        size_t local_elems = static_cast<size_t>(time_len / nprocs) * static_cast<size_t>(mArgs.mNumLocs);
        t2m = (int *) malloc(total_elems * sizeof(int));
        t2m_local = (int *) malloc(local_elems * sizeof(int));

        MPI_Offset index[] = {(MPI_Offset) (time_len / nprocs) * rank,
                              (MPI_Offset) min_loc, 0};
        MPI_Offset count[] = {(MPI_Offset) (time_len / nprocs), 1, (MPI_Offset) mArgs.mNumLocs};
        
        double x_read = 0;
        START_TIMING(x_read)
        ncmpi_get_vara_int_all(ncid, t2m_varid, index, count, t2m_local);
        STOP_TIMING(x_read)

        double gather_start = MPI_Wtime();
        MPI_Allgather(t2m_local, (int) local_elems, MPI_INT, t2m,
                      (int) local_elems, MPI_INT, MPI_COMM_WORLD);
        double gather_time = MPI_Wtime() - gather_start;
        
        fprintf(stderr, "[StageZero] Year %d: Read time: %.3f sec, Gather time: %.3f sec\n", y, x_read, gather_time);

        for (int lu = 0; lu < mArgs.mNumLocs; lu++) {
            int r = 0;
            for (size_t k = lu; k < (size_t) time_len * mArgs.mNumLocs; k += mArgs.mNumLocs) {
                if (!(IsLeapYear(y) && (r >= 1416 && r <= 1439))) {
                    mArgs.mT2mHourlyPerYear[lu][mArgs.mT2mHourlyPerYearCount[lu]++] =
                            (static_cast<double>(t2m[k]) * scaling_var) + offset_var - 273.15;
                }
                r++;
            }
        }
        closeFileNCmpi(ncid);
        free(t2m);
        free(t2m_local);
        fprintf(stderr, "[StageZero] Completed processing year %d\n\n", y);
    }
    fprintf(stderr, "[StageZero] NetCDF loading completed for years %d-%d (%d files)\n\n", 
            start_year, end_year, end_year - start_year + 1);
}

template<typename T>
void StageZeroGeneratorParsec<T>::ReadForcingData() {
    std::string forcing_path = mArgs.mConfigs->GetForcingDataPath();
    if (!forcing_path.empty() && forcing_path.back() == '/') forcing_path.pop_back();
    
    fprintf(stderr, "[StageZero] Loading forcing data from: %s\n", forcing_path.c_str());
    mArgs.mForcing = ReadObsFile((char *) mArgs.mConfigs->GetForcingDataPath().c_str(), mArgs.mNoYears);
    fprintf(stderr, "[StageZero] Successfully loaded forcing data (%d years)\n\n", mArgs.mNoYears);
}

template<typename T>
void StageZeroGeneratorParsec<T>::Allocate() {

    mArgs.mStartingTheta = new double[mArgs.mNumParams];
    mArgs.mTargetTheta = new double[mArgs.mNumParams];
    mArgs.mInitialTheta = new double[mArgs.mNumParams];
    mArgs.mLb = new double[mArgs.mNumParams];
    mArgs.mUp = new double[mArgs.mNumParams];

    // Copy configuration values to allocated arrays
    auto lb_vec = mArgs.mConfigs->GetLowerBounds();
    auto up_vec = mArgs.mConfigs->GetUpperBounds();
    auto starting_theta_vec = mArgs.mConfigs->GetStartingTheta();

    for(size_t i = 0; i < mArgs.mNumParams; ++i) {
        mArgs.mLb[i] = lb_vec[i];
        mArgs.mUp[i] = up_vec[i];
        mArgs.mStartingTheta[i] = starting_theta_vec[i];
    }

    mArgs.mT2mHourlyPerYear = (double **) malloc(mArgs.mNumLocs * sizeof(double *));
    mArgs.mT2mHourlyPerYearCount = (int *) calloc(mArgs.mNumLocs, sizeof(int));
    for (int k = 0; k < mArgs.mNumLocs; k++) {
        mArgs.mT2mHourlyPerYear[k] = (double *) malloc(mArgs.mN * sizeof(double));
    }
}

template<typename T>
void StageZeroGeneratorParsec<T>::RunMeanTrend() {
    const int no_locs = mArgs.mNumLocs;

    // Print configuration summary
    const std::string kernel = mArgs.mConfigs->GetKernelName();
    const std::string data_path = mArgs.mConfigs->GetDataPath();
    const std::string forcing_path = mArgs.mConfigs->GetForcingDataPath();
    const std::string results_path = mArgs.mConfigs->GetResultsPath();
    const int start_year = mArgs.mConfigs->GetStartYear();
    const int end_year = mArgs.mConfigs->GetEndYear();
    const int years = mArgs.mNoYears;
    const int period_hours = mArgs.mT;
    const int M = mArgs.mM;
    const size_t N = mArgs.mN;
    const int num_locs = mArgs.mNumLocs;
    const double tol = mArgs.mConfigs->GetTolerance();
    const int maxeval = mArgs.mConfigs->GetMaxMleIterations();
    const double st = mArgs.mStartingTheta[0];
    const double lb0 = mArgs.mLb[0];
    const double ub0 = mArgs.mUp[0];
    const int dts = mArgs.mConfigs->GetDenseTileSize();
    const int lts_val = mArgs.mConfigs->GetLowTileSize();
    
    fprintf(stderr, "----- StageZero Arguments -----\n");
    fprintf(stderr, "kernel: %s\n", kernel.c_str());
    fprintf(stderr, "data_path: %s\n", data_path.c_str());
    fprintf(stderr, "forcing_data_path: %s\n", forcing_path.c_str());
    fprintf(stderr, "results_path: %s\n", results_path.c_str());
    fprintf(stderr, "start_year: %d, end_year: %d, years: %d\n", start_year, end_year, end_year-start_year+1);
    fprintf(stderr, "num_locs: %d\n", num_locs);
    fprintf(stderr, "period_hours(T): %d, harmonics(M): %d\n", period_hours, M);
    fprintf(stderr, "N (observations): %zu\n", N);
    fprintf(stderr, "tolerance(ftol_abs): %.10e, maxeval: %d\n", tol, maxeval);
    fprintf(stderr, "-------------------------------\n");

    // Key bounds info
    fprintf(stderr, "Starting theta[0]: %f\n", mArgs.mStartingTheta[0]);
    fprintf(stderr, "Lower bound: %f, Upper bound: %f\n", mArgs.mLb[0], mArgs.mUp[0]);
    
    this->SetupMLEComponents();
    
    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, mArgs.mNumParams);
    nlopt_set_lower_bounds(opt, mArgs.mLb);
    nlopt_set_upper_bounds(opt, mArgs.mUp);
    nlopt_set_max_objective(opt, &StageZeroGeneratorParsec<T>::StageZeroObjectiveCallback, (void *)this);
    
    nlopt_set_ftol_abs(opt, mArgs.mConfigs->GetTolerance());
    nlopt_set_maxeval(opt, mArgs.mConfigs->GetMaxMleIterations());
    
    std::vector<double> theta(mArgs.mNumParams, 0.0);
    
    for (int l = 0; l < no_locs; ++l) {
        mArgs.mIterCount = 0;
        mArgs.mCurrentLocation = l;

        this->ConvertT2MToZForLocation(l);

        fprintf(stderr, "[StageZero PaRSEC] Starting NLopt for location %d/%d\n", l + 1, no_locs);
        theta[0] = mArgs.mStartingTheta[0];
        double opt_f = 0.0;

        nlopt_set_max_objective(opt, &StageZeroGeneratorParsec<T>::StageZeroObjectiveCallback, (void *)this);
        nlopt_result nlres = nlopt_optimize(opt, theta.data(), &opt_f);
        fprintf(stderr, "[StageZero PaRSEC] NLopt finished (res=%d), theta=%0.10f, f=%0.10f\n", (int) nlres, theta[0], opt_f);

        {
            std::vector<double> grad;
            this->MLEAlgorithm(theta, grad, nullptr);
        }
    }
    
    nlopt_destroy(opt);
}

template<typename T>
void StageZeroGeneratorParsec<T>::CleanUp() {

    delete[] mArgs.mStartingTheta;
    delete[] mArgs.mTargetTheta;
    delete[] mArgs.mInitialTheta;
    delete[] mArgs.mLb;
    delete[] mArgs.mUp;
    delete[] mArgs.mForcing;

    for (int i = 0; i < mArgs.mNumLocs; ++i) {
        free(mArgs.mT2mHourlyPerYear[i]);
    }
    free(mArgs.mT2mHourlyPerYear);
    free(mArgs.mT2mHourlyPerYearCount);
    
    // Clean up PaRSEC descriptors
    if (mArgs.mpDescZ) {
        parsec_data_free(mArgs.mpDescZ->mat);
        free(mArgs.mpDescZ);
        mArgs.mpDescZ = nullptr;
    }
    if (mArgs.mpX) {
        parsec_data_free(mArgs.mpX->mat);
        free(mArgs.mpX);
        mArgs.mpX = nullptr;
    }
    if (mArgs.mpXtX) {
        parsec_data_free(mArgs.mpXtX->mat);
        free(mArgs.mpXtX);
        mArgs.mpXtX = nullptr;
    }
    if (mArgs.mpDescPart1) {
        parsec_data_free(mArgs.mpDescPart1->mat);
        free(mArgs.mpDescPart1);
        mArgs.mpDescPart1 = nullptr;
    }
    if (mArgs.mpDescPart2) {
        parsec_data_free(mArgs.mpDescPart2->mat);
        free(mArgs.mpDescPart2);
        mArgs.mpDescPart2 = nullptr;
    }
    if (mArgs.mpPart2Vector) {
        parsec_data_free(mArgs.mpPart2Vector->mat);
        free(mArgs.mpPart2Vector);
        mArgs.mpPart2Vector = nullptr;
    }
    if (mArgs.mpEstimatedMeanTrend) {
        parsec_data_free(mArgs.mpEstimatedMeanTrend->mat);
        free(mArgs.mpEstimatedMeanTrend);
        mArgs.mpEstimatedMeanTrend = nullptr;
    }
}

template<typename T>
void StageZeroGeneratorParsec<T>::SetupMLEComponents() {
    
    int N = mArgs.mN;
    int nparams = 3 + 2*mArgs.mM;
    int dts = mArgs.mConfigs->GetDenseTileSize();
    int p_grid = 1, q_grid = 1;
    
    fprintf(stderr, "Creating PaRSEC descriptors: N=%d, nparams=%d, dts=%d\n", N, nparams, dts);
    
    mArgs.mIterCount = 0;
    
    // Get PaRSEC context from hardware layer
    try {
        mArgs.mpParsecContext = static_cast<parsec_context_t*>(ExaGeoStatHardware::GetParsecContext());
        if (mArgs.mpParsecContext) {
            fprintf(stderr, "PaRSEC context obtained successfully\n");
        } else {
            fprintf(stderr, "Warning: PaRSEC context is null, using fallback computations\n");
        }
    } catch (const std::exception& e) {
        fprintf(stderr, "Warning: Could not get PaRSEC context: %s, using fallback computations\n", e.what());
        mArgs.mpParsecContext = nullptr;
    }
    
    // Create PaRSEC block-cyclic descriptors using proper initialization
    int rank = 0; // Single process for now
    int nodes = 1; // Single node for now
    
    // Use same tile size as CHAMELEON version for exact equivalence
    int tile_size = dts; // Use full dense tile size like CHAMELEON
    
    // Z vector (observations) - N x 1
    mArgs.mpDescZ = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpDescZ, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank, 
                                   tile_size, tile_size, N, 1, 0, 0, N, 1, 1, nodes, 1, 1, 0, 0);
    mArgs.mpDescZ->mat = parsec_data_allocate((size_t)mArgs.mpDescZ->super.nb_local_tiles *
                                             (size_t)mArgs.mpDescZ->super.bsiz *
                                             (size_t)parsec_datadist_getsizeoftype(mArgs.mpDescZ->super.mtype));
    
    // X matrix (design matrix) - N x nparams
    mArgs.mpX = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpX, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   tile_size, tile_size, N, nparams, 0, 0, N, nparams, 1, nodes, 1, 1, 0, 0);
    mArgs.mpX->mat = parsec_data_allocate((size_t)mArgs.mpX->super.nb_local_tiles *
                                         (size_t)mArgs.mpX->super.bsiz *
                                         (size_t)parsec_datadist_getsizeoftype(mArgs.mpX->super.mtype));
    
    // XtX matrix - nparams x nparams
    mArgs.mpXtX = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpXtX, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   tile_size, tile_size, nparams, nparams, 0, 0, nparams, nparams, 1, nodes, 1, 1, 0, 0);
    mArgs.mpXtX->mat = parsec_data_allocate((size_t)mArgs.mpXtX->super.nb_local_tiles *
                                           (size_t)mArgs.mpXtX->super.bsiz *
                                           (size_t)parsec_datadist_getsizeoftype(mArgs.mpXtX->super.mtype));
    
    // part1 scalar - 1 x 1
    mArgs.mpDescPart1 = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpDescPart1, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   tile_size, tile_size, 1, 1, 0, 0, 1, 1, 1, nodes, 1, 1, 0, 0);
    mArgs.mpDescPart1->mat = parsec_data_allocate((size_t)mArgs.mpDescPart1->super.nb_local_tiles *
                                                 (size_t)mArgs.mpDescPart1->super.bsiz *
                                                 (size_t)parsec_datadist_getsizeoftype(mArgs.mpDescPart1->super.mtype));
    
    // part2 scalar - 1 x 1
    mArgs.mpDescPart2 = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpDescPart2, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   tile_size, tile_size, 1, 1, 0, 0, 1, 1, 1, nodes, 1, 1, 0, 0);
    mArgs.mpDescPart2->mat = parsec_data_allocate((size_t)mArgs.mpDescPart2->super.nb_local_tiles *
                                                 (size_t)mArgs.mpDescPart2->super.bsiz *
                                                 (size_t)parsec_datadist_getsizeoftype(mArgs.mpDescPart2->super.mtype));
    
    // part2_vector - nparams x 1
    mArgs.mpPart2Vector = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpPart2Vector, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   tile_size, tile_size, nparams, 1, 0, 0, nparams, 1, 1, nodes, 1, 1, 0, 0);
    mArgs.mpPart2Vector->mat = parsec_data_allocate((size_t)mArgs.mpPart2Vector->super.nb_local_tiles *
                                                   (size_t)mArgs.mpPart2Vector->super.bsiz *
                                                   (size_t)parsec_datadist_getsizeoftype(mArgs.mpPart2Vector->super.mtype));
    
    // estimated_mean_trend vector - N x 1
    mArgs.mpEstimatedMeanTrend = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
    parsec_matrix_block_cyclic_init(mArgs.mpEstimatedMeanTrend, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   tile_size, tile_size, N, 1, 0, 0, N, 1, 1, nodes, 1, 1, 0, 0);
    mArgs.mpEstimatedMeanTrend->mat = parsec_data_allocate((size_t)mArgs.mpEstimatedMeanTrend->super.nb_local_tiles *
                                                          (size_t)mArgs.mpEstimatedMeanTrend->super.bsiz *
                                                          (size_t)parsec_datadist_getsizeoftype(mArgs.mpEstimatedMeanTrend->super.mtype));
    
    fprintf(stderr, "All PaRSEC descriptors created successfully\n");
}

template<typename T>
void StageZeroGeneratorParsec<T>::ConvertT2MToZForLocation(int location_index) {
    
    parsec_matrix_block_cyclic_t *DESC_Z = mArgs.mpDescZ;
    int N = mArgs.mN;

    if (location_index < 0 || location_index >= mArgs.mNumLocs) {
        throw std::runtime_error("ConvertT2MToZForLocation: invalid location index");
    }

    // Copy data directly to PaRSEC descriptor
    double *z_data = (double*)DESC_Z->mat;
    int count = std::min(N, mArgs.mT2mHourlyPerYearCount[location_index]);
    int mb = DESC_Z->super.mb;
    int nb = DESC_Z->super.nb;
    int mt = (N + mb - 1) / mb;
    int nt = (1 + nb - 1) / nb;
    int bsiz = DESC_Z->super.bsiz;
    for (int tj = 0; tj < nt; ++tj) {
        for (int ti = 0; ti < mt; ++ti) {
            double *tile_ptr = z_data + (tj * mt + ti) * bsiz;
            int row_offset = ti * mb;
            int col_offset = tj * nb; // vector has only one column
            int rows = std::min(mb, N - row_offset);
            int cols = std::min(nb, 1 - col_offset);
            if (cols <= 0) continue;
            for (int cj = 0; cj < cols; ++cj) {
                for (int ri = 0; ri < rows; ++ri) {
                    int global_row = row_offset + ri;
                    double val = (global_row < count) ? mArgs.mT2mHourlyPerYear[location_index][global_row] : 0.0;
                    tile_ptr[ri + cj * mb] = val;
                }
            }
        }
    }
}

template<typename T>
void StageZeroGeneratorParsec<T>::GenerateDesignMatrixExact(double *matrix, int m, int n, int m0, int n0, double *localtheta) {
    // localtheta layout: [theta, T, M, forcing...]
    const double theta = localtheta[0];
    const int    T_val = static_cast<int>(localtheta[1]);
    const int    M_val = static_cast<int>(localtheta[2]);
    double *forcing = &localtheta[3];

    int i0 = m0;
    int j0_start = n0;
    double i_x = static_cast<double>(i0) + 1.0;

    // Precompute the AR-term recursively: r[0] = 0; r[n+1] = theta*r[n] + (1-theta)*forcing[n]
    std::vector<double> ar_acc(static_cast<size_t>(mArgs.mNoYears) + 1, 0.0);
    for (int n = 0; n < mArgs.mNoYears; ++n) {
        ar_acc[static_cast<size_t>(n + 1)] = theta * ar_acc[static_cast<size_t>(n)] + (1.0 - theta) * forcing[n];
    }

    for (int i = 0; i < m; i++) {
        int j0 = j0_start;
        for (int j = 0; j < n; j++) {
            const int idx = i + j * m; // column-major
            if (j0 == 0) {
                // First column: constant 1
                matrix[idx] = 1.0;
            } else if (j0 == 1) {
                // Second column: forcing[(i0/T) + OFFSET]
                const int ty = (i0 / T_val);
                const int f_idx = ty + kForcingOffset;
                matrix[idx] = forcing[std::min(f_idx, mArgs.mNoYears - 1)];
            } else if (j0 == 2) {
                // Third column: use recursive accumulation r[n]
                const int ty = (i0 / T_val);
                const int n_idx = std::min(ty + kForcingOffset, mArgs.mNoYears);
                matrix[idx] = ar_acc[static_cast<size_t>(n_idx)];
            } else {
                // Remaining 2*M columns: alternate sin/cos per C core (j even -> sin, j odd -> cos)
                // frequency index = floor((j0-3)/2) + 1  in [1..M]
                const int freq = static_cast<int>(std::floor((j0 - 3.0) / 2.0)) + 1;
                const double angle = 2.0 * M_PI * (i_x) * static_cast<double>(freq) / static_cast<double>(T_val);
                if ((j % 2) == 0) {
                    matrix[idx] = std::sin(angle);
                } else {
                    matrix[idx] = std::cos(angle);
                }
            }
            j0++;
        }
        i0++;
        i_x += 1.0;
    }
}

template<typename T>
double StageZeroGeneratorParsec<T>::MLEAlgorithm(const std::vector<double> &aThetaVec,
                                               std::vector<double> &aGrad, void *apObj) {
    
    parsec_matrix_block_cyclic_t *Zobs = mArgs.mpDescZ;
    parsec_matrix_block_cyclic_t *X = mArgs.mpX;
    parsec_matrix_block_cyclic_t *XtX = mArgs.mpXtX;
    parsec_matrix_block_cyclic_t *part1 = mArgs.mpDescPart1;
    parsec_matrix_block_cyclic_t *part2 = mArgs.mpDescPart2;
    parsec_matrix_block_cyclic_t *part2_vector = mArgs.mpPart2Vector;
    int N = static_cast<int>(mArgs.mN);
    double value = 0.0;
    const bool is_final_call = (apObj == nullptr);
    
    if (!X) { throw std::runtime_error("X descriptor is null"); }
    if (!Zobs) { throw std::runtime_error("Zobs descriptor is null"); }
    
    double* localtheta = (double *) malloc((3+mArgs.mNoYears) * sizeof(double));
    localtheta[0] = aThetaVec[0];
    localtheta[1] = mArgs.mT;
    localtheta[2] = mArgs.mM;
    
    for(int ii=0; ii<mArgs.mNoYears; ii++)
        localtheta[3+ii] = mArgs.mForcing[ii];
    
    try {
        // Step 1: generate design matrix X
        {
            int Nrows = static_cast<int>(mArgs.mN);
            int Ncols = 3 + 2*mArgs.mM;
            
            // Create LAPACK buffer like CHAMELEON version
            double *x_lap = (double*)calloc((size_t)Nrows * (size_t)Ncols, sizeof(double));
            if (!x_lap) { throw std::runtime_error("Allocation failed for design matrix buffer"); }
            
            this->GenerateDesignMatrixExact(x_lap, Nrows, Ncols, 0, 0, localtheta);
            
            // Tile-aware copy from LAPACK (column-major) to PaRSEC tiled layout
            double *x_data = (double*)X->mat;
            int mb = X->super.mb;
            int nb = X->super.nb;
            int mt = (Nrows + mb - 1) / mb;
            int nt = (Ncols + nb - 1) / nb;
            int bsiz = X->super.bsiz;
            for (int tj = 0; tj < nt; ++tj) {
                for (int ti = 0; ti < mt; ++ti) {
                    double *tile_ptr = x_data + (tj * mt + ti) * bsiz;
                    int row_offset = ti * mb;
                    int col_offset = tj * nb;
                    int rows = std::min(mb, Nrows - row_offset);
                    int cols = std::min(nb, Ncols - col_offset);
                    for (int cj = 0; cj < cols; ++cj) {
                        for (int ri = 0; ri < rows; ++ri) {
                            int global_row = row_offset + ri;
                            int global_col = col_offset + cj;
                            tile_ptr[ri + cj * mb] = x_lap[global_row + global_col * Nrows];
                        }
                    }
                }
            }
            
            free(x_lap);
        }
        
        // Ensure workspaces are clean before computations (equivalent to CHAMELEON_dlaset_Tile)
        dplasma_dlaset(mArgs.mpParsecContext, dplasmaUpperLower, 0.0, 0.0, (parsec_tiled_matrix_t*)part1);
        dplasma_dlaset(mArgs.mpParsecContext, dplasmaUpperLower, 0.0, 0.0, (parsec_tiled_matrix_t*)part2);
        dplasma_dlaset(mArgs.mpParsecContext, dplasmaUpperLower, 0.0, 0.0, (parsec_tiled_matrix_t*)part2_vector);
        dplasma_dlaset(mArgs.mpParsecContext, dplasmaUpperLower, 0.0, 0.0, (parsec_tiled_matrix_t*)XtX);

        // Step 2: part1 = Z^T * Z (using DPLASMA)
        dplasma_dgemm(mArgs.mpParsecContext, dplasmaTrans, dplasmaNoTrans,
                     1.0, (parsec_tiled_matrix_t*)Zobs, (parsec_tiled_matrix_t*)Zobs, 0.0, (parsec_tiled_matrix_t*)part1);
        
        // Optional debug removed
        
        // Step 3: part2_vector = X^T * Z (using DPLASMA)
        dplasma_dgemm(mArgs.mpParsecContext, dplasmaTrans, dplasmaNoTrans,
                         1.0, (parsec_tiled_matrix_t*)X, (parsec_tiled_matrix_t*)Zobs, 0.0, (parsec_tiled_matrix_t*)part2_vector);
        
        // Optional debug removed
    
        // Step 4: XtX = X^T * X (using DPLASMA)
            dplasma_dgemm(mArgs.mpParsecContext, dplasmaTrans, dplasmaNoTrans,
                         1.0, (parsec_tiled_matrix_t*)X, (parsec_tiled_matrix_t*)X, 0.0, (parsec_tiled_matrix_t*)XtX);

        // Step 5: Cholesky decomposition (using DPLASMA)
        int info = dplasma_dpotrf(mArgs.mpParsecContext, dplasmaLower, (parsec_tiled_matrix_t*)XtX);
        
        if(info != 0) {
            if (apObj) { free(localtheta); return -1e18; }
            else {
                fprintf(stderr, "[StageZero PaRSEC] ERROR: Cholesky failed. Aborting.\n");
                free(localtheta);
                exit(1);
            }
        }
    
        // Step 6: Triangular solve
            dplasma_dtrsm(mArgs.mpParsecContext, dplasmaLeft, dplasmaLower, dplasmaNoTrans, dplasmaNonUnit,
                         1.0, (parsec_tiled_matrix_t*)XtX, (parsec_tiled_matrix_t*)part2_vector);
        
        if (!apObj) {
            dplasma_dtrsm(mArgs.mpParsecContext, dplasmaLeft, dplasmaLower, dplasmaTrans, dplasmaNonUnit,
                         1.0, (parsec_tiled_matrix_t*)XtX, (parsec_tiled_matrix_t*)part2_vector);
        }

        // If this is an optimization callback, compute and return the C objective using y^T y
        if (apObj) {
            double part1_val_cb = 0.0;
            if (part1 && part1->mat) part1_val_cb = static_cast<double *>(part1->mat)[0];
            // Compute part2 = y^T y (equivalent to CHAMELEON_dgemm_Tile)
            dplasma_dlaset(mArgs.mpParsecContext, dplasmaUpperLower, 0.0, 0.0, (parsec_tiled_matrix_t*)part2);
            dplasma_dgemm(mArgs.mpParsecContext, dplasmaTrans, dplasmaNoTrans, 1.0, 
                         (parsec_tiled_matrix_t*)part2_vector, (parsec_tiled_matrix_t*)part2_vector, 0.0, (parsec_tiled_matrix_t*)part2);
            double part2_val_cb = 0.0;
            if (part2 && part2->mat) part2_val_cb = static_cast<double *>(part2->mat)[0];
            
            // If numerical issues cause part1 <= part2, penalize
            if (part1_val_cb <= part2_val_cb) {
                free(localtheta);
                return -1e18;
            }
            value = (-1.0) * std::log(part1_val_cb - part2_val_cb);
            
            mArgs.mIterCount++;
            free(localtheta);
            return value;
        }

        // Step 7: Final computation with CSV output
        // Create estimated_mean_trend = X * beta (using part2_vector as beta)
        parsec_matrix_block_cyclic_t *estimated_mean_trend = mArgs.mpEstimatedMeanTrend;
        
            dplasma_dgemm(mArgs.mpParsecContext, dplasmaNoTrans, dplasmaNoTrans,
                         1.0, (parsec_tiled_matrix_t*)X, (parsec_tiled_matrix_t*)part2_vector, 0.0, (parsec_tiled_matrix_t*)estimated_mean_trend);
        
        // Residuals = Z - trend (stored back into estimated_mean_trend) - equivalent to CHAMELEON_dgeadd_Tile
        dplasma_dgeadd(mArgs.mpParsecContext, dplasmaNoTrans, 1.0, (parsec_tiled_matrix_t*)Zobs, -1.0, (parsec_tiled_matrix_t*)estimated_mean_trend);
        
        // Compute sigma² = residuals^T * residuals (equivalent to CHAMELEON_dgemm_Tile)
        dplasma_dgemm(mArgs.mpParsecContext, dplasmaTrans, dplasmaNoTrans, 1.0, 
                     (parsec_tiled_matrix_t*)estimated_mean_trend, (parsec_tiled_matrix_t*)estimated_mean_trend, 0.0, (parsec_tiled_matrix_t*)part2);
        
        // Get sigma² value and divide by N (like CHAMELEON version)
        double sigma_squared_raw = static_cast<double *>(part2->mat)[0];
        double sigma_squared = sigma_squared_raw / N;
        
        // Standardize residuals: residuals = residuals / sqrt(sigma²)
        double sqrt_sigma = std::sqrt(sigma_squared);
        double *trend_data = (double*)estimated_mean_trend->mat;
        for (int i = 0; i < N; ++i) {
            trend_data[i] /= sqrt_sigma;
        }
        
        // Set objective value (maximize -sigma^2 -> minimize sigma^2)
        value = -sigma_squared;

        // Static global arrays for multi-location storage (like CHAMELEON version)
        static std::vector<std::vector<double>> Z_new(mArgs.mNumLocs, std::vector<double>(N));
        static std::vector<std::vector<double>> params(mArgs.mNumLocs, std::vector<double>(3 + 2*mArgs.mM + 2));
        
        // Store normalized residuals for this location
        int location_index = mArgs.mCurrentLocation;  // Current location being processed (0-based)
        for (int i = 0; i < N; ++i) {
            Z_new[location_index][i] = trend_data[i];
        }
        
        // Store parameters
        params[location_index][0] = aThetaVec[0];  // optimized theta
        params[location_index][1] = sigma_squared;  // sigma²
        double *part2_vector_data_local = (double*)part2_vector->mat;
        for(int i = 0; i < 3 + 2*mArgs.mM; i++) {
            params[location_index][i+2] = part2_vector_data_local[i];
        }
        
        // Write CSV files only on the final call and when processing the last location
        if (is_final_call && (location_index == mArgs.mNumLocs - 1)) {
            std::string results_path;
            try { 
                results_path = mArgs.mConfigs->GetResultsPath(); 
            } catch (...) { 
                results_path.clear(); 
            }
            
            if (results_path.empty()) { 
                fprintf(stderr, "[StageZero PaRSEC] ERROR: ResultsPath is required. Please set --resultspath.\n");
                throw std::runtime_error("ResultsPath is required. Please set --resultspath.");
            }
            
            if (results_path.back() != '/') results_path += "/";
            
            // Create results directory if it doesn't exist
            try {
                if (!std::filesystem::exists(results_path)) {
                    fprintf(stderr, "[StageZero PaRSEC] Creating results directory: %s\n", results_path.c_str());
                    std::filesystem::create_directories(results_path);
                }
            } catch (const std::exception &e) {
                fprintf(stderr, "[StageZero PaRSEC] ERROR: Failed to create output directory: %s (%s)\n", results_path.c_str(), e.what());
                throw std::runtime_error("Failed to create output directory: " + std::string(e.what()));
            }
            
            fprintf(stderr, "[StageZero PaRSEC] Writing outputs to: %s\n", results_path.c_str());
            
            // Write Z files for each time slot (replace mode instead of append)
            #pragma omp parallel for
            for (int time_slot = 0; time_slot < N; time_slot++) {
                char file_path[256];
                std::snprintf(file_path, sizeof(file_path), "%sz_%d.csv", results_path.c_str(), time_slot);
                
                // Retry loop for file locking with replace mode
                while (true) {
                    FILE *fp = std::fopen(file_path, "w");  // Replace mode instead of append
                    if (!fp) {
                        fprintf(stderr, "[StageZero PaRSEC] WARNING: Failed to open file %s, retrying...\n", file_path);
                        std::this_thread::sleep_for(std::chrono::seconds(1));
                        continue;
                    }
                    
                    // Write all locations for this time slot
                    for (int loc = 0; loc < mArgs.mNumLocs; ++loc) {
                        std::fprintf(fp, "%.14f\n", Z_new[loc][time_slot]);
                    }
                    std::fclose(fp);
                    break;
                }
            }
            
            // Write params file (replace mode)
            char params_file_path[256];
            std::snprintf(params_file_path, sizeof(params_file_path), "%sparams.csv", results_path.c_str());
            FILE* fp = std::fopen(params_file_path, "w");  // Replace mode
            if(fp == NULL){ 
                fprintf(stderr, "[StageZero PaRSEC] ERROR: Failed to create params file: %s\n", params_file_path);
                exit(1); 
            }
            
            for(int k = 0; k < mArgs.mNumLocs; k++) {
                for(int i = 0; i < 3 + 2*mArgs.mM + 2; i++) {
                    fprintf(fp, "%.14f ", params[k][i]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            
            fprintf(stderr, "[StageZero PaRSEC] Successfully wrote %d Z files and params.csv to %s\n", N, results_path.c_str());
        }

        // Set summary information for final output
        if (is_final_call) {
            Results::GetInstance()->SetGeneratedLocationsNumber(mArgs.mNumLocs);
            Results::GetInstance()->SetIsLogger(mArgs.mConfigs->GetLogger());
            Results::GetInstance()->SetLoggerPath(mArgs.mConfigs->GetLoggerPath());
        }

        mArgs.mIterCount++;
        
    } catch (const std::exception& e) {
        free(localtheta);
        throw;
    }
    free(localtheta);
    return value;
}

template<typename T>
double StageZeroGeneratorParsec<T>::StageZeroObjectiveCallback(unsigned aN, const double *aTheta, double *aGrad, void *aData) {
    auto *generator = static_cast<StageZeroGeneratorParsec<T> *>(aData);
    std::vector<double> theta_vec(aTheta, aTheta + aN);
    std::vector<double> grad_vec(aN);
    double result = generator->MLEAlgorithm(theta_vec, grad_vec, aData);
    return result;
}

template<typename T>
bool StageZeroGeneratorParsec<T>::IsLeapYear(const int &aYear) {
    if (aYear % 400 == 0) {
        return true;
    } else if (aYear % 100 == 0) {
        return false;
    } else if (aYear % 4 == 0) {
        return true;
    } else {
        return false;
    }
}

template<typename T>
double * StageZeroGeneratorParsec<T>::ReadObsFile(char *aFileName, const int &aNumLoc) {
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    double *z_vec = new double[aNumLoc];

    fp = fopen(aFileName, "r");
    if (fp == NULL) {
        throw std::runtime_error("readObsFile:cannot open observations file");
    }

    while ((read = getline(&line, &len, fp)) != -1 && count < aNumLoc) {
        z_vec[count++] = atof(line);
    }

    fclose(fp);
    if (line) free(line);
    return z_vec;
}

template<typename T> StageZeroGeneratorParsec<T> *StageZeroGeneratorParsec<T>::mpInstance = nullptr; 
