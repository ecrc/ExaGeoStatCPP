// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StageZeroGenerator.cpp
 * @brief Implementation of the StageZeroGenerator class
 * @version 1.1.0
**/

#include <data-generators/concrete/StageZeroGenerator.hpp>
#include <data-generators/LocationGenerator.hpp>
#if !DEFAULT_RUNTIME
#include <data-loader/concrete/ParsecLoader.hpp>
#else
#include <data-loader/concrete/CSVLoader.hpp>
#endif
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
#include <chameleon.h>
#include <chameleon/runtime.h>

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
StageZeroGenerator<T> *StageZeroGenerator<T>::GetInstance() {
    if (mpInstance == nullptr) {
        mpInstance = new StageZeroGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
StageZeroGenerator<T>::CreateData(Configurations &aConfigurations,
                                  exageostat::kernels::Kernel<T> &aKernel) {
    this->Runner(aConfigurations);
    return std::move(this->mData);
}

template<typename T>
void StageZeroGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

template<typename T>
void StageZeroGenerator<T>::Runner(Configurations &aConfigurations) {
    mArgs.mConfigs = &aConfigurations;
    this->ConfigureGenerator();
    this->Allocate();
    this->ReadForcingData();
    this->ReadNetCDFFiles();
    this->RunMeanTrend();
    
    this->CleanUp();
}

template<typename T>
void StageZeroGenerator<T>::ConfigureGenerator() {
    
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
    
    // Derive observation years and N from configuration (keep forcing length separate)
    int start_year = mArgs.mConfigs->GetStartYear();
    int end_year   = mArgs.mConfigs->GetEndYear();
    if (end_year < start_year) {
        throw std::runtime_error("EndYear must be >= StartYear");
    }
    int obs_years = (end_year - start_year + 1);
    mArgs.mN = static_cast<size_t>(mArgs.mT) * static_cast<size_t>(obs_years);
    
    // Number of locations (prefer explicit arg; fallback to configured ProblemSize)
    try { mArgs.mNumLocs = mArgs.mConfigs->GetNumLocs(); }
    catch (...) { mArgs.mNumLocs = mArgs.mConfigs->GetProblemSize(); }
    
    
}

template<typename T>
void StageZeroGenerator<T>::ReadNetCDFFiles() {

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
        double start_time = MPI_Wtime();
        if ((ret = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &id))) {
            double end_time = MPI_Wtime();
            fprintf(stderr, "[StageZero] FAILED to open NetCDF file: %s (%.3f sec) - Error: %s\n", 
                    filename, end_time - start_time, ncmpi_strerror(ret));
            throw std::runtime_error("Error opening NetCDF file: " + std::string(ncmpi_strerror(ret)));
        }
        double end_time = MPI_Wtime();
        fprintf(stderr, "[StageZero] SUCCESS opening NetCDF file: %s (%.3f sec)\n", filename, end_time - start_time);
        return id;
    };

    auto closeFileNCmpi = [&](int id) {
        int ret;
        if ((ret = ncmpi_close(id))) {
            throw std::runtime_error("Error closing NetCDF file: " + std::string(ncmpi_strerror(ret)));
        }
    };

    
    char path[256];
    std::string data_path = mArgs.mConfigs->GetDataPath();
    if (!data_path.empty() && data_path.back() != '/') data_path += '/';
    snprintf(path, sizeof(path), "%sdata_%d.nc", data_path.c_str(), start_year);

    fprintf(stderr, "[StageZero] Loading initial NetCDF file for dimension setup: %s\n", path);
    ncid = openFileNCmpi(path);

    // Dimension IDs and lengths with error checks (exactly like reference)
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
        snprintf(path2, sizeof(path2), "%sdata_%d.nc", data_path.c_str(), y);

        fprintf(stderr, "[StageZero] Loading NetCDF data for year %d: %s\n", y, path2);
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
void StageZeroGenerator<T>::ReadForcingData() {
    std::string forcing_path = mArgs.mConfigs->GetForcingDataPath();
    if (!forcing_path.empty() && forcing_path.back() == '/') forcing_path.pop_back();
    
    fprintf(stderr, "[StageZero] Loading forcing data from: %s\n", forcing_path.c_str());
    mArgs.mForcing = ReadObsFile((char *) forcing_path.c_str(), mArgs.mNoYears);
    fprintf(stderr, "[StageZero] Successfully loaded forcing data (%d years)\n\n", mArgs.mNoYears);
}

template<typename T>
void StageZeroGenerator<T>::Allocate() {

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
void StageZeroGenerator<T>::RunMeanTrend() {
    const int     lts      = mArgs.mConfigs->GetLowTileSize();
    const int     no_locs  = mArgs.mNumLocs;

#if defined(CHAMELEON_USE_MPI)
    if (CHAMELEON_Comm_rank() == 0) {
#endif
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
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    // Key bounds info
    fprintf(stderr, "Starting theta[0]: %f\n", mArgs.mStartingTheta[0]);
    fprintf(stderr, "Lower bound: %f, Upper bound: %f\n", mArgs.mLb[0], mArgs.mUp[0]);
    
    // Initialize CHAMELEON and create descriptors ONCE
    this->SetupMLEComponents();
    
    // Set up NLOPT optimization ONCE
    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, mArgs.mNumParams);
    nlopt_set_lower_bounds(opt, mArgs.mLb);
    nlopt_set_upper_bounds(opt, mArgs.mUp);
    nlopt_set_max_objective(opt, &StageZeroGenerator<T>::StageZeroObjectiveCallback, (void *)this);
    
    // Use configured tolerance and max iterations (defaults preserved in config)
    nlopt_set_ftol_abs(opt, mArgs.mConfigs->GetTolerance());
    nlopt_set_maxeval(opt, mArgs.mConfigs->GetMaxMleIterations());
    
    // Working theta vector (mutable for NLopt)
    std::vector<double> theta(mArgs.mNumParams, 0.0);
    
    // Process locations
    for (int l = 0; l < no_locs; ++l) {
        mArgs.mIterCount = 0;   // reset iteration counter
        mArgs.mCurrentLocation = l;  // Set current location for MLEAlgorithm

        // Fill Z for current location
        this->ConvertT2MToZForLocation(l);

        // Optimization
        fprintf(stderr, "[StageZero] Starting NLopt for location %d/%d\n", l + 1, no_locs);
        // Respect configured starting theta
        theta[0] = mArgs.mStartingTheta[0];
        double opt_f = 0.0;

        // Re-register the objective before each optimize call
        nlopt_set_max_objective(opt, &StageZeroGenerator<T>::StageZeroObjectiveCallback, (void *)this);
        nlopt_result nlres = nlopt_optimize(opt, theta.data(), &opt_f);
        fprintf(stderr, "[StageZero] NLopt finished (res=%d), theta=%0.10f, f=%0.10f\n\n", (int) nlres, theta[0], opt_f);

        // Recompute full pipeline once with optimal theta to generate CSV (apObj == nullptr)
        {
            std::vector<double> grad; // unused
            this->MLEAlgorithm(theta, grad, nullptr);
        }
    }
    
    // Cleanup
    nlopt_destroy(opt);
    
}

template<typename T>
void StageZeroGenerator<T>::CleanUp() {

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
    
    // Clean up CHAMELEON descriptors
    if (mArgs.mpDescZ) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpDescZ);
        mArgs.mpDescZ = nullptr;
    }
    if (mArgs.mpX) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpX);
        mArgs.mpX = nullptr;
    }
    if (mArgs.mpXtX) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpXtX);
        mArgs.mpXtX = nullptr;
    }
    if (mArgs.mpDescPart1) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpDescPart1);
        mArgs.mpDescPart1 = nullptr;
    }
    if (mArgs.mpDescPart2) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpDescPart2);
        mArgs.mpDescPart2 = nullptr;
    }
    if (mArgs.mpPart2Vector) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpPart2Vector);
        mArgs.mpPart2Vector = nullptr;
    }
    if (mArgs.mpEstimatedMeanTrend) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t**)&mArgs.mpEstimatedMeanTrend);
        mArgs.mpEstimatedMeanTrend = nullptr;
    }
}

template<typename T>
void StageZeroGenerator<T>::SetupMLEComponents() {
    
    int N = mArgs.mN;
    int nparams = 3 + 2*mArgs.mM;
    int dts = mArgs.mConfigs->GetDenseTileSize();
    int p_grid = 1, q_grid = 1; // Single process grid for now
    
    fprintf(stderr, "Creating descriptors: N=%d, nparams=%d, dts=%d\n", N, nparams, dts);
    
    // Initialize iteration counter
    mArgs.mIterCount = 0;
    
    // Allocate CHAMELEON descriptors - Use UNIFORM tile size for compatibility
    // Z vector (observations) - N x 1
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpDescZ, NULL, ChamRealDouble,
                          dts, dts, dts*dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    
    // X matrix (design matrix) - N x nparams
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpX, NULL, ChamRealDouble,
                          dts, dts, dts*dts, N, nparams, 0, 0, N, nparams, p_grid, q_grid);
    
    // XtX matrix - nparams x nparams (use same tile size as others)
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpXtX, NULL, ChamRealDouble,
                          dts, dts, dts*dts, nparams, nparams, 0, 0, nparams, nparams, p_grid, q_grid);
    
    // part1 scalar - 1 x 1 (use same tile size pattern)
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpDescPart1, NULL, ChamRealDouble,
                          dts, dts, dts*dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    
    // part2 scalar - 1 x 1 (use same tile size pattern)
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpDescPart2, NULL, ChamRealDouble,
                          dts, dts, dts*dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    
    // part2_vector - nparams x 1 (use same tile size pattern)
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpPart2Vector, NULL, ChamRealDouble,
                          dts, dts, dts*dts, nparams, 1, 0, 0, nparams, 1, p_grid, q_grid);
    
    // estimated_mean_trend vector - N x 1 (for final mean trend calculation)
    CHAMELEON_Desc_Create((CHAM_desc_t**)&mArgs.mpEstimatedMeanTrend, NULL, ChamRealDouble,
                          dts, dts, dts*dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    
    fprintf(stderr, "All CHAMELEON descriptors created successfully\n\n");
}

template<typename T>
void StageZeroGenerator<T>::ConvertT2MToZForLocation(int location_index) {
    
    auto *DESC_Z = static_cast<CHAM_desc_t *>(mArgs.mpDescZ);
    int N = mArgs.mN;

    if (location_index < 0 || location_index >= mArgs.mNumLocs) {
        throw std::runtime_error("ConvertT2MToZForLocation: invalid location index");
    }

    // Copy data from specified location to a LAPACK buffer, then tile it into Z
    std::vector<double> z_lap(static_cast<size_t>(N), 0.0);
    int count = std::min(N, mArgs.mT2mHourlyPerYearCount[location_index]);
    for (int i = 0; i < count; i++) {
        z_lap[static_cast<size_t>(i)] = mArgs.mT2mHourlyPerYear[location_index][i];
    }
    CHAMELEON_Lapack_to_Tile(static_cast<void*>(z_lap.data()), N, DESC_Z);
}

template<typename T>
void StageZeroGenerator<T>::GenerateDesignMatrixExact(double *matrix, int m, int n, int m0, int n0, double *localtheta) {
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
double StageZeroGenerator<T>::MLEAlgorithm(const std::vector<double> &aThetaVec,
                                           std::vector<double> &aGrad, void *apObj) {
    
    // Get descriptors
    auto *Zobs = static_cast<CHAM_desc_t *>(mArgs.mpDescZ);
    auto *X = static_cast<CHAM_desc_t *>(mArgs.mpX);
    auto *XtX = static_cast<CHAM_desc_t *>(mArgs.mpXtX);
    auto *part1 = static_cast<CHAM_desc_t *>(mArgs.mpDescPart1);
    auto *part2 = static_cast<CHAM_desc_t *>(mArgs.mpDescPart2);
    auto *part2_vector = static_cast<CHAM_desc_t *>(mArgs.mpPart2Vector);
    int N = X->m;
    double value = 0.0;
    const bool is_final_call = (apObj == nullptr); // final call after optimization
    
    // Null pointer checks
    if (!X) { throw std::runtime_error("X descriptor is null"); }
    if (!Zobs) { throw std::runtime_error("Zobs descriptor is null"); }
    
    
    
    // Compose local parameter vector
    double* localtheta = (double *) malloc((3+mArgs.mNoYears) * sizeof(double));
    localtheta[0] = aThetaVec[0];  // theta[0]
    localtheta[1] = mArgs.mT;      // Direct access: 8760
    localtheta[2] = mArgs.mM;      // Direct access: 10
    
    for(int ii=0; ii<mArgs.mNoYears; ii++)
        localtheta[3+ii] = mArgs.mForcing[ii];  // Direct access to forcing array
    
    try {
        // Step 1: generate design matrix X
        {
            int Nrows = X->m;
            int Ncols = X->n;
            double *x_lap = (double*)calloc((size_t)Nrows * (size_t)Ncols, sizeof(double));
            if (!x_lap) { throw std::runtime_error("Allocation failed for design matrix buffer"); }
            this->GenerateDesignMatrixExact(x_lap, Nrows, Ncols, 0, 0, localtheta);
            CHAMELEON_Lapack_to_Tile((void*)x_lap, Nrows, X);
            free(x_lap);
        }
        
        // Ensure workspaces are clean before computations
        CHAMELEON_dlaset_Tile(ChamUpperLower, 0.0, 0.0, part1);
        CHAMELEON_dlaset_Tile(ChamUpperLower, 0.0, 0.0, part2);
        CHAMELEON_dlaset_Tile(ChamUpperLower, 0.0, 0.0, part2_vector);
        CHAMELEON_dlaset_Tile(ChamUpperLower, 0.0, 0.0, XtX);

        // Step 2: part1 = Z^T * Z
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, Zobs, Zobs, 0.0, part1);
        
        // Step 3: part2_vector = X^T * Z
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, X, Zobs, 0.0, part2_vector);
    
        // Step 4: XtX = X^T * X
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, X, X, 0.0, XtX);

        // Step 5: Cholesky decomposition of XtX
        {
            int info = CHAMELEON_dpotrf_Tile(ChamLower, XtX);
            if(info != 0) {
                if (apObj) { free(localtheta); return -1e18; }
                else {
                    fprintf(stderr, "[StageZero] ERROR: Cholesky(X^T X) failed in final pass. Aborting.\n");
                    free(localtheta);
                    exit(1);
                }
            }
        }
    
        // Step 6: First triangular solve: y = L^{-1} (X^T Z)
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1.0, XtX, part2_vector);
        // For the final pass compute beta = L^{-T} y before forming the trend
        if (!apObj) {
            CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1.0, XtX, part2_vector);
        }

        // If this is an optimization callback, compute and return the C objective using y^T y
        if (apObj) {
            double part1_val_cb = 0.0;
            if (part1 && part1->mat) part1_val_cb = static_cast<double *>(part1->mat)[0];
            // Compute part2 = y^T y
            CHAMELEON_dlaset_Tile(ChamUpperLower, 0.0, 0.0, part2);
            CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, part2_vector, part2_vector, 0.0, part2);
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

        // Optional: second triangular solve is only needed for final beta; handled below

        // Step 7: Final computation with CSV output
        // Create estimated_mean_trend = X * beta (using part2_vector as beta)
        auto *estimated_mean_trend = static_cast<CHAM_desc_t *>(mArgs.mpEstimatedMeanTrend);
        CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1.0, X, part2_vector, 0.0, estimated_mean_trend);
        
        // Residuals = Z - trend (stored back into estimated_mean_trend)
        CHAMELEON_dgeadd_Tile(ChamNoTrans, 1.0, Zobs, -1.0, estimated_mean_trend);
        
        // Compute sigma² = residuals^T * residuals
        auto *sigma_desc = static_cast<CHAM_desc_t *>(mArgs.mpDescPart2);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1.0, estimated_mean_trend, estimated_mean_trend, 0.0, sigma_desc);
        
        // Get sigma² value and divide by N (like C version)
        double sigma_squared_raw = static_cast<double *>(sigma_desc->mat)[0];
        double sigma_squared = sigma_squared_raw / N;
        
        // Pull residuals to LAPACK layout for correct order and normalization
        double *emt_lap = (double*)calloc((size_t)N, sizeof(double));
        if (!emt_lap) { throw std::runtime_error("Allocation failed for emt_lap"); }
        CHAMELEON_Tile_to_Lapack(estimated_mean_trend, emt_lap, N);
        
        // Standardize residuals: residuals = residuals / sqrt(sigma²)
        double sqrt_sigma = std::sqrt(sigma_squared);
        for (int i = 0; i < N; ++i) {
            emt_lap[i] /= sqrt_sigma;
        }
        
        // Set objective value (maximize -sigma^2 -> minimize sigma^2)
        value = -sigma_squared;

        // Static global arrays for multi-location storage (like C version)
        static std::vector<std::vector<double>> Z_new(mArgs.mNumLocs, std::vector<double>(N));
        static std::vector<std::vector<double>> params(mArgs.mNumLocs, std::vector<double>(3 + 2*mArgs.mM + 2));
        
        // Store normalized residuals for this location
        int location_index = mArgs.mCurrentLocation;  // Current location being processed (0-based)
        for (int i = 0; i < N; ++i) {
            Z_new[location_index][i] = emt_lap[i];
        }
        free(emt_lap);
        
        // Pull beta to LAPACK layout and store parameters
        int vlen = part2_vector->m;
        double *p2v_lap = (double*)calloc((size_t)vlen, sizeof(double));
        if (!p2v_lap) { throw std::runtime_error("Allocation failed for p2v_lap"); }
        CHAMELEON_Tile_to_Lapack(part2_vector, p2v_lap, vlen);
        params[location_index][0] = aThetaVec[0];  // optimized theta
        params[location_index][1] = sigma_squared;  // sigma²
        for(int i = 0; i < 3 + 2*mArgs.mM; i++) {
            params[location_index][i+2] = p2v_lap[i];
        }
        free(p2v_lap);
        
        // Print debug info matching C version exactly
        // remove noisy debug prints
        
        // Write CSV files only on the final call and when processing the last location
        if (is_final_call && (location_index == mArgs.mNumLocs - 1)) {
            std::string results_path;
            try { 
                results_path = mArgs.mConfigs->GetResultsPath(); 
            } catch (...) { 
                results_path.clear(); 
            }
            
            if (results_path.empty()) { 
                fprintf(stderr, "[StageZero] ERROR: ResultsPath is required. Please set --resultspath.\n");
                throw std::runtime_error("ResultsPath is required. Please set --resultspath.");
            }
            
            if (results_path.back() != '/') results_path += "/";
            
            // Create results directory if it doesn't exist
            try {
                if (!std::filesystem::exists(results_path)) {
                    fprintf(stderr, "[StageZero] Creating results directory: %s\n", results_path.c_str());
                    std::filesystem::create_directories(results_path);
                }
            } catch (const std::exception &e) {
                fprintf(stderr, "[StageZero] ERROR: Failed to create output directory: %s (%s)\n", results_path.c_str(), e.what());
                throw std::runtime_error("Failed to create output directory: " + std::string(e.what()));
            }
            
            fprintf(stderr, "[StageZero] Writing outputs to: %s\n", results_path.c_str());
            
            // Write Z files for each time slot (replace mode instead of append)
            #pragma omp parallel for
            for (int time_slot = 0; time_slot < N; time_slot++) {
                char file_path[256];
                std::snprintf(file_path, sizeof(file_path), "%sz_%d.csv", results_path.c_str(), time_slot);
                
                // Retry loop for file locking with replace mode
                while (true) {
                    FILE *fp = std::fopen(file_path, "w");  // Replace mode instead of append
                    if (!fp) {
                        fprintf(stderr, "[StageZero] WARNING: Failed to open file %s, retrying...\n", file_path);
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
                fprintf(stderr, "[StageZero] ERROR: Failed to create params file: %s\n", params_file_path);
                exit(1); 
            }
            
            for(int k = 0; k < mArgs.mNumLocs; k++) {
                for(int i = 0; i < 3 + 2*mArgs.mM + 2; i++) {
                    fprintf(fp, "%.14f ", params[k][i]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            
            fprintf(stderr, "[StageZero] Successfully wrote %d Z files and params.csv to %s\n", N, results_path.c_str());
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
double StageZeroGenerator<T>::StageZeroObjectiveCallback(unsigned aN, const double *aTheta, double *aGrad, void *aData) {
    auto *generator = static_cast<StageZeroGenerator<T> *>(aData);
    std::vector<double> theta_vec(aTheta, aTheta + aN);
    std::vector<double> grad_vec(aN);
    double result = generator->MLEAlgorithm(theta_vec, grad_vec, aData);
    return result;
}

template<typename T>
bool StageZeroGenerator<T>::IsLeapYear(const int &aYear) {
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
double * StageZeroGenerator<T>::ReadObsFile(char *aFileName, const int &aNumLoc) {
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    double *z_vec = new double[aNumLoc];

    double start_time = MPI_Wtime();
    fp = fopen(aFileName, "r");
    if (fp == NULL) {
        double end_time = MPI_Wtime();
        fprintf(stderr, "[StageZero] FAILED to open file: %s (%.3f sec)\n", aFileName, end_time - start_time);
        throw std::runtime_error("readObsFile:cannot open observations file: " + std::string(aFileName));
    }
    double end_time = MPI_Wtime();
    fprintf(stderr, "[StageZero] SUCCESS opening file: %s (%.3f sec)\n", aFileName, end_time - start_time);

    while ((read = getline(&line, &len, fp)) != -1 && count < aNumLoc) {
        z_vec[count++] = atof(line);
    }

    fclose(fp);
    if (line) free(line);
    return z_vec;
}

template<typename T> StageZeroGenerator<T> *StageZeroGenerator<T>::mpInstance = nullptr;

// Explicit instantiation of StageZeroGenerator for supported types
template class StageZeroGenerator<float>;
template class StageZeroGenerator<double>;
