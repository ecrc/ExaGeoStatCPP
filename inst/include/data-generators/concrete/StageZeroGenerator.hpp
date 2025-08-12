#ifndef EXAGEOSTAT_STAGEZEROGENERATOR_HPP
#define EXAGEOSTAT_STAGEZEROGENERATOR_HPP

#include "../DataGenerator.hpp"
#include <configurations/Configurations.hpp>
#include <cstddef>
#include <memory>
#include <vector>
#include <string>

// Forward declarations to avoid heavy includes
namespace exageostat { namespace kernels { template<typename T> class Kernel; } }

// Fallback define for template instantiation macro
#ifndef EXAGEOSTAT_INSTANTIATE_CLASS
#define EXAGEOSTAT_INSTANTIATE_CLASS(...)
#endif

namespace exageostat::generators::stagezero {

    /**
     * @brief Holds runtime configuration and state for Stage Zero.
     */
    struct StageZeroArgs {
        // Model/config
        int mM = 10;                        // harmonics for mean-trend
        int mT = 365 * 24;                  // period in hours (used to compute N)
        int mNoYears = 751;                 // forcing length (does not change obs N)
        int mNumParams = 1;                 // number of optimized parameters (theta)
        size_t mN = 365 * 24 * 3;           // observations count default (T * years)
        int mNumLocs = 2;                   // number of locations

        // Optimization bounds and vectors
        double *mStartingTheta = nullptr;
        double *mTargetTheta = nullptr;
        double *mInitialTheta = nullptr;
        double *mLb = nullptr;
        double *mUp = nullptr;

        // Forcing and input buffers
        double *mForcing = nullptr;
        double **mT2mHourlyPerYear = nullptr;
        int *mT2mHourlyPerYearCount = nullptr;

        // Runtime grid and flags
        int mPGrid = 1;
        int mQGrid = 1;
        int mZMiss = 0;
        int mLog = 0;
        exageostat::configurations::Configurations *mConfigs = nullptr;
        int mNcid = 0;
        int mAsync = 0;
        int mDiagThick = 1;
        int mCheck = 0;
        int mHicmaMaxRank = 0;

        // CHAMELEON descriptors
        void *mpDescZ = nullptr;            // Z observations vector
        void *mpX = nullptr;                // Design matrix X
        void *mpXtX = nullptr;              // X^T * X matrix
        void *mpDescPart1 = nullptr;        // part1 scalar
        void *mpDescPart2 = nullptr;        // part2 scalar
        void *mpPart2Vector = nullptr;      // part2_vector
        void *mpEstimatedMeanTrend = nullptr; // estimated mean trend

        // Scalars/counters
        double mPart1 = 0.0;
        double mPart2 = 0.0;
        int mIterCount = 0;
        int mCurrentLocation = 0;
    };

    /**
     * @class StageZeroGenerator
     * @brief Stage Zero pipeline: read inputs, build mean-trend, optimize, and write CSV outputs.
     * @tparam T float or double
     */
    template<typename T>
    class StageZeroGenerator : public DataGenerator<T> {

    public:

        /**
         * @brief Get a pointer to the singleton instance.
         * @return Pointer to the `StageZeroGenerator<T>` instance.
         */
        static StageZeroGenerator<T> *GetInstance();

        /**
         * @brief Configure and run Stage Zero, returning the generated data container.
         * @param aConfigurations Global configuration object.
         * @param aKernel Kernel (unused by Stage Zero mean-trend; retained for API compatibility).
         * @return Unique pointer to `ExaGeoStatData<T>`.
         */
        std::unique_ptr<ExaGeoStatData<T>>
        CreateData(configurations::Configurations &aConfigurations,
                   exageostat::kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Release the singleton instance.
         */
        static void ReleaseInstance();

        /**
         * @brief Core MLE objective (apObj != nullptr) and final compute+write path (apObj == nullptr).
         * @param aThetaVec Current theta vector.
         * @param aGrad Gradient vector (unused by BOBYQA).
         * @param apObj When non-null, acts as NLopt objective; when null, runs final compute and writes outputs.
         * @return Objective value (negative log-likelihood component in objective mode; -sigma^2 in final path).
         */
        double MLEAlgorithm(const std::vector<double> &aThetaVec,
                            std::vector<double> &aGrad, void *apObj);

    private:
        /**
         * @brief NLopt objective callback.
         * @param aN Number of parameters.
         * @param aTheta Parameter vector.
         * @param aGrad Gradient vector.
         * @param aData Pointer to StageZeroGenerator instance.
         * @return Objective value.
         */
        static double StageZeroObjectiveCallback(unsigned aN, const double *aTheta, double *aGrad, void *aData);

        /**
         * @brief Orchestrate Stage Zero pipeline (configure, allocate, load, optimize, cleanup).
         * @param aConfigurations Global configuration object.
         */
        void Runner(configurations::Configurations &aConfigurations);

        /**
         * @brief Initialize internal arguments from configuration (years, N, locations, bounds...).
         */
        void ConfigureGenerator();

        /**
         * @brief Read NetCDF observations for the configured year range into per-location buffers.
         */
        void ReadNetCDFFiles();

        /**
         * @brief Load forcing series from file (length `mNoYears`).
         */
        void ReadForcingData();

        /**
         * @brief Allocate arrays and copy configuration vectors (bounds, theta, forcing, buffers).
         */
        void Allocate();

        /**
         * @brief Execute mean-trend optimization per location and write CSV outputs.
         */
        void RunMeanTrend();

        /**
         * @brief Leap year helper.
         * @param aYear Year integer.
         * @return true if leap year, else false.
         */
        bool IsLeapYear(const int &aYear);

        /**
         * @brief Read a text file of doubles (one per line) into a heap array.
         * @param aFileName Path to input file.
         * @param aNumLoc Number of values to read (capacity).
         * @return Pointer to newly allocated array (caller owns).
         */
        double * ReadObsFile(char *aFileName, const int &aNumLoc);

        /**
         * @brief Free allocated arrays and CHAMELEON descriptors.
         */
        void CleanUp();

        /**
         * @brief Create CHAMELEON descriptors for Z, X, XtX, part scalars and buffers.
         */
        void SetupMLEComponents();
        
        /**
         * @brief Copy all locations to Z (batch path; per-location path preferred).
         */
        void ConvertT2MToZ();

        /**
         * @brief Copy a single location's observations into the Z descriptor.
         * @param location_index Index of the location to copy.
         */
        void ConvertT2MToZForLocation(int location_index);
        
        /**
         * @brief Generate design matrix X (column-major) for the mean-trend model.
         * @param matrix LAPACK-style buffer to write (size m*n).
         * @param m Number of rows (observations).
         * @param n Number of columns (parameters).
         * @param m0 Row offset.
         * @param n0 Column offset.
         * @param localtheta [theta, T, M, forcing...].
         */
        void GenerateDesignMatrixExact(double *matrix, int m, int n, int m0, int n0, double *localtheta);
        
        /**
         * @brief Constructor for the SyntheticGenerator class.
         * @return void
         */
        StageZeroGenerator() = default;

        /**
         * @brief Default destructor.
         */
        ~StageZeroGenerator() override = default;

        // Pointer to the singleton instance
        static StageZeroGenerator<T> *mpInstance;

    private:

        std::unique_ptr<ExaGeoStatData<T>> mData;
        StageZeroArgs mArgs;

    };

} // namespace exageostat

#endif //EXAGEOSTAT_STAGEZEROGENERATOR_HPP
