#ifndef EXAGEOSTAT_STAGEZEROGENERATORPARSEC_HPP
#define EXAGEOSTAT_STAGEZEROGENERATORPARSEC_HPP

#include "../DataGenerator.hpp"
#include <configurations/Configurations.hpp>
#include <cstddef>
#include <memory>
#include <vector>
#include <string>

// Forward declarations to avoid heavy includes
namespace exageostat { 
    namespace kernels { template<typename T> class Kernel; }
    namespace configurations { class Configurations; }
}

// Fallback define for template instantiation macro
#ifndef EXAGEOSTAT_INSTANTIATE_CLASS
#define EXAGEOSTAT_INSTANTIATE_CLASS(...)
#endif

// PaRSEC includes
extern "C" {
#include <parsec.h>
#include <parsec/data_dist/matrix/matrix.h>
#include <dplasma.h>
#include <dplasma/types.h>
}

namespace exageostat::generators::stagezero {

    /**
     * @brief Holds runtime configuration and state for Stage Zero with PaRSEC.
     */
    struct StageZeroArgsParsec {
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

        // PaRSEC descriptors (replacing CHAMELEON descriptors)
        parsec_matrix_block_cyclic_t *mpDescZ = nullptr;            // Z observations vector
        parsec_matrix_block_cyclic_t *mpX = nullptr;                // X design matrix
        parsec_matrix_block_cyclic_t *mpXtX = nullptr;              // X^T * X matrix
        parsec_matrix_block_cyclic_t *mpDescPart1 = nullptr;        // part1 scalar
        parsec_matrix_block_cyclic_t *mpDescPart2 = nullptr;        // part2 scalar
        parsec_matrix_block_cyclic_t *mpPart2Vector = nullptr;      // part2 vector
        parsec_matrix_block_cyclic_t *mpEstimatedMeanTrend = nullptr; // estimated mean trend

        // PaRSEC context
        parsec_context_t *mpParsecContext = nullptr;

        // Iteration tracking
        int mIterCount = 0;
        int mCurrentLocation = 0;
    };

    /**
     * @brief Stage Zero data generator using PaRSEC/DPLASMA for climate data preprocessing.
     * @details This class implements mean trend removal from time series climate data
     * using PaRSEC runtime and DPLASMA linear algebra operations.
     * @tparam T Data Type: float or double
     */
    template<typename T>
    class StageZeroGeneratorParsec : public DataGenerator<T> {
    public:
        /**
         * @brief Creates data using the Stage Zero mean trend removal pipeline.
         * @param[in] aConfigurations Reference to Configurations object.
         * @param[in] aKernel Reference to Kernel object.
         * @return Unique pointer to ExaGeoStatData object.
         */
        std::unique_ptr<ExaGeoStatData<T>> CreateData(exageostat::configurations::Configurations &aConfigurations,
                                                      exageostat::kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Gets the singleton instance of StageZeroGeneratorParsec.
         * @return Pointer to StageZeroGeneratorParsec instance.
         */
        static StageZeroGeneratorParsec<T> *GetInstance();

        /**
         * @brief Releases the singleton instance.
         */
        static void ReleaseInstance();

    private:
        /// Singleton instance
        static StageZeroGeneratorParsec<T> *mpInstance;

        /// Runtime arguments and state
        StageZeroArgsParsec mArgs;
        
        /// Data object (for compatibility with base class)
        std::unique_ptr<ExaGeoStatData<T>> mData;

        /**
         * @brief Main execution pipeline for Stage Zero.
         * @param[in] aConfigurations Reference to Configurations object.
         */
        void Runner(exageostat::configurations::Configurations &aConfigurations);

        /**
         * @brief Configures the generator with parameters.
         */
        void ConfigureGenerator();

        /**
         * @brief Allocates memory for arrays and descriptors.
         */
        void Allocate();

        /**
         * @brief Reads NetCDF files for climate data.
         */
        void ReadNetCDFFiles();

        /**
         * @brief Reads forcing data from CSV file.
         */
        void ReadForcingData();

        /**
         * @brief Runs the mean trend removal algorithm using PaRSEC/DPLASMA.
         */
        void RunMeanTrend();

        /**
         * @brief Sets up PaRSEC components and descriptors.
         */
        void SetupMLEComponents();

        /**
         * @brief Converts T2M data to Z vector for a specific location.
         * @param[in] location_index Index of the location to process.
         */
        void ConvertT2MToZForLocation(int location_index);

        /**
         * @brief Generates the design matrix for mean trend modeling.
         * @param[in] matrix Pointer to matrix data.
         * @param[in] m Number of rows.
         * @param[in] n Number of columns.
         * @param[in] m0 Starting row index.
         * @param[in] n0 Starting column index.
         * @param[in] localtheta Local theta parameters.
         */
        void GenerateDesignMatrixExact(double *matrix, int m, int n, int m0, int n0, double *localtheta);

        /**
         * @brief MLE algorithm implementation using PaRSEC/DPLASMA.
         * @param[in] aThetaVec Vector of theta parameters.
         * @param[out] aGrad Gradient vector (unused in this implementation).
         * @param[in] apObj Pointer to objective function data.
         * @return Objective function value.
         */
        double MLEAlgorithm(const std::vector<double> &aThetaVec,
                           std::vector<double> &aGrad, void *apObj);

        /**
         * @brief Objective function callback for NLOPT optimization.
         * @param[in] aN Number of parameters.
         * @param[in] aTheta Parameter vector.
         * @param[out] aGrad Gradient vector.
         * @param[in] aData Pointer to generator instance.
         * @return Objective function value.
         */
        static double StageZeroObjectiveCallback(unsigned aN, const double *aTheta, double *aGrad, void *aData);

        /**
         * @brief Checks if a year is a leap year.
         * @param[in] aYear Year to check.
         * @return True if leap year, false otherwise.
         */
        bool IsLeapYear(const int &aYear);

        /**
         * @brief Reads observation data from file.
         * @param[in] aFileName File name to read from.
         * @param[in] aNumLoc Number of locations.
         * @return Pointer to data array.
         */
        double *ReadObsFile(char *aFileName, const int &aNumLoc);

        /**
         * @brief Cleans up allocated memory and descriptors.
         */
        void CleanUp();
    };

    /**
     * @brief Instantiates the StageZeroGeneratorParsec class for float and double types.
     * @tparam T Data Type: float or double
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(StageZeroGeneratorParsec)

}//namespace exageostat

#endif //EXAGEOSTAT_STAGEZEROGENERATORPARSEC_HPP 
