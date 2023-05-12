
/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <kernels/Kernel.hpp>
#include <common/PluginRegistry.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace exageostat::linearAlgebra;
using namespace exageostat::generators::Synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;

SyntheticGenerator::SyntheticGenerator(SyntheticDataConfigurations *apConfigurations) {
    // Set configuration map and init locations.
    this->mpConfigurations = apConfigurations;
    this->mpLocations = new Locations();

    // Set selected Kernel
    std::string kernel_name = this->mpConfigurations->GetKernel();
    this->mpKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
            this->mpConfigurations->GetKernel());
    this->mpKernel->SetPValue(this->mpConfigurations->GetTimeSlot());

    int parameters_number = this->mpKernel->GetParametersNumbers();
    this->mpConfigurations->SetParametersNumber(parameters_number);

    // Check if any theta is not initialized, This means that the user didn't send it as an argument
    if(this->mpConfigurations->GetLowerBounds() == nullptr){
        this->mpConfigurations->SetLowerBounds(InitTheta(this->mpConfigurations->GetLowerBounds(), parameters_number));
    }
    if(this->mpConfigurations->GetUpperBounds() == nullptr){
        this->mpConfigurations->SetUpperBounds(InitTheta(this->mpConfigurations->GetUpperBounds(), parameters_number));
    }
    if(this->mpConfigurations->GetInitialTheta() == nullptr){
        this->mpConfigurations->SetInitialTheta(InitTheta(this->mpConfigurations->GetInitialTheta(), parameters_number));
    }
    if(this->mpConfigurations->GetTargetTheta() == nullptr){
        this->mpConfigurations->SetTargetTheta(InitTheta(this->mpConfigurations->GetTargetTheta(), parameters_number));
    }

    // Set starting theta with the lower bounds values
    this->mpConfigurations->SetStartingTheta(this->mpConfigurations->GetLowerBounds());

    for (int i = 0; i < parameters_number; i++) {
        if (this->mpConfigurations->GetTargetTheta()[i] != -1) {
            this->mpConfigurations->GetLowerBounds()[i] = this->mpConfigurations->GetTargetTheta()[i];
            this->mpConfigurations->GetUpperBounds()[i] = this->mpConfigurations->GetTargetTheta()[i];
            this->mpConfigurations->GetStartingTheta()[i] = this->mpConfigurations->GetTargetTheta()[i];
        }
    }
}

void SyntheticGenerator::GenerateDescriptors() {

    // Create and initialize linear algebra solvers for different precision types.
    if (this->mpConfigurations->GetPrecision() == SINGLE) {
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(
                this->mpConfigurations->GetComputation());
        linearAlgebraSolver->SetConfigurations(this->mpConfigurations);
        linearAlgebraSolver->InitiateDescriptors();
    } else if (this->mpConfigurations->GetPrecision() == DOUBLE) {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(
                this->mpConfigurations->GetComputation());
        linearAlgebraSolver->SetConfigurations(this->mpConfigurations);
        linearAlgebraSolver->InitiateDescriptors();
    } else if (this->mpConfigurations->GetPrecision() == MIXED) {
        // TODO: Add implementation for mixed-precision linear algebra solver.
    }

}

void SyntheticGenerator::GenerateLocations() {

    int p;
    if (this->mpKernel) {
        p = this->mpKernel->GetPValue();
    } else {
        throw std::runtime_error("Error in Allocating Kernel plugin");
    }
    int N = this->mpConfigurations->GetProblemSize();

    int index = 0;
    Dimension dimension = this->mpConfigurations->GetDimension();
    int time_slots = this->mpConfigurations->GetTimeSlot();

    //Allocate memory
    this->mpLocations->SetLocationX((double *) malloc(N * time_slots * sizeof(double)));
    this->mpLocations->SetLocationY((double *) malloc(N * time_slots * sizeof(double)));

    if (dimension != Dimension2D) {
        this->mpLocations->SetLocationZ((double *) malloc(N * time_slots * sizeof(double)));
    }

    int rootN;
    if (dimension == Dimension3D) {
        //Cubic root.
        rootN = ceil(cbrt(N));
    } else {
        //Square root.
        rootN = ceil(sqrt(N));
    }

    int *grid = (int *) calloc((int) rootN, sizeof(int));
    for (auto i = 0; i < rootN; i++) {
        grid[i] = i + 1;
    }

    for (auto i = 0; i < rootN && index < N; i++) {
        for (auto j = 0; j < rootN && index < N; j++) {
            if (dimension == Dimension3D) {
                for (auto k = 0; k < rootN && index < N; k++) {
                    this->mpLocations->GetLocationX()[index] = (grid[i] - 0.5 + UniformDistribution(-0.4, 0.4)) / rootN;
                    this->mpLocations->GetLocationY()[index] = (grid[j] - 0.5 + UniformDistribution(-0.4, 0.4)) / rootN;
                    this->mpLocations->GetLocationZ()[index] = (grid[k] - 0.5 + UniformDistribution(-0.4, 0.4)) / rootN;
                    index++;
                }
            } else {
                this->mpLocations->GetLocationX()[index] = (grid[i] - 0.5 + UniformDistribution(-0.4, 0.4)) / rootN;
                this->mpLocations->GetLocationY()[index] = (grid[j] - 0.5 + UniformDistribution(-0.4, 0.4)) / rootN;
                if (dimension == DimensionST) {
                    this->mpLocations->GetLocationZ()[index] = 1.0;
                }
                index++;
            }
        }
    }
    free(grid);
    if (dimension != DimensionST) {
        SortLocations(N);
    } else {
        for (auto j = 1; j < time_slots; j++) {
            for (auto i = 0; i < N; i++) {
                this->mpLocations->GetLocationX()[i + j * N] = this->mpLocations->GetLocationX()[i];
                this->mpLocations->GetLocationY()[i + j * N] = this->mpLocations->GetLocationY()[i];
                this->mpLocations->GetLocationZ()[i + j * N] = (double) (j + 1);
            }
        }
    }
}

double SyntheticGenerator::UniformDistribution(double aRangeLow, double aRangeHigh) {

    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = aRangeHigh - aRangeLow;
    return (myRand * range) + aRangeLow;
}


void SyntheticGenerator::SortLocations(int aN) {

    // Some sorting, required by spatial statistics code
    uint16_t x, y, z;
    Dimension dimension = this->mpConfigurations->GetDimension();
    uint64_t vectorZ[aN];

    // Encode data into vector z
    for (auto i = 0; i < aN; i++) {
        x = (uint16_t) (this->mpLocations->GetLocationX()[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (this->mpLocations->GetLocationY()[i] * (double) UINT16_MAX + .5);
        if (dimension != Dimension2D) {
            z = (uint16_t) (this->mpLocations->GetLocationZ()[i] * (double) UINT16_MAX + .5);
        } else {
            z = (uint16_t) 0.0;
        }
        vectorZ[i] = (SpreadBits(z) << 2) + (SpreadBits(y) << 1) + SpreadBits(x);
    }
    // Sort vector z
    std::sort(vectorZ, vectorZ + aN, CompareUint64);

    // Decode data from vector z
    for (auto i = 0; i < aN; i++) {
        x = ReverseSpreadBits(vectorZ[i] >> 0);
        y = ReverseSpreadBits(vectorZ[i] >> 1);
        z = ReverseSpreadBits(vectorZ[i] >> 2);
        this->mpLocations->GetLocationX()[i] = (double) x / (double) UINT16_MAX;
        this->mpLocations->GetLocationY()[i] = (double) y / (double) UINT16_MAX;
        if (dimension == Dimension3D) {
            this->mpLocations->GetLocationZ()[i] = (double) z / (double) UINT16_MAX;
        }
    }
}

uint64_t SyntheticGenerator::SpreadBits(uint64_t aInputByte) {
    aInputByte &= 0x000000000000ffff;
    // aInputByte = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    aInputByte = (aInputByte ^ (aInputByte << 24)) & 0x000000ff000000ff;
    // aInputByte = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    aInputByte = (aInputByte ^ (aInputByte << 12)) & 0x000f000f000f000f; //000 7000f000f000f
    // aInputByte = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    aInputByte = (aInputByte ^ (aInputByte << 6)) & 0x0303030303030303; //0 0001 0 0011 0 0011 0 0011 0
    // aInputByte = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    aInputByte = (aInputByte ^ (aInputByte << 3)) & 0x1111111111111111;
    // aInputByte = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    return aInputByte;
}

uint64_t SyntheticGenerator::ReverseSpreadBits(uint64_t aInputByte) {

    aInputByte &= 0x1111111111111111;
    // aInputByte = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    aInputByte = (aInputByte ^ (aInputByte >> 3)) & 0x0303030303030303;
    // aInputByte = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    aInputByte = (aInputByte ^ (aInputByte >> 6)) & 0x000f000f000f000f;
    // aInputByte = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    aInputByte = (aInputByte ^ (aInputByte >> 12)) & 0x000000ff000000ff;
    // aInputByte = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    aInputByte = (aInputByte ^ (aInputByte >> 24)) & 0x000000000000ffff;
    // aInputByte = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    return aInputByte;

}

bool SyntheticGenerator::CompareUint64(const uint64_t &aFirstValue, const uint64_t &aSecondValue) {
    return aFirstValue < aSecondValue;
}

void SyntheticGenerator::GenerateObservations() {

    const auto& configurations = this->mpConfigurations;
    auto descriptorC = this->mpConfigurations->GetDescriptorC()[0];
    const auto& l1 = this->GetLocations();

    std::unique_ptr<double[]> initial_theta(new double[this->mpKernel->GetParametersNumbers()]);
    initial_theta[0] = 1.0;
    initial_theta[1] = 0.1;
    initial_theta[2] = 0.5;

    switch (configurations->GetPrecision()) {
        case SINGLE: {
            auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(configurations->GetComputation());
            linearAlgebraSolver->SetConfigurations(configurations);
            auto* A = (double* ) linearAlgebraSolver->EXAGEOSTAT_DATA_GET_ADDRESS(descriptorC, 0, 0);
            mpKernel->GenerateCovarianceMatrix(A, 5, 5, 0, 0, l1, l1, nullptr, initial_theta.get(), 0);
            break;
        }
        case DOUBLE: {
            auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(configurations->GetComputation());
            linearAlgebraSolver->SetConfigurations(configurations);
            auto* A = (double* ) linearAlgebraSolver->EXAGEOSTAT_DATA_GET_ADDRESS(descriptorC, 0, 0);
            mpKernel->GenerateCovarianceMatrix(A, 5, 5, 0, 0, l1, l1, nullptr, initial_theta.get(), 0);
            break;
        }
        case MIXED: {
            // TODO: Add implementation for mixed-precision linear algebra solver.
            break;
        }
    }

}

double* SyntheticGenerator::InitTheta(double *apTheta, int size) {
    apTheta = (double*) malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        apTheta[i] = -1;
    }
    return apTheta;
}
