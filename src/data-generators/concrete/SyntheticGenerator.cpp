
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

template<typename T>
SyntheticGenerator<T>::SyntheticGenerator(SyntheticDataConfigurations *apConfigurations) {

    // Set configuration map and init locations.
    this->mpConfigurations = apConfigurations;
    this->mpLocations = new Locations();

    // Set selected Kernel
    std::string kernel_name = this->mpConfigurations->GetKernel();
    this->mpKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
            this->mpConfigurations->GetKernel());
    this->mpKernel->SetPValue(this->mpConfigurations->GetTimeSlot());

    this->mpConfigurations->SetProblemSize(this->mpConfigurations->GetProblemSize() * this->mpKernel->GetPValue());

    int parameters_number = this->mpKernel->GetParametersNumbers();
    this->mpConfigurations->SetParametersNumber(parameters_number);

    this->mpConfigurations->SetLowerBounds(InitTheta(this->mpConfigurations->GetLowerBounds(), parameters_number));
    this->mpConfigurations->SetUpperBounds(InitTheta(this->mpConfigurations->GetUpperBounds(), parameters_number));
    this->mpConfigurations->SetInitialTheta(InitTheta(this->mpConfigurations->GetInitialTheta(), parameters_number));
    this->mpConfigurations->SetTargetTheta(InitTheta(this->mpConfigurations->GetTargetTheta(), parameters_number));

    // Set starting theta with the lower bounds values
    this->mpConfigurations->SetStartingTheta(this->mpConfigurations->GetLowerBounds());

    for (int i = 0; i < parameters_number; i++) {
        if (this->mpConfigurations->GetTargetTheta()[i] != -1) {
            this->mpConfigurations->GetLowerBounds()[i] = this->mpConfigurations->GetTargetTheta()[i];
            this->mpConfigurations->GetUpperBounds()[i] = this->mpConfigurations->GetTargetTheta()[i];
            this->mpConfigurations->GetStartingTheta()[i] = this->mpConfigurations->GetTargetTheta()[i];
        }
    }

    // Set linear Algebra solver
    this->mpLinearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            this->mpConfigurations->GetComputation());
    this->mpLinearAlgebraSolver->SetConfigurations(this->mpConfigurations);

}

template<typename T>
void SyntheticGenerator<T>::GenerateDescriptors() {

    this->mpLinearAlgebraSolver->InitiateDescriptors();
}

template<typename T>
void SyntheticGenerator<T>::GenerateLocations() {

    int p;
    if (this->mpKernel) {
        p = this->mpKernel->GetPValue();
    } else {
        throw std::runtime_error("Error in Allocating Kernel plugin");
    }
    int N = this->mpConfigurations->GetProblemSize() / p;

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

template<typename T>
double SyntheticGenerator<T>::UniformDistribution(double aRangeLow, double aRangeHigh) {

    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = aRangeHigh - aRangeLow;
    return (myRand * range) + aRangeLow;
}

template<typename T>
void SyntheticGenerator<T>::SortLocations(int aN) {

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

template<typename T>
uint64_t SyntheticGenerator<T>::SpreadBits(uint64_t aInputByte) {
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

template<typename T>
uint64_t SyntheticGenerator<T>::ReverseSpreadBits(uint64_t aInputByte) {

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

template<typename T>
bool SyntheticGenerator<T>::CompareUint64(const uint64_t &aFirstValue, const uint64_t &aSecondValue) {
    return aFirstValue < aSecondValue;
}

template<typename T>
void SyntheticGenerator<T>::GenerateObservations() {

    auto descriptorC = this->mpConfigurations->GetDescriptorC()[0];
    const auto &l1 = this->GetLocations();

    this->mpLinearAlgebraSolver->GenerateObservationsVector(descriptorC, l1, l1, nullptr,
                                                            this->mpConfigurations->GetInitialTheta(), 0,
                                                            this->mpKernel);

}

template<typename T>
std::vector<double> SyntheticGenerator<T>::InitTheta(std::vector<double> apTheta, int size) {

    // If null, this mean user have not passed the values arguments, Make values equal -1
    if (apTheta.empty()) {
        for (int i = 0; i < size; i++) {
            apTheta.push_back(-1);
        }
    } else if (apTheta.size() < size) {

        // Also allocate new memory as maybe they are not the same size.
        for (size_t i = apTheta.size(); i < size; i++) {
            apTheta.push_back(0);
        }
    }
    return apTheta;
}
