
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
#include <configurations/Configurations.hpp>

using namespace exageostat::linearAlgebra;
using namespace exageostat::generators::synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations;

using namespace std;

template<typename T>
SyntheticGenerator<T> *
SyntheticGenerator<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new SyntheticGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
SyntheticGenerator<T>::SyntheticGenerator() {

    auto configurations = Configurations::GetConfigurations();
    this->mpLocations = new Locations((configurations->GetProblemSize() * configurations->GetTimeSlot()),
                                      configurations->GetDimension());

    // Set selected Kernel
    std::string kernel_name = configurations->GetKernelName();
    this->mpKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
            configurations->GetKernelName());
    this->mpKernel->SetPValue(configurations->GetTimeSlot());

    configurations->SetProblemSize(configurations->GetProblemSize() * this->mpKernel->GetPValue());
    configurations->SetP(this->mpKernel->GetPValue());

    int parameters_number = this->mpKernel->GetParametersNumbers();
    configurations->SetParametersNumber(parameters_number);

    configurations->SetLowerBounds(InitTheta(configurations->GetLowerBounds(), parameters_number));
    configurations->SetUpperBounds(InitTheta(configurations->GetUpperBounds(), parameters_number));
    configurations->SetInitialTheta(InitTheta(configurations->GetInitialTheta(), parameters_number));
    configurations->SetTargetTheta(InitTheta(configurations->GetTargetTheta(), parameters_number));

    // Set starting theta with the lower bounds values
    configurations->SetStartingTheta(configurations->GetLowerBounds());

    for (int i = 0; i < parameters_number; i++) {
        if (configurations->GetTargetTheta()[i] != -1) {
            configurations->GetLowerBounds()[i] = configurations->GetTargetTheta()[i];
            configurations->GetUpperBounds()[i] = configurations->GetTargetTheta()[i];
            configurations->GetStartingTheta()[i] = configurations->GetTargetTheta()[i];
        }
    }

    // Set linear Algebra solver
    this->mpLinearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(configurations->GetComputation());
}

template<typename T>
void SyntheticGenerator<T>::GenerateDescriptors() {
    this->mpLinearAlgebraSolver->InitiateDescriptors();
}

template<typename T>
void SyntheticGenerator<T>::GenerateLocations() {

    auto configurations = Configurations::GetConfigurations();
    int p;
    if (this->mpKernel) {
        p = this->mpKernel->GetPValue();
    } else {
        throw std::runtime_error("Error in Allocating Kernel plugin");
    }
    int N = configurations->GetProblemSize() / p;
    this->mpLocations->SetSize(N);

    int index = 0;
    Dimension dimension = configurations->GetDimension();
    int time_slots = configurations->GetTimeSlot();
    this->mpLocations->SetDimension(dimension);


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

    double range_low = -0.4, range_high = 0.4;

    for (auto i = 0; i < rootN && index < N; i++) {
        for (auto j = 0; j < rootN && index < N; j++) {
            if (dimension == Dimension3D) {
                for (auto k = 0; k < rootN && index < N; k++) {
                    this->mpLocations->GetLocationX()[index] =
                            (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    this->mpLocations->GetLocationY()[index] =
                            (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    this->mpLocations->GetLocationZ()[index] =
                            (grid[k] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;

                    index++;
                }
            } else {
                this->mpLocations->GetLocationX()[index] =
                        (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                this->mpLocations->GetLocationY()[index] =
                        (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
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
double SyntheticGenerator<T>::UniformDistribution(const double &aRangeLow, const double &aRangeHigh) {
    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = aRangeHigh - aRangeLow;
    return (myRand * range) + aRangeLow;
}

template<typename T>
void SyntheticGenerator<T>::SortLocations(int &aN) {

    // Some sorting, required by spatial statistics code
    uint16_t x, y, z;
    Dimension dimension = Configurations::GetConfigurations()->GetDimension();
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

    void *descriptorC = Configurations::GetConfigurations()->GetDescriptorC()[0];
    const auto &l1 = this->GetLocations();

    this->mpLinearAlgebraSolver->GenerateObservationsVector(descriptorC, l1, l1, nullptr,
                                                            Configurations::GetConfigurations()->GetInitialTheta(), 0,
                                                            this->mpKernel);
}

template<typename T>
std::vector<double> &SyntheticGenerator<T>::InitTheta(std::vector<double> &apTheta, int &size) {

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

template<typename T>
void SyntheticGenerator<T>::DestroyDescriptors() {
    this->mpLinearAlgebraSolver->DestroyDescriptors();
}

template<typename T>
SyntheticGenerator<T>::~SyntheticGenerator() {
    delete this->mpLinearAlgebraSolver;
    delete this->mpKernel;
    delete this->mpLocations;
}

template<typename T>
void SyntheticGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> SyntheticGenerator<T> *SyntheticGenerator<T>::mpInstance = nullptr;
