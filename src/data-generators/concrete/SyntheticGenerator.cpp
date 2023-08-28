
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <kernels/Kernel.hpp>
#include <common/PluginRegistry.hpp>
#include <configurations/Configurations.hpp>

using namespace std;

using namespace exageostat::generators::synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::kernels;

template<typename T>
SyntheticGenerator<T> *
SyntheticGenerator<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new SyntheticGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
Locations<T> *SyntheticGenerator<T>::CreateLocationsData(configurations::Configurations &aConfigurations) {

    // Allocated new Locations object.
    auto *locations = new Locations<T>((aConfigurations.GetProblemSize() * aConfigurations.GetTimeSlot()),
                                       aConfigurations.GetDimension());

    // Register and create a kernel object
    Kernel<T> *kernel = exageostat::plugins::PluginRegistry<Kernel<T>>::Create(aConfigurations.GetKernelName());

    // Set some kernel and arguments values.
    kernel->SetPValue(aConfigurations.GetTimeSlot());
    aConfigurations.SetProblemSize(aConfigurations.GetProblemSize() * kernel->GetPValue());
    aConfigurations.SetP(kernel->GetPValue());
    int parameters_number = kernel->GetParametersNumbers();

    // Set theta's values.
    aConfigurations.SetLowerBounds(InitTheta(aConfigurations.GetLowerBounds(), parameters_number));
    aConfigurations.SetUpperBounds(InitTheta(aConfigurations.GetUpperBounds(), parameters_number));
    aConfigurations.SetInitialTheta(InitTheta(aConfigurations.GetInitialTheta(), parameters_number));
    aConfigurations.SetTargetTheta(InitTheta(aConfigurations.GetTargetTheta(), parameters_number));

    // Set starting theta with the lower bounds values
    aConfigurations.SetStartingTheta(aConfigurations.GetLowerBounds());
    for (int i = 0; i < parameters_number; i++) {
        if (aConfigurations.GetTargetTheta()[i] != -1) {
            aConfigurations.GetLowerBounds()[i] = aConfigurations.GetTargetTheta()[i];
            aConfigurations.GetUpperBounds()[i] = aConfigurations.GetTargetTheta()[i];
            aConfigurations.GetStartingTheta()[i] = aConfigurations.GetTargetTheta()[i];
        }
    }

    // Generate Locations.
    int N = aConfigurations.GetProblemSize() / kernel->GetPValue();
    GenerateLocations(N, aConfigurations.GetTimeSlot(), aConfigurations.GetDimension(), *locations);

    delete kernel;
    return locations;
}

template<typename T>
void SyntheticGenerator<T>::GenerateLocations(const int &aN, const int &aTimeSlot, const Dimension &aDimension,
                                              Locations<T> &aLocations) {

    aLocations.SetSize(aN);

    int index = 0;
    aLocations.SetDimension(aDimension);


    int rootN;
    if (aDimension == Dimension3D) {
        //Cubic root.
        rootN = ceil(cbrt(aN));
    } else {
        //Square root.
        rootN = ceil(sqrt(aN));
    }

    int *grid = new int[rootN]();
    for (auto i = 0; i < rootN; i++) {
        grid[i] = i + 1;
    }

    double range_low = -0.4, range_high = 0.4;

    for (auto i = 0; i < rootN && index < aN; i++) {
        for (auto j = 0; j < rootN && index < aN; j++) {
            if (aDimension == Dimension3D) {
                for (auto k = 0; k < rootN && index < aN; k++) {
                    aLocations.GetLocationX()[index] =
                            (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    aLocations.GetLocationY()[index] =
                            (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    aLocations.GetLocationZ()[index] =
                            (grid[k] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                    index++;
                }
            } else {
                aLocations.GetLocationX()[index] =
                        (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                aLocations.GetLocationY()[index] =
                        (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                if (aDimension == DimensionST) {
                    aLocations.GetLocationZ()[index] = 1.0;
                }
                index++;
            }
        }
    }
    delete[] grid;
    if (aDimension != DimensionST) {
        SortLocations(aN, aDimension, aLocations);
    } else {
        for (auto j = 1; j < aTimeSlot; j++) {
            for (auto i = 0; i < aN; i++) {
                aLocations.GetLocationX()[i + j * aN] = aLocations.GetLocationX()[i];
                aLocations.GetLocationY()[i + j * aN] = aLocations.GetLocationY()[i];
                aLocations.GetLocationZ()[i + j * aN] = (double) (j + 1);
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
void SyntheticGenerator<T>::SortLocations(const int &aN, const Dimension &aDimension, Locations<T> &aLocations) {

    // Some sorting, required by spatial statistics code
    uint16_t x, y, z;
    uint64_t vectorZ[aN];

    // Encode data into vector z
    for (auto i = 0; i < aN; i++) {
        x = (uint16_t) (aLocations.GetLocationX()[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (aLocations.GetLocationY()[i] * (double) UINT16_MAX + .5);
        if (aDimension != Dimension2D) {
            z = (uint16_t) (aLocations.GetLocationZ()[i] * (double) UINT16_MAX + .5);
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
        aLocations.GetLocationX()[i] = (double) x / (double) UINT16_MAX;
        aLocations.GetLocationY()[i] = (double) y / (double) UINT16_MAX;
        if (aDimension == Dimension3D) {
            aLocations.GetLocationZ()[i] = (double) z / (double) UINT16_MAX;
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
std::vector<double> &SyntheticGenerator<T>::InitTheta(std::vector<double> &aTheta, const int &size) {

    // If null, this mean user have not passed the values arguments, Make values equal -1
    if (aTheta.empty()) {
        for (int i = 0; i < size; i++) {
            aTheta.push_back(-1);
        }
    } else if (aTheta.size() < size) {

        // Also allocate new memory as maybe they are not the same size.
        for (size_t i = aTheta.size(); i < size; i++) {
            aTheta.push_back(0);
        }
    }
    return aTheta;
}

template<typename T>
void SyntheticGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> SyntheticGenerator<T> *SyntheticGenerator<T>::mpInstance = nullptr;
