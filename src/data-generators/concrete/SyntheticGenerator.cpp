
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>

using namespace exageostat::generators::synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::hardware;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::linearAlgebra;

template<typename T>
SyntheticGenerator<T> *SyntheticGenerator<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new SyntheticGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
SyntheticGenerator<T>::CreateData(exageostat::configurations::Configurations &aConfigurations,
                                  const exageostat::hardware::ExaGeoStatHardware &aHardware,
                                  exageostat::kernels::Kernel<T> &aKernel) {

    int n = aConfigurations.GetProblemSize() * aConfigurations.GetTimeSlot();
    auto data = std::make_unique<ExaGeoStatData<T>>(n, aConfigurations.GetDimension());
    // Allocated new Locations object.
    auto *locations = new Locations<T>(n, aConfigurations.GetDimension());

    int parameters_number = aKernel.GetParametersNumbers();

    // Set initial theta values.
    Configurations::InitTheta(aConfigurations.GetInitialTheta(), parameters_number);
    aConfigurations.SetInitialTheta(aConfigurations.GetInitialTheta());

    // Generate Locations phase
    GenerateLocations(n, aConfigurations.GetTimeSlot(), aConfigurations.GetDimension(), *locations);
    data->SetLocations(*locations);

    // Generate Descriptors phase
    auto linear_algebra_solver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(EXACT_DENSE);
    linear_algebra_solver->GenerateSyntheticData(aConfigurations, aHardware, data, aKernel);

    if (aConfigurations.GetLogger()) {
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        auto pMatrix = new T[aConfigurations.GetProblemSize()];
        std::string path = aConfigurations.GetLoggerPath();

        CHAMELEON_Desc2Lap(ChamUpperLower, data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                               DESCRIPTOR_Z).chameleon_desc,  pMatrix, aConfigurations.GetProblemSize());
        if ( CHAMELEON_Comm_rank == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(*((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                               DESCRIPTOR_Z).chameleon_desc->mat),
                                              aConfigurations.GetProblemSize(), aConfigurations.GetVariablesNumber(), path,
                                              *data->GetLocations());
        }
        delete[] pMatrix;
#else
        std::string path = aConfigurations.GetLoggerPath();
        DiskWriter<T>::WriteVectorsToDisk(*((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                           DESCRIPTOR_Z).chameleon_desc->mat),
                                          aConfigurations.GetProblemSize(), aKernel.GetVariablesNumber(), path,
                                          *data->GetLocations());
#endif
        VERBOSE("Done.")
    }
    results::Results::GetInstance()->SetGeneratedLocationsNumber(n);
    results::Results::GetInstance()->SetIsLogger(aConfigurations.GetLogger());
    results::Results::GetInstance()->SetLoggerPath(aConfigurations.GetLoggerPath());

    return data;
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

    T range_low = -0.4, range_high = 0.4;

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
                aLocations.GetLocationX()[index] = (grid[i] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
                aLocations.GetLocationY()[index] = (grid[j] - 0.5 + UniformDistribution(range_low, range_high)) / rootN;
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
        for (auto i = 0; i < aN; i++) {
            aLocations.GetLocationX()[i] = aLocations.GetLocationX()[i];
            aLocations.GetLocationY()[i] = aLocations.GetLocationY()[i];
            aLocations.GetLocationZ()[i] = (T) (i / aTimeSlot + 1);
        }
    }
}

template<typename T>
T SyntheticGenerator<T>::UniformDistribution(const T &aRangeLow, const T &aRangeHigh) {
    T myRand = (T) rand() / (T) (1.0 + RAND_MAX);
    T range = aRangeHigh - aRangeLow;
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
void SyntheticGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> SyntheticGenerator<T> *SyntheticGenerator<T>::mpInstance = nullptr;
