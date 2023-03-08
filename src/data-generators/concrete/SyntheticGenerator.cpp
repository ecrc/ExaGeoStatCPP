
/**
 * @file SyntheticGenerator.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <cmath>

using namespace exageostat::generators::Synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::configurations::data_configurations;

Locations SyntheticGenerator::InitializeLocations(Locations aLocations) {

    int p = 1;
    int N = this->mpConfigurations->GetProblemSize() / p;

    this->GenerateLocations(N);
    return aLocations;
}

void SyntheticGenerator::GenerateLocations(int aN, int aTimeSlots) {

    int index = 0, aSeed = 0;
    std::string dimension = this->mpConfigurations->GetDimension();
    srand(aSeed);

    //Allocate memory
    this->mLocations.SetLocationX((double *) malloc(aN * aTimeSlots * sizeof(double)));
    this->mLocations.SetLocationY((double *) malloc(aN * aTimeSlots * sizeof(double)));

    if (dimension != "2D") {
        this->mLocations.SetLocationZ((double *) malloc(aN * aTimeSlots * sizeof(double)));
    }

    int rootN;

    if (dimension == "3D"){
        // Cubic root.
        rootN = ceil(cbrt(aN));
    }
    else{
        // Square root.
        rootN = ceil(sqrt(aN));
    }

    int *grid = (int *) calloc((int) rootN, sizeof(int));

    for (auto i = 0; i < rootN; i++) {
        grid[i] = i + 1;
    }

    for (auto i = 0; i < rootN && index < aN; i++) {
        for (auto j = 0; j < rootN && index < aN; j++) {
            if (dimension == "3D"){
                for (auto k = 0; k < rootN && index < aN; k++) {
                    this->mLocations.GetLocationX()[index] = (grid[i] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                    this->mLocations.GetLocationY()[index] = (grid[j] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                    this->mLocations.GetLocationZ()[index] = (grid[k] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                    index++;
                }
            }
            else{
                this->mLocations.GetLocationX()[index] = (grid[i] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                this->mLocations.GetLocationY()[index] = (grid[j] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                if (dimension == "ST"){
                    this->mLocations.GetLocationZ()[index] = 1.0;
                }
                index++;
            }
        }
    }
    free(grid);
    if (dimension != "ST"){
        SortLocations(aN, this->mLocations);
    }
    else{
        for (auto j = 1; j < aTimeSlots; j++) {
            for (auto i = 0; i < aN; i++) {
                this->mLocations.GetLocationX()[i + j * aN] = this->mLocations.GetLocationX()[i];
                this->mLocations.GetLocationY()[i + j * aN] = this->mLocations.GetLocationY()[i];
                this->mLocations.GetLocationZ()[i + j * aN] = (double) (j + 1);
            }
        }
    }
}

void SyntheticGenerator::Print() {
    std::cout << "HELLO YOU'RE USING SYNTHETIC DATA GENERATION" << std::endl;
    std::cout << "N: " << mpConfigurations->GetProblemSize() << std::endl;
}

SyntheticGenerator::SyntheticGenerator(configurations::data_configurations::SyntheticDataConfigurations *apConfigurations) {
    // Set configuration map.
    this->SetConfigurations(apConfigurations);
}

double SyntheticGenerator::UniformDistribution(double aRangeLow, double aRangeHigh) {

    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = aRangeHigh - aRangeLow;

    return (myRand * range) + aRangeLow;
}

uint32_t SyntheticGenerator::SpreadBits(uint32_t aInputByte)
{
    /// TODO: Ask sameh if we need this for 2D or it works with 3D implementation.
//    aInputByte &= 0x0000ffff;
//    // aInputByte = ---- ---- ---- ---- fedc ba98 7654 3210
//    aInputByte = (aInputByte ^ (aInputByte << 8)) & 0x00ff00ff;
//    // aInputByte = ---- ---- fedc ba98 ---- ---- 7654 3210
//    aInputByte = (aInputByte ^ (aInputByte << 4)) & 0x0f0f0f0f;
//    // aInputByte = ---- fedc ---- ba98 ---- 7654 ---- 3210
//    aInputByte = (aInputByte ^ (aInputByte << 2)) & 0x33333333;
//    // aInputByte = --fe --dc --ba --98 --76 --54 --32 --10
//    aInputByte = (aInputByte ^ (aInputByte << 1)) & 0x55555555;
//    // aInputByte = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
//    return aInputByte;


    aInputByte &= 0x000000000000ffff;
    // aInputByte = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    aInputByte = (aInputByte ^ (aInputByte << 24)) & 0x000000ff000000ff;
    // aInputByte = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    aInputByte = (aInputByte ^ (aInputByte << 12)) & 0x000f000f000f000f;
    // aInputByte = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    aInputByte = (aInputByte ^ (aInputByte << 6)) & 0x0303030303030303;
    // aInputByte = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    aInputByte = (aInputByte ^ (aInputByte << 3)) & 0x1111111111111111;
    // aInputByte = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    return aInputByte;
}

uint32_t SyntheticGenerator::ReverseSpreadBits(uint32_t aInputByte)
{
//    // Collect every second bit into lower part of input
//    aInputByte &= 0x55555555;
//    // aInputByte = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
//    aInputByte = (aInputByte ^ (aInputByte >> 1)) & 0x33333333;
//    // aInputByte = --fe --dc --ba --98 --76 --54 --32 --10
//    aInputByte = (aInputByte ^ (aInputByte >> 2)) & 0x0f0f0f0f;
//    // aInputByte = ---- fedc ---- ba98 ---- 7654 ---- 3210
//    aInputByte = (aInputByte ^ (aInputByte >> 4)) & 0x00ff00ff;
//    // aInputByte = ---- ---- fedc ba98 ---- ---- 7654 3210
//    aInputByte = (aInputByte ^ (aInputByte >> 8)) & 0x0000ffff;
//    // aInputByte = ---- ---- ---- ---- fedc ba98 7654 3210
//    return aInputByte;

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

int SyntheticGenerator::CompareUint32(const void *apFirstValue, const void *apSecondValue) {

    uint32_t firstValue = *(uint32_t *) apFirstValue;
    uint32_t secondValue = *(uint32_t *) apSecondValue;

    if (firstValue == secondValue){
        return 0;
    }
    else if (firstValue > secondValue) {
        return 1;
    }
    // SecondValue > firstValue
    return -1;
}

void SyntheticGenerator::SortLocations(int aN, Locations aLocations) {

    // Some sorting, required by spatial statistics code
    uint16_t x, y, z;
    std::string dimension = this->mpConfigurations->GetDimension();
    uint64_t vectorZ[aN];

    // Encode data into vector z
    for (auto i = 0; i < aN; i++) {
        x = (uint16_t) (this->mLocations.GetLocationX()[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (this->mLocations.GetLocationY()[i] * (double) UINT16_MAX + .5);
        if (dimension == "3D"){
            z = (uint16_t) (this->mLocations.GetLocationZ()[i] * (double) UINT16_MAX + .5);
        }
        else{
            z = (uint16_t) 0.0;
        }
        vectorZ[i] = (this->SpreadBits(z) << 2) + (this->SpreadBits(y) << 1) + this->SpreadBits(x);
    }
    // Sort vector z
    qsort(vectorZ, aN, sizeof(uint64_t), CompareUint32);

    // Decode data from vector z
    for (auto i = 0; i < aN; i++) {
        x = ReverseSpreadBits(vectorZ[i] >> 0);
        y = ReverseSpreadBits(vectorZ[i] >> 1);
        z = ReverseSpreadBits(vectorZ[i] >> 2);
        this->mLocations.GetLocationX()[i] = (double) x / (double) UINT16_MAX;
        this->mLocations.GetLocationY()[i] = (double) y / (double) UINT16_MAX;
        if (dimension == "3D"){
            this->mLocations.GetLocationZ()[i] = (double) z / (double) UINT16_MAX;
        }
    }
}
