
/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <cmath>
#include <algorithm>

using namespace exageostat::generators::Synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;

SyntheticGenerator::SyntheticGenerator(SyntheticDataConfigurations *apConfigurations) {
    // Set configuration map.
    this->SetConfigurations(apConfigurations);
    this->InitLocationsClass();
}

void SyntheticGenerator::InitializeLocations() {

    int p = 1;
    int N = this->mpConfigurations->GetProblemSize() / p;
    this->GenerateLocations(N);
}

void SyntheticGenerator::GenerateLocations(int aN) {

    int index = 0;
    Dimension dimension = this->mpConfigurations->GetDimension();
    int time_slots = this->mpConfigurations->GetTimeSlot();

    //Allocate memory
    this->mpLocations->SetLocationX((double *) malloc(aN * time_slots * sizeof(double)));
    this->mpLocations->SetLocationY((double *) malloc(aN * time_slots * sizeof(double)));

    if (dimension != Dimension2D) {
        this->mpLocations->SetLocationZ((double *) malloc(aN * time_slots * sizeof(double)));
    }

    int rootN;
    if (dimension == Dimension3D){
        //Cubic root.
        rootN = ceil(cbrt(aN));
    }
    else{
         //Square root.
        rootN = ceil(sqrt(aN));
    }

    int *grid = (int *) calloc((int) rootN, sizeof(int));
    for (auto i = 0; i < rootN; i++) {
        grid[i] = i + 1;
    }

    for (auto i = 0; i < rootN && index < aN; i++) {
        for (auto j = 0; j < rootN && index < aN; j++) {
            if (dimension == Dimension3D){
                for (auto k = 0; k < rootN && index < aN; k++) {
                    this->mpLocations->GetLocationX()[index] = (grid[i] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                    this->mpLocations->GetLocationY()[index] = (grid[j] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                    this->mpLocations->GetLocationZ()[index] = (grid[k] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                    index++;
                }
            }
            else{
                this->mpLocations->GetLocationX()[index] = (grid[i] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                this->mpLocations->GetLocationY()[index] = (grid[j] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / rootN;
                if (dimension == DimensionST){
                    this->mpLocations->GetLocationZ()[index] = 1.0;
                }
                index++;
            }
        }
    }
    free(grid);
    if (dimension != DimensionST){
        SortLocations(aN);
    }
    else{
        for (auto j = 1; j < time_slots; j++) {
            for (auto i = 0; i < aN; i++) {
                this->mpLocations->GetLocationX()[i + j * aN] = this->mpLocations->GetLocationX()[i];
                this->mpLocations->GetLocationY()[i + j * aN] = this->mpLocations->GetLocationY()[i];
                this->mpLocations->GetLocationZ()[i + j * aN] = (double) (j + 1);
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
        if (dimension != Dimension2D){
            z = (uint16_t) (this->mpLocations->GetLocationZ()[i] * (double) UINT16_MAX + .5);
        }
        else{
            z = (uint16_t) 0.0;
        }
        vectorZ[i] = (this->SpreadBits(z) << 2) + (this->SpreadBits(y) << 1) + this->SpreadBits(x);
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
        if (dimension == Dimension3D){
            this->mpLocations->GetLocationZ()[i] = (double) z / (double) UINT16_MAX;
        }
    }
}

uint64_t SyntheticGenerator::SpreadBits(uint64_t aInputByte)
{
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

uint64_t SyntheticGenerator::ReverseSpreadBits(uint64_t aInputByte)
{

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
