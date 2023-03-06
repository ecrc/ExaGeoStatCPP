
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
    /// TODO: Is seed always equal to zero?
    int seed = 0;

    this->GenerateLocations(N, seed);

    return aLocations;
}

void SyntheticGenerator::GenerateLocations(int aN, int aSeed) {

    int index = 0;
    srand(aSeed);

    //Allocate memory
    this->mLocations.SetLocationX((double *) malloc(aN * sizeof(double)));
    this->mLocations.SetLocationY((double *) malloc(aN * sizeof(double)));

    if (this->mpConfigurations->GetDimension() != "2D") {
        this->mLocations.SetLocationX((double *) malloc(aN * sizeof(double)));
    }

    int squareRootN = ceil(sqrt(aN));
    int *pGrid = (int *) calloc((int) squareRootN, sizeof(int));

    for (auto i = 0; i < squareRootN; i++) {
        pGrid[i] = i + 1;
    }

    for (auto i = 0; i < squareRootN && index < aN; i++) {
        for (auto j = 0; j < squareRootN && index < aN; j++) {
            this->mLocations.GetLocationX()[index] = (pGrid[i] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / squareRootN;
            this->mLocations.GetLocationY()[index] = (pGrid[j] - 0.5 + this->UniformDistribution(-0.4, 0.4)) / squareRootN;
            index++;
        }
    }
    free(pGrid);
    SortLocations(aN, this->mLocations);
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

uint32_t SyntheticGenerator::EncodeMorton2(uint32_t x, uint32_t y)
//! Encode two inputs into one
{
    return (this->Part1By1(y) << 1) + this->Part1By1(x);
}
uint32_t SyntheticGenerator::Part1By1(uint32_t x)
//! Spread lower bits of input
{
    x &= 0x0000ffff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x << 8)) & 0x00ff00ff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x << 4)) & 0x0f0f0f0f;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x << 2)) & 0x33333333;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x << 1)) & 0x55555555;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

int SyntheticGenerator::compare_uint32(const void *a, const void *b)
//! Compare two uint32_t
{
    uint32_t _a = *(uint32_t *) a;
    uint32_t _b = *(uint32_t *) b;
    if (_a < _b) return -1;
    if (_a == _b) return 0;
    return 1;
}
uint32_t SyntheticGenerator::Compact1By1(uint32_t x)
//! Collect every second bit into lower part of input
{
    x &= 0x55555555;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

void SyntheticGenerator::SortLocations(int aN, Locations aLocations) {

    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y;
    uint32_t z[aN];
    // Encode data into vector z
    for (i = 0; i < aN; i++) {
        x = (uint16_t) (this->mLocations.GetLocationX()[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (this->mLocations.GetLocationY()[i] * (double) UINT16_MAX + .5);
        z[i] = EncodeMorton2(x, y);
    }
    // Sort vector z
    qsort(z, aN, sizeof(uint32_t), compare_uint32);
    // Decode data from vector z
    for (i = 0; i < aN; i++) {
        x = Compact1By1(z[i] >> 0);
        y = Compact1By1(z[i] >> 1);
        this->mLocations.GetLocationX()[i] = (double) x / (double) UINT16_MAX;
        this->mLocations.GetLocationY()[i] = (double) y / (double) UINT16_MAX;
    }
}
