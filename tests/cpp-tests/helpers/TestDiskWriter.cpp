
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestDiskWriter.cpp
 * @brief Unit tests for the DiskWriter in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the class DiskWriter.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2024-01-24
**/

#include <fstream>
#include <catch2/catch_all.hpp>

#include <common/Definitions.hpp>
#include <helpers/DiskWriter.hpp>

using namespace exageostat::dataunits;
using namespace exageostat::common;

void TEST_3D_VECTORS_WRITING() {
    // Initialize data vectors
    int N = 8;
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916};

    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011};

    auto *location_z = new double[N]{1, 1, 1, 1,
                                     1, 1, 1, 1};

    auto *measurements_matrix = new double[N]{-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
                                              -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
                                              0.161705025505231914, 0.623389205185149065};

    // Initialize Locations object, that will be written in file
    Locations<double> locations(N, Dimension3D);
    locations.SetLocationX(*location_x, N);
    locations.SetLocationY(*location_y, N);
    locations.SetLocationZ(*location_z, N);

    const int p = 1;

    std::string write_path = PROJECT_SOURCE_DIR;
    std::string expectedFilePath = write_path + "/synthetic_ds/SYN_" + std::to_string(N / p) + "_1";

    // Write the data into file
    exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N, p, write_path, locations);
    REQUIRE(std::filesystem::exists(expectedFilePath));

    std::ifstream file(expectedFilePath);
    REQUIRE(file.is_open());

    // Read the contents from the file
    std::string line;
    for (int i = 0; i < N / p; ++i) {
        std::getline(file, line);

        std::ostringstream oss;
        oss << std::setprecision(15) << locations.GetLocationX()[i] << ','
            << locations.GetLocationY()[i] << ',' << locations.GetLocationZ()[i] << ","<< std::setprecision(15) << measurements_matrix[i] ;
        std::string expectedLine = oss.str();

        REQUIRE(line == expectedLine);
    }

    delete[] location_x;
    delete[] location_y;
    delete[] location_z;
    delete[] measurements_matrix;
}
TEST_CASE("Disk Writer Tests"){
    TEST_3D_VECTORS_WRITING();
}
