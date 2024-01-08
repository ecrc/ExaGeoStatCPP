
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestCSVDataGenerator.cpp
 * @brief Unit tests for the CSVDataGenerator class in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the CSVDataGenerator class
 * in the ExaGeoStat software package. The tests cover various aspects of data generation, including spreading
 * and reversing bits, generating locations for different dimensions, and testing helper functions.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-03-08
**/

#include <iostream>

#include <catch2/catch_all.hpp>
#include <data-generators/DataGenerator.hpp>
#include <configurations/Configurations.hpp>

using namespace std;

using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::kernels;

void TEST_CSV_P_1() {
    int N = 16;
    string write_path = PROJECT_SOURCE_DIR;
    write_path = write_path + "tests/cpp-tests/data-generators/concrete";
    string read_path = PROJECT_SOURCE_DIR;
    read_path = read_path + +"tests/cpp-tests/data-generators/concrete/synthetic_ds/SYN_16_1";

    Configurations configurations;
    configurations.SetIsCSV(true);
    configurations.SetIsSynthetic(false);
    configurations.SetProblemSize(16);
    configurations.SetDenseTileSize(8);
    configurations.SetKernelName("UnivariateMaternStationary");
    configurations.SetComputation(exageostat::common::EXACT_DENSE);
    configurations.SetDataPath(read_path);

    vector<double> initial_theta{1, 0.1, 0.5};
    configurations.SetInitialTheta(initial_theta);

    Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<Kernel<double>>::Create(configurations.GetKernelName(), configurations.GetTimeSlot());

    auto hardware = exageostat::hardware::ExaGeoStatHardware(configurations.GetComputation(),
                                                             configurations.GetCoresNumber(),
                                                             configurations.GetGPUsNumbers());

    //creating locations x and y.
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                     0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                     0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                     0.877592126344701295};

    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                     0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                     0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                     0.942824444953078489};

    //creating measurements matrix.
    auto *measurements_matrix = new double[N]{-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
                                              -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
                                              0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
                                              -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                                              0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
                                              0.290822066007430102};

    SECTION("Test CSV 2 Dimensions locations, p = 1.") {
        int p = 1;
        configurations.SetDimension(Dimension2D);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, Dimension2D);

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N, p, write_path, locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
        }
    }

    SECTION("Test CSV 3 Dimensions locations, p = 1.") {
        int p = 1;
        configurations.SetDimension(Dimension3D);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(configurations);

        Locations<double> locations(N, Dimension3D);

        auto *location_z = new double[N]{1, 1, 1, 1,
                                         1, 1, 1, 1,
                                         1, 1, 1, 1,
                                         1, 1, 1, 1};

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);
        locations.SetLocationZ(*location_z, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N, p, write_path, locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();
        auto data_loc_z = data->GetLocations()->GetLocationZ();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_z[i] - location_z[i] == Catch::Approx(0.0).margin(1e-6));
        }
        delete[] location_z;
    }

    SECTION("Test CSV ST Dimensions locations, p = 1.") {
        configurations.SetDimension(DimensionST);
        int p = 1;
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(configurations);

        Locations<double> locations(N, DimensionST);

        auto *location_time = new double[N]{1, 1, 1, 1,
                                            1, 1, 1, 1,
                                            1, 1, 1, 1,
                                            1, 1, 1, 1};

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);
        locations.SetLocationZ(*location_time, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N, p, write_path, locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();
        auto data_loc_time = data->GetLocations()->GetLocationZ();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_time[i] - location_time[i] == Catch::Approx(0.0).margin(1e-6));
        }
        delete[] location_time;
    }

    delete[] location_x;
    delete[] location_y;
    delete[] measurements_matrix;
    delete pKernel;
}

void TEST_CSV_P_2() {
    int N = 16;
    string write_path = PROJECT_SOURCE_DIR;
    write_path = write_path + "tests/cpp-tests/data-generators/concrete";
    string read_path = PROJECT_SOURCE_DIR;
    read_path = read_path + +"tests/cpp-tests/data-generators/concrete/synthetic_ds/SYN_16_1";

    Configurations configurations;
    configurations.SetIsCSV(true);
    configurations.SetIsSynthetic(false);
    configurations.SetProblemSize(16);
    configurations.SetDenseTileSize(8);
    configurations.SetKernelName("BivariateMaternParsimonious");
    configurations.SetComputation(exageostat::common::EXACT_DENSE);
    configurations.SetDataPath(read_path);

    Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<Kernel<double>>::Create(configurations.GetKernelName(), configurations.GetTimeSlot());
    int p = pKernel->GetP();

    auto hardware = exageostat::hardware::ExaGeoStatHardware(configurations.GetComputation(),
                                                             configurations.GetCoresNumber(),
                                                             configurations.GetGPUsNumbers());

    //creating locations x and y.
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                     0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                     0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                     0.877592126344701295};

    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                     0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                     0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                     0.942824444953078489};

    auto *measurements_matrix = new double[N * p]{-1.27233614036019, -2.35150374494046,
                                                  .623639758366895, -0.066869877091063,
                                                  .41584737961021, -1.55886176806236,
                                                  .222741609165895, 0.776901795224594,
                                                  1.52805745908099, -0.807621433158528,
                                                  1.62936768085431, 0.291111981187859,
                                                  12028165177046, 0.506286682956441,
                                                  0.930042552942533, 0.303787789424809,
                                                  .579026608633258, -1.44558492413659,
                                                  .27477774798935, -0.644413978717501,
                                                  0.795653357849832, -0.948531131756453,
                                                  0.407532390364584, 0.508624646611594,
                                                  .837631413486745, -1.94405298285362,
                                                  .245007169878023, -1.17842247843059,
                                                  .835589664598579, -2.30510437980008,
                                                  0.545185925535392, -1.08908827912583};

    SECTION("Test CSV 2 Dimensions locations, p = 2.") {
        vector<double> initial_theta_p_2{1, 1, 0.1, 0.5, 0.5, 0.1};
        configurations.SetInitialTheta(initial_theta_p_2);

        configurations.SetDimension(Dimension2D);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, Dimension2D);

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N * p, p, write_path,
                                                                    locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
        }
    }

    SECTION("Test CSV 3 Dimensions locations, p = 2.") {
        vector<double> initial_theta_p_2{1, 1, 0.1, 0.5, 0.5, 0.1};
        configurations.SetInitialTheta(initial_theta_p_2);

        configurations.SetDimension(Dimension3D);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, Dimension3D);

        auto *location_z = new double[N]{1, 1, 1, 1,
                                         1, 1, 1, 1,
                                         1, 1, 1, 1,
                                         1, 1, 1, 1};

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);
        locations.SetLocationZ(*location_z, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N * p, p, write_path,
                                                                    locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();
        auto data_loc_z = data->GetLocations()->GetLocationZ();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_z[i] - location_z[i] == Catch::Approx(0.0).margin(1e-6));
        }
        delete[] location_z;
    }

    SECTION("Test CSV ST Dimensions locations, p = 2.") {
        vector<double> initial_theta_p_2{1, 1, 0.1, 0.5, 0.5, 0.1};
        configurations.SetInitialTheta(initial_theta_p_2);

        configurations.SetDimension(DimensionST);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, DimensionST);

        auto *location_time = new double[N]{1, 1, 1, 1,
                                            1, 1, 1, 1,
                                            1, 1, 1, 1,
                                            1, 1, 1, 1};

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);
        locations.SetLocationZ(*location_time, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N * p, p, write_path,
                                                                    locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();
        auto data_loc_time = data->GetLocations()->GetLocationZ();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_time[i] - location_time[i] == Catch::Approx(0.0).margin(1e-6));
        }
        delete[] location_time;
    }

    delete[] location_x;
    delete[] location_y;
    delete[] measurements_matrix;
    delete pKernel;
}

void TEST_CSV_P_3() {
    int N = 16;
    string write_path = PROJECT_SOURCE_DIR;
    write_path = write_path + "tests/cpp-tests/data-generators/concrete";
    string read_path = PROJECT_SOURCE_DIR;
    read_path = read_path + +"tests/cpp-tests/data-generators/concrete/synthetic_ds/SYN_16_1";

    Configurations configurations;
    configurations.SetIsCSV(true);
    configurations.SetIsSynthetic(false);
    configurations.SetProblemSize(16);
    configurations.SetDenseTileSize(3);
    configurations.SetKernelName("TrivariateMaternParsimonious");
    configurations.SetComputation(exageostat::common::EXACT_DENSE);
    configurations.SetDataPath(read_path);
    vector<double> initial_theta{1, 1, 1, 0.1, 0.5, 1, 1.5, 0.1, 0.1, 0};
    configurations.SetInitialTheta(initial_theta);

    Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<Kernel<double>>::Create(configurations.GetKernelName(), configurations.GetTimeSlot());
    int p = pKernel->GetP();

    auto hardware = exageostat::hardware::ExaGeoStatHardware(configurations.GetComputation(),
                                                             configurations.GetCoresNumber(),
                                                             configurations.GetGPUsNumbers());

    //creating locations x and y.
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                     0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                     0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                     0.877592126344701295};

    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                     0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                     0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                     0.942824444953078489};

    auto *measurements_matrix = new double[N * p]{-1.27233614036019, -2.46098629036509, 0.530373965839942,
                                                  -0.388093615978819, -0.72087019951912, -0.938307444110507,
                                                  0.168800205967935, 0.588775175666892, -1.51480692753486,
                                                  -0.834396737107814, -1.52008593884929, -0.643173921670331,
                                                  1.07615959493385, 0.491034244827904, -1.0359774407964,
                                                  0.489992246076744, 0.749032032301894, -1.72043252385337,
                                                  2.17036429629696, -0.819513652833177, 0.479739022782627,
                                                  -0.636218843975581, -0.277047000397385, -0.259928555173475,
                                                  -0.819522989966449, -1.58696558364785, -0.788818078506206,
                                                  -1.33356599226017, -0.251556901407877, -2.19320893826632,
                                                  -0.584864781619095, -1.21161512739778, 0.285565941432617,
                                                  -2.40856161553517, -0.449549224795429, -1.28344730821537,
                                                  .638514678459714, 0.422008673882837, -0.434734699353173,
                                                  -0.54244903139338, 1.87311710862092, -2.24573298002902,
                                                  -0.493980914795387, -0.0396955273481424, -1.48622486098038,
                                                  0.988614043949225, 1.30567487689688, -2.85832585549073};

    SECTION("Test CSV 2 Dimensions locations, p = 3.") {
        configurations.SetDimension(Dimension2D);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, Dimension2D);

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N * p, p, write_path,
                                                                    locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
        }

    }

    SECTION("Test CSV 3 Dimensions locations, p = 3.") {
        configurations.SetDimension(Dimension3D);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, Dimension3D);

        auto *location_z = new double[N]{1, 1, 1, 1,
                                         1, 1, 1, 1,
                                         1, 1, 1, 1,
                                         1, 1, 1, 1};

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);
        locations.SetLocationZ(*location_z, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N * p, p, write_path,
                                                                    locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();
        auto data_loc_z = data->GetLocations()->GetLocationZ();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_z[i] - location_z[i] == Catch::Approx(0.0).margin(1e-6));
        }

        delete[] location_z;
    }

    SECTION("Test CSV ST Dimensions locations, p = 3.") {
        configurations.SetDimension(DimensionST);
        unique_ptr<DataGenerator<double>> csv_reader = DataGenerator<double>::CreateGenerator(
                configurations);

        Locations<double> locations(N, DimensionST);

        auto *location_time = new double[N]{1, 1, 1, 1,
                                            1, 1, 1, 1,
                                            1, 1, 1, 1,
                                            1, 1, 1, 1};

        locations.SetLocationX(*location_x, N);
        locations.SetLocationY(*location_y, N);
        locations.SetLocationZ(*location_time, N);

        exageostat::helpers::DiskWriter<double>::WriteVectorsToDisk(*measurements_matrix, N * p, p, write_path,
                                                                    locations);
        auto data = csv_reader->CreateData(configurations, hardware, *pKernel);

        auto z_desc_mat = data->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR,
                                                                         data->GetDescriptorData()->GetDescriptor(
                                                                                 CHAMELEON_DESCRIPTOR,
                                                                                 DESCRIPTOR_Z).chameleon_desc);
        auto data_loc_x = data->GetLocations()->GetLocationX();
        auto data_loc_y = data->GetLocations()->GetLocationY();
        auto data_loc_time = data->GetLocations()->GetLocationZ();

        for (int i = 0; i < N; i++) {
            REQUIRE(z_desc_mat[i] - measurements_matrix[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_x[i] - location_x[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_y[i] - location_y[i] == Catch::Approx(0.0).margin(1e-6));
            REQUIRE(data_loc_time[i] - location_time[i] == Catch::Approx(0.0).margin(1e-6));
        }
        delete[] location_time;
    }

    delete[] location_x;
    delete[] location_y;
    delete[] measurements_matrix;
    delete pKernel;
}


TEST_CASE("CSV Data Generation tests") {
    TEST_CSV_P_1();
    TEST_CSV_P_2();
    TEST_CSV_P_3();
}