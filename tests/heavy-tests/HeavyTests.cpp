
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatApi.cpp
 * @brief Test suite for the ExaGeoStat APIs data generation functionality.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-12-20
**/

#include <iostream>
#include <thread>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#include <catch2/catch_all.hpp>

#include <kernels/Kernel.hpp>

using namespace std;

using namespace exageostat::kernels;

void GenerateCommandLineArguments(const string &aKernelName, const string &aComputation, vector<string> &aDimensions,
                                  const string &aPrecision, vector<string> &arguments_vector,
                                  const string &aDistanceType) {

    // Search for the substring "TimeSpace" in the kernel name, As ST dimension only works with kernels having Spacetime in their name.
    // This is a required convention described in the user manual
    size_t found = aKernelName.find("Spacetime");
    // Check if the substring was found
    if (found != std::string::npos) {
        aDimensions = {"ST"};
    } else {
        aDimensions = {"2D", "3D"};
    }

    // Set up a random number generator
    random_device rd;
    mt19937 gen(rd());
    int LOWER_SIZE = 6;
    int MAX_SIZE = 25;
    int MAX_LTS = 18;
    int lts = INT16_MAX;
    uniform_int_distribution<int> problem_size_distribution(LOWER_SIZE, MAX_SIZE);
    int max_threads_size = static_cast<int>(std::thread::hardware_concurrency());
    uniform_int_distribution<int> cpu_size_distribution(1, max_threads_size - 5);
    uniform_int_distribution<int> max_iteration_distribution(1, 4);
    uniform_int_distribution<int> band_distribution(1, 6);
    uniform_int_distribution<int> tolerance_distribution(1, 5);
    uniform_int_distribution<int> max_rank_distribution(200, 600);
    uniform_int_distribution<int> seed_distribution(0, 3);

    int N_value = problem_size_distribution(gen);

    if (aComputation == "tlr") {
        uniform_int_distribution<int> lts_distribution(LOWER_SIZE, MAX_LTS);
        lts = lts_distribution(gen);
        N_value = lts * lts;
        if (aKernelName.find("Bivariate") != string::npos || aKernelName.find("Trivariate") != string::npos) {
            lts = 18;
        }
        arguments_vector.push_back("--lts=" + to_string(lts));
    }

    // tile size need to be smaller than N or lts in case of HiCMA.
    uniform_int_distribution<int> dts_distribution(LOWER_SIZE, min(N_value, lts));
    int dts = dts_distribution(gen);

    // Since both Bivariate and Trivariate requires specific values.
    if (aKernelName.find("Bivariate") != string::npos || aKernelName.find("Trivariate") != string::npos) {
        dts = 18;
        N_value = 162;
    }
    if (aKernelName.find("Trivariate") != string::npos) {
        N_value = 108;
    }
    arguments_vector.push_back("--dts=" + to_string(dts));
    arguments_vector.push_back("--N=" + to_string(N_value));

    uniform_int_distribution<int> zMiss_distribution(2, min(N_value / 2, 100));
    // Create a kernel object to get the number of parameters for each kernel.
    Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<Kernel<double >>::Create(aKernelName, 1);
    int param_number = pKernel->GetParametersNumbers();
    delete pKernel;

    string initial_theta_arguments, lb_arguments, ub_arguments;
    int zmiss;
    zmiss = zMiss_distribution(gen);
    if (aKernelName.find("Bivariate") != string::npos) {
        do {
            zmiss = zMiss_distribution(gen);
        } while (zmiss % 2 != 0);
    }

    if (aKernelName == "BivariateMaternFlexible") {
        initial_theta_arguments = "0.3:0.6:0.01:0.3:0.9:0.9:0.05:0.3:1.5:0.9:0.99";
        lb_arguments = "0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01";
        ub_arguments = "50:50:50:50:50:50:50:50:50:50:50";
    } else if (aKernelName == "BivariateMaternParsimonious") {
        initial_theta_arguments = "1:1:0.1:0.5:0.5:0.1";
        lb_arguments = "0.1:0.1:0.1:0.1:0.1:0.1";
        ub_arguments = "5:5:5:5:5:5";
    } else if (aKernelName == "BivariateSpacetimeMaternStationary" || aKernelName == "TrivariateMaternParsimonious") {
        initial_theta_arguments = "1:1:0.1:0.5:0.5:0.1:0:0:0:0";
        lb_arguments = "0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1";
        ub_arguments = "5:5:5:5:5:5:5:5:5:5";
    } else {
        uniform_real_distribution<double> initial_theta_distribution(0.01, 2);
        // Upper bond need to be larger than the initial theta
        uniform_real_distribution<double> ub_distribution(2, 5);

        for (int n = 0; n < param_number; n++) {
            double initial_theta_value = initial_theta_distribution(gen);
            initial_theta_arguments += to_string(initial_theta_value);
            double lb_temp;
            do {
                lb_temp = initial_theta_distribution(gen);
            } while (lb_temp >= initial_theta_value);
            lb_arguments += to_string(lb_temp);
            ub_arguments += to_string(ub_distribution(gen));
            if (n < param_number - 1) {
                initial_theta_arguments += ":";
                ub_arguments += ":";
                lb_arguments += ":";
            }
        }
    }

    arguments_vector.push_back("--cores=" + to_string(cpu_size_distribution(gen)));
#ifdef USE_CUDA
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    uniform_int_distribution<int> gpu_size_distribution(1, nDevices);
    arguments_vector.push_back("--gpus=" + to_string(gpu_size_distribution(gen)));
#endif

    arguments_vector.push_back("--kernel=" + aKernelName);
    arguments_vector.push_back("--computation=" + aComputation);
    arguments_vector.push_back("--precision=" + aPrecision);
    arguments_vector.push_back("--initial_theta=" + initial_theta_arguments);
    arguments_vector.push_back("--ub=" + ub_arguments);
    arguments_vector.push_back("--lb=" + lb_arguments);
    arguments_vector.push_back("--max_mle_iterations=" + to_string(max_iteration_distribution(gen)));
    arguments_vector.push_back("--tolerance=" + to_string(tolerance_distribution(gen)));
    arguments_vector.push_back("--max_rank=" + to_string(max_rank_distribution(gen)));
    arguments_vector.push_back("--band=" + to_string(band_distribution(gen)));
    arguments_vector.push_back("--ZMiss=" + to_string(zmiss));
    arguments_vector.push_back("--seed=" + to_string(seed_distribution(gen)));
    arguments_vector.emplace_back("--idw");
    arguments_vector.emplace_back("--distance_metric=" + aDistanceType);
    if (aKernelName.find("Bivariate") == string::npos && aKernelName.find("Trivariate") == string::npos) {
        arguments_vector.emplace_back("--fisher");
    }
    if (aKernelName.find("Trivariate") == string::npos) {
        arguments_vector.emplace_back("--mspe");
        arguments_vector.emplace_back("--mloe-mmom");
    }
}

TEST_CASE("END-TO-END") {

    // Add all the possible combination of code.
    vector<string> kernels_vector = {
            "BivariateMaternFlexible",
            "UnivariateMaternStationary",
            "BivariateMaternParsimonious",
            "BivariateSpacetimeMaternStationary",
            "TrivariateMaternParsimonious",
            "UnivariateExpNonGaussian",
            "UnivariateMaternNonGaussian",
            "UnivariateMaternNuggetsStationary",
            "UnivariateSpacetimeMaternStationary"
    };

    vector<string> computation_vector = {"exact", "diag_approx"};
#ifdef USE_HICMA
    computation_vector.emplace_back("tlr");
#endif
    vector<string> dimensions_vector;
    vector<string> precision_vector = {"single", "double"};
    vector<string> distance_vector = {"eg", "gcd"};

    size_t combination_number = 1;
    size_t number_of_iterations = 1;
    for (int i = 0; i < number_of_iterations; i++) {
        cout << "**** END-TO_END TESTS -- ITERATION NUMBER (" << i + 1 << ") ****" << endl;

        // Generate combinations.
        for (const auto &current_kernel: kernels_vector) {
            for (const auto &computation: computation_vector) {
                for (const auto &precision: precision_vector) {
                    for (const auto &distance_type: distance_vector) {

                        vector<string> arguments_vector;
                        // This helper function will fill the arguments vector with the software arguments commands
                        GenerateCommandLineArguments(current_kernel, computation, dimensions_vector, precision,
                                                     arguments_vector, distance_type);
                        // Generate combination of dimensions.
                        for (const auto &dimension: dimensions_vector) {
                            // Add the dimension into the arguments vector
                            arguments_vector.push_back("--dimension=" + dimension);

                            // Great Circle (GC) distance is only valid for 2D!
                            if (distance_type == "gcd" && (dimension == "3D" || dimension == "ST")) {
                                continue;
                            }
                            if (dimension == "ST" && computation == "tlr") {
                                continue;
                            }
                            if (dimension == "ST") {
                                random_device rd;
                                mt19937 gen(rd());
                                uniform_int_distribution<int> time_slot_distribution(1, 5);
                                arguments_vector.push_back("--time_slot=" + to_string(time_slot_distribution(gen)));
                            }

                            string arguments_string;
                            for (const auto &j: arguments_vector) {
                                arguments_string += " " + j;
                            }

                            cout << "(" << combination_number
                                 << ") Testing the software with arguments: " + arguments_string << endl << flush;

                            // Specify the executable path as the first argument
                            std::filesystem::path currentPath =
                                    std::filesystem::current_path().parent_path().parent_path() /
                                    "examples/end-to-end/Example_Data_Generation_Modeling_and_Prediction";
                            std::string fullCommand = std::string(currentPath) + arguments_string;
                            if (std::system(fullCommand.c_str())) {
                                throw runtime_error("This test failed " + fullCommand);
                            }
                            // To remove the previous dimension.
                            arguments_vector.pop_back();

                            if (dimension == "ST") {
                                // To remove the time slot.
                                arguments_vector.pop_back();
                            }
                            combination_number += 1;
                        }
                    }
                }
            }
        }
    }
}

