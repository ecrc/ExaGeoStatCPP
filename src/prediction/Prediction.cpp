
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Prediction.cpp
 * @brief Contains the implementation of the Prediction class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <prediction/Prediction.hpp>
#include <prediction/PredictionHelpers.hpp>
#include <prediction/PredictionAuxiliaryFunctions.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace exageostat::prediction;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;

template<typename T>
void Prediction<T>::PredictMissingData(const hardware::ExaGeoStatHardware &aHardware,
                                       std::unique_ptr<dataunits::ExaGeoStatData<T>> &aData,
                                       Configurations &aConfigurations, T *apMeasurementsMatrix,
                                       const kernels::Kernel<T> &aKernel) {

    int i, j;
    bool can_predict = true;
    int num_params = aKernel.GetParametersNumbers();
    for (i = 0; i < num_params; i++) {
        if (aConfigurations.GetEstimatedTheta()[i] == -1) {
            can_predict = false;
            break;
        }
    }
    if (!can_predict &&
        (aConfigurations.GetIsMLOEMMOM() || aConfigurations.GetIsMSPE() || aConfigurations.GetIsFisher())) {
        throw std::runtime_error(
                "Can't predict without an estimated theta, please either pass --etheta or run the modeling module before prediction");
    }

    int number_of_mspe = 3;
    int p = aKernel.GetP();
    int z_miss_number = aConfigurations.GetUnknownObservationsNb();
    int n_z_obs = aConfigurations.CalculateZObsNumber();
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(common::EXACT_DENSE);

    //FISHER Prediction Function Call
    if (aConfigurations.GetIsFisher()) {
        LOGGER("---- Using Prediction Function Fisher ----")
        T *fisher_results;
        fisher_results = linear_algebra_solver->ExaGeoStatFisherTile(aConfigurations, aData, aHardware,
                                                                     (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                     aKernel);

        LOGGER("- Sd of sigma2, alpha, nu: " << sqrt(fisher_results[0]) << " " << sqrt(fisher_results[4]) << " "
                                             << sqrt(fisher_results[8]))
        LOGGER("- CI for sigma2: " << aConfigurations.GetEstimatedTheta()[0] - Q_NORM * sqrt(fisher_results[0]) << " "
                                   << aConfigurations.GetEstimatedTheta()[0] + Q_NORM * sqrt(fisher_results[0]))

        LOGGER("- CI for alpha: " << aConfigurations.GetEstimatedTheta()[1] - Q_NORM * sqrt(fisher_results[4]) << " "
                                  << aConfigurations.GetEstimatedTheta()[1] + Q_NORM * sqrt(fisher_results[4]))
        LOGGER("- CI for nu: " << aConfigurations.GetEstimatedTheta()[2] - Q_NORM * sqrt(fisher_results[8]) << " "
                               << aConfigurations.GetEstimatedTheta()[2] + Q_NORM * sqrt(fisher_results[8]))

        LOGGER("- Fisher Matrix:")
        for (i = 0; i < num_params; i++) {
            LOGGER("  ", true)
            for (j = 0; j < num_params; j++) {
                LOGGER_PRECISION(fisher_results[i * num_params + j], 18)
                if (j != num_params - 1) {
                    LOGGER_PRECISION(", ")
                }
            }
            LOGGER("")
        }
        LOGGER("")
        delete[] fisher_results;
    }

    if (z_miss_number <= 0) {
        return;
    }

    LOGGER("- Total number of Z: " << aConfigurations.GetProblemSize())
    LOGGER("- Number of Z Miss: " << z_miss_number)
    LOGGER("- Number of Z observations: " << n_z_obs)

    results::Results::GetInstance()->SetZMiss(z_miss_number);
    T *z_obs = new T[n_z_obs * p];
    T *z_miss = new T[z_miss_number];
    T *z_actual = new T[z_miss_number * p];
    std::vector<T> avg_pred_value(number_of_mspe);
    auto miss_locations = new Locations<T>(z_miss_number, aConfigurations.GetDimension());
    // Prediction is only supported with 2D.
    auto obs_locations = new Locations<T>(n_z_obs, aConfigurations.GetDimension());

    // We Predict date with only Exact computation. This is a pre-request.
    InitializePredictionArguments(aConfigurations, aData, linear_algebra_solver, z_obs, z_actual, *miss_locations,
                                  *obs_locations, apMeasurementsMatrix, p);

    // MLOE MMOM Auxiliary Function Call
    if (aConfigurations.GetIsMLOEMMOM()) {
        LOGGER("---- Using Auxiliary Function MLOE MMOM ----")
        linear_algebra_solver->ExaGeoStatMLETileMLOEMMOM(aConfigurations, aData, aHardware,
                                                         (T *) aConfigurations.GetInitialTheta().data(),
                                                         (T *) aConfigurations.GetEstimatedTheta().data(),
                                                         *miss_locations, *obs_locations, aKernel);
    }

    // IDW Auxiliary Function Call
    if (aConfigurations.GetIsIDW()) {
        LOGGER("---- Using Auxiliary Function IDW ----")
        T *mspe = new T[number_of_mspe];
        PredictionAuxiliaryFunctions<T>::PredictIDW(z_miss, z_actual, z_obs, z_miss_number, n_z_obs, *miss_locations,
                                                    *obs_locations, mspe);

        LOGGER("- Average prediction Error (IDW): ", true)
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += mspe[i];
            LOGGER_PRECISION(avg_pred_value[i], 9)
            if (i != number_of_mspe - 1) {
                LOGGER_PRECISION(", ")
            }
        }
        results::Results::GetInstance()->SetIDWError(
                std::vector<double>(avg_pred_value.data(), avg_pred_value.data() + number_of_mspe));
        LOGGER("")
        delete[] mspe;
    }

    // MSPE Prediction Function Call
    if (aConfigurations.GetIsMSPE()) {
        LOGGER("---- Using MSPE ----")
        T *prediction_error_mspe;
        if (aConfigurations.GetIsNonGaussian()) {
            prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLENonGaussianPredictTile(aData,
                                                                                               (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                               z_miss_number, n_z_obs,
                                                                                               z_obs,
                                                                                               z_actual, z_miss,
                                                                                               aHardware,
                                                                                               aConfigurations,
                                                                                               *miss_locations,
                                                                                               *obs_locations, aKernel);
        } else {
            prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLEPredictTile(aData,
                                                                                    (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                    z_miss_number, n_z_obs, z_obs,
                                                                                    z_actual, z_miss, aHardware,
                                                                                    aConfigurations, *miss_locations,
                                                                                    *obs_locations, aKernel);
        }
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += prediction_error_mspe[i];
        }
        LOGGER("- MSPE: " << avg_pred_value[0])
        LOGGER("- Average prediction Error (mspe): ", true)
        for (i = 0; i < number_of_mspe; i++) {
            LOGGER_PRECISION(avg_pred_value[i] << "\t", 8)
        }
        LOGGER("")
        delete[] prediction_error_mspe;
    }

    delete[] z_obs;
    delete[] z_miss;
    delete[] z_actual;
    delete miss_locations;
    delete obs_locations;
}

template<typename T>
void Prediction<T>::InitializePredictionArguments(Configurations &aConfigurations,
                                                  std::unique_ptr<dataunits::ExaGeoStatData<T>> &aData,
                                                  std::unique_ptr<linearAlgebra::LinearAlgebraMethods<T>> &aLinearAlgebraSolver,
                                                  T *apZObs, T *apZActual, Locations<T> &aMissLocation,
                                                  Locations<T> &aObsLocation, T *apMeasurementsMatrix, const int &aP) {

    int full_problem_size = aConfigurations.GetProblemSize() * aP;
    T *z = new T[full_problem_size];

    aLinearAlgebraSolver->ExaGeoStatGetZObs(aConfigurations, z, full_problem_size, *aData->GetDescriptorData(),
                                            apMeasurementsMatrix, aP);
    PredictionHelpers<T>::PickRandomPoints(aConfigurations, aData, apZObs, apZActual, z, aMissLocation, aObsLocation,
                                           aP);
    delete[] z;
}