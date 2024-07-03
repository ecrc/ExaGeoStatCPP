#!/bin/bash
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file Benchmarking.sh
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-10-10

# Set the output file name
output_file="40CORE_EXACT.csv"

# Create or truncate the output file
echo "N,computation,kernel,ncores,ngpus,dts,Zmiss,maximum_likelihood,maximum_theta,mspe_value,time_generation,time_modeling,time_prediction,time_modeling_per_iter,Gflops_data_generation,Gflops_data_modeling,Gflops_data_modeling_iter,Gflops_data_mspe" > "$output_file"

CORES=40
GPUS=0
MLE_ITERATIONS=1
DTS=960

# Define the desired values for N
desired_N=(55225 63001 71289 79524 87616 96100 104329 112225 120081 130889 150000 200000)

# Loop over N values
for N in "${desired_N[@]}"; do
  # Repeat each iteration 10 times
    for iteration in {1..10}; do
      # Run the command and capture the output
      command_output=$(./bin/examples/end-to-end/Example_Data_Generation_Modeling_and_Prediction --ncores=$CORES --gpus=$GPUS --computation=exact --itheta=1:0.1:0.5 --etheta=1:0.1:\? --olb=0.1:0.1:0.1 --oub=5:5:5 --dts="$DTS" --verbose=detailed --N="$N" --max_mle_iterations=$MLE_ITERATIONS --kernel=univariate_matern_stationary --tolerance=4 --Zmiss="$(($N/10))" --mspe)
      # Extract the desired values from the command output
      COMPUTATION=$(echo "$command_output" | awk -F "#Computation:" '{print $2}' | awk -F "," '{print $1}')
      COMPUTATION=$(echo "$COMPUTATION" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      KERNEL=$(echo "$command_output" | awk -F "#Kernel:" '{print $2}' | awk -F "," '{print $1}')
      KERNEL=$(echo "$KERNEL" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      logli_result=$(echo "$command_output" | awk -F "#Final Log Likelihood value: " '{print $2}' | awk -F "," '{print $1}')
      logli_result=$(echo "$logli_result" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      maximum_theta=$(echo "$command_output" | awk -F "#Found Maximum Theta at: " '{print $2}' | awk -F "," '{print $1}')
      maximum_theta=$(echo "$maximum_theta" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      mspe_value=$(echo "$command_output" | awk -F "#Mean Square Error MSPE: " '{print $2}' | awk -F "," '{print $1}')
      mspe_value=$(echo "$mspe_value" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      time_generation=$(echo "$command_output" | awk -F "#Total Data Generation Execution Time: " '{print $2}' | awk -F "," '{print $1}')
      time_generation=$(echo "$time_generation" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      time_modeling=$(echo "$command_output" | awk -F "#Total MLE Execution time: " '{print $2}' | awk -F "," '{print $1}')
      time_modeling=$(echo "$time_modeling" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      time_prediction=$(echo "$command_output" | awk -F "#MSPE Prediction Execution Time: " '{print $2}' | awk -F "," '{print $1}')
      time_prediction=$(echo "$time_prediction" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      time_modeling_iteration=$(echo "$command_output" | awk -F "#Average Time Modeling per Iteration: " '{print $2}' | awk -F "," '{print $1}')
      time_modeling_iteration=$(echo "$time_modeling_iteration" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      flops_generation=$(echo "$command_output" | awk -F "#Total Data Generation GFlop/s: " '{print $2}' | awk -F "," '{print $1}')
      flops_generation=$(echo "$flops_generation" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      flops_modeling=$(echo "$command_output" | awk -F "#Total MLE GFlop/s: " '{print $2}' | awk -F "," '{print $1}')
      flops_modeling=$(echo "$flops_modeling" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      flops_modeling_iteration=$(echo "$command_output" | awk -F "#Average Flops per Iteration: " '{print $2}' | awk -F "," '{print $1}')
      flops_modeling_iteration=$(echo "$flops_modeling_iteration" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      flops_mspe=$(echo "$command_output" | awk -F "#MSPE GFlop/s: " '{print $2}' | awk -F "," '{print $1}')
      flops_mspe=$(echo "$flops_mspe" | tr -d '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')

    echo "$command_output"
    echo "Iteration number $iteration done."
      # Append the values to the output file
      echo '**Results**'
      echo "$N,$COMPUTATION,$KERNEL,$CORES,$GPUS,$DTS,$((N/10)),$logli_result,$maximum_theta,$mspe_value,$time_generation,$time_modeling,$time_prediction,$time_modeling_iteration,$flops_generation,$flops_modeling,$flops_modeling_iteration,$flops_mspe"
      echo ''
      echo "$N,$COMPUTATION,$KERNEL,$CORES,$GPUS,$DTS,$((N/10)),$logli_result,$maximum_theta,$mspe_value,$time_generation,$time_modeling,$time_prediction,$time_modeling_iteration,$flops_generation,$flops_modeling,$flops_modeling_iteration,$flops_mspe" >> "$output_file"
    done
done
