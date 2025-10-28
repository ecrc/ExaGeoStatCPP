#!/bin/bash
################################################################################
# Synchronized StageZero Processing
# Processes each latitude band independently but writes to SHARED output files
################################################################################

#==============================================================================
# CONFIGURATION - Edit these values for your setup OR override via command line
#==============================================================================

# Default values (can be overridden by command-line arguments)
# Stage Zero Configuration
total_latitudes=72
lon=144
startyear=2000
endyear=2002   
data_path="/path/to/ERA_data/"
forcing_data_path="/path/to/forcing_new.csv"
resultspath="/path/to/results/"

# Pipeline Configuration
STAGEZERO_BIN="./bin/examples/stage-zero/Example_Stage_Zero"
mpi_processes=2
parallel_jobs=10

# Climate Emulator Configuration
run_climate_emulator=false
EMULATOR_BIN="./bin/examples/climate-emulator/Example_Climate_Emulator"
N=5184
dts=72
timeslot=100
ObjectsNumber=72
add_diagonal=10
Accuracy=0
band_dense=1000
verbose="detailed"
gpus=0

# Parse command-line arguments to override defaults
while [[ $# -gt 0 ]]; do
    case $1 in
        --total-latitudes=*)
            total_latitudes="${1#*=}"
            shift
            ;;
        --lon=*)
            lon="${1#*=}"
            shift
            ;;
        --startyear=*)
            startyear="${1#*=}"
            shift
            ;;
        --endyear=*)
            endyear="${1#*=}"
            shift
            ;;
        --data-path=*)
            data_path="${1#*=}"
            shift
            ;;
        --forcing-data-path=*)
            forcing_data_path="${1#*=}"
            shift
            ;;
        --resultspath=*)
            resultspath="${1#*=}"
            shift
            ;;
        --mpi-processes=*)
            mpi_processes="${1#*=}"
            shift
            ;;
        --parallel-jobs=*)
            parallel_jobs="${1#*=}"
            shift
            ;;
        --cores=*)
            cores="${1#*=}"
            shift
            ;;
        --run-climate-emulator)
            run_climate_emulator=true
            shift
            ;;
        --N=*)
            N="${1#*=}"
            shift
            ;;
        --dts=*)
            dts="${1#*=}"
            shift
            ;;
        --timeslot=*)
            timeslot="${1#*=}"
            shift
            ;;
        --ObjectsNumber=*)
            ObjectsNumber="${1#*=}"
            shift
            ;;
        --cores=*)
            cores="${1#*=}"
            shift
            ;;
        --add-diagonal=*)
            add_diagonal="${1#*=}"
            shift
            ;;
        --Accuracy=*)
            Accuracy="${1#*=}"
            shift
            ;;
        --band-dense=*)
            band_dense="${1#*=}"
            shift
            ;;
        --verbose=*)
            verbose="${1#*=}"
            shift
            ;;
        --gpus=*)
            gpus="${1#*=}"
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "===== REQUIRED Stage Zero Options ====="
            echo "  --total-latitudes=N      Number of latitude bands to process"
            echo "  --lon=N                  Number of longitudes per latitude"
            echo "  --startyear=YYYY         Starting year for data processing"
            echo "  --endyear=YYYY           Ending year for data processing"
            echo "  --data-path=PATH         Path to NetCDF data directory"
            echo "  --forcing-data-path=PATH Path to forcing data file"
            echo "  --resultspath=PATH       Output directory path"
            echo ""
            echo "===== OPTIONAL Pipeline Options ====="
            echo "  --mpi-processes=N        MPI processes per job (default: $mpi_processes)"
            echo "  --parallel-jobs=N        Number of parallel jobs (default: $parallel_jobs)"
            echo "  --cores=N                Cores for Stage Zero and Emulator (default: auto-calculated)"
            echo ""
            echo "===== Climate Emulator Options ====="
            echo "  --run-climate-emulator   Run climate emulator after Stage Zero"
            echo "  --N=N                    Spatial problem size (default: $N)"
            echo "  --dts=N                  Dense tile size (default: $dts)"
            echo "  --timeslot=N             Number of z files to process (default: $timeslot)"
            echo "  --ObjectsNumber=N        Objects number (default: $ObjectsNumber)"
            echo "  --add-diagonal=N         Diagonal value (default: $add_diagonal)"
            echo "  --Accuracy=N             Accuracy (default: $Accuracy)"
            echo "  --band-dense=N           Band dense (default: $band_dense)"
            echo "  --verbose=MODE           Verbose mode (default: $verbose)"
            echo "  --gpus=N                 Number of GPUs (default: $gpus)"
            echo ""
            echo "===== Parameter Relationships ====="
            echo "Stage Zero Output:"
            echo "  - Creates: total_latitudes × lon spatial locations"
            echo "  - Generates: 365 × 24 × (endyear-startyear+1) z_*.csv files"
            echo "  - Each z file contains one timeslot for all spatial locations"
            echo "  - Processing time scales with: locations × years × optimization_iterations"
            echo ""
            echo "Climate Emulator Constraints:"
            echo "  - CRITICAL: N must equal dts² (e.g., N=5184 requires dts=72)"
            echo "  - N should be ≤ (total_latitudes × lon) from Stage Zero output"
            echo "  - timeslot: number of z files to read (max = Stage Zero output count)"
            echo ""
            echo "===== Examples ====="
            echo "  1. Stage Zero only (3 years of data):"
            echo "    $0 --total-latitudes=72 --lon=144 --startyear=2020 --endyear=2022 \\"
            echo "       --data-path=/path/to/ERA_data/ --forcing-data-path=/path/to/forcing.csv \\"
            echo "       --resultspath=/path/to/results/ --parallel-jobs=10 --cores=4"
            echo ""
            echo "  2. Full Pipeline (Stage Zero + Climate Emulator):"
            echo "    $0 --total-latitudes=72 --lon=144 --startyear=2020 --endyear=2022 \\"
            echo "       --data-path=/path/to/ERA_data/ --forcing-data-path=/path/to/forcing.csv \\"
            echo "       --resultspath=/path/to/results/ --parallel-jobs=10 --cores=4 \\"
            echo "       --run-climate-emulator --N=5184 --dts=72 --timeslot=3000"
            echo ""
            echo "  Note: In example 2, Stage Zero creates 26,280 z files (3 years × 365 × 24),"
            echo "        but emulator processes only first 3,000 for faster testing."
            echo "        For full processing: --timeslot=26280"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Calculate optimal cores per job to avoid oversubscription (if not manually set)
cpu_cores=$(nproc)
if [[ -z "$cores" ]]; then
    cores=$((cpu_cores / parallel_jobs))
    
    # Ensure at least 1 core per job, but not more than total cores
    if [ $cores -lt 1 ]; then
        cores=1
    elif [ $cores -gt $cpu_cores ]; then
        cores=$cpu_cores
    fi
fi

#==============================================================================
# YOU USUALLY DON'T NEED TO EDIT BELOW THIS LINE
#==============================================================================

echo ""
echo "=========================================="
echo "  StageZero Pipeline"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Years:        ${startyear}-${endyear} ($((endyear - startyear + 1)) years)"
echo "  Grid:         ${total_latitudes} latitudes × ${lon} longitudes"
echo "  Data path:    ${data_path}"
echo "  Forcing:      ${forcing_data_path}"
echo "  Results:      ${resultspath}"
echo "  Binary:       ${STAGEZERO_BIN}"
echo ""
echo "Parallelization:"
echo "  System cores: ${cpu_cores}"
echo "  Parallel jobs: ${parallel_jobs} (1 latitude per job)"
echo "  MPI ranks/job: ${mpi_processes}"
echo "  Cores/job:     ${cores}"
echo ""
years=$((endyear - startyear + 1))
time_slots=$((365 * 24 * years))
echo "Expected output: ${time_slots} timeslot files (z_*.csv) + 1 params.csv"
echo "=========================================="
echo ""

# Check for Climate Emulator auxiliary files if emulator is enabled
if [ "$run_climate_emulator" = true ]; then
    echo "⚠ Climate Emulator Check:"
    echo "  The Climate Emulator requires pre-computed auxiliary files in the results directory."
    echo "  Required files (example for dts=720): 720_Et1.csv, 720_Et2.csv, 720_Ep.csv, etc."
    echo "  Please copy auxiliary files to the results directory before running."
    echo ""
fi

# Create shared output directory
mkdir -p "${resultspath}"

# Clear existing output files for clean start (optional - comment out if you want to append to existing)
echo "Clearing existing output files for clean start..."
rm -f "${resultspath}"/z_*.csv
rm -f "${resultspath}"/params.csv
rm -f "${resultspath}"/lat_*_output.log
rm -f "${resultspath}"/climate_emulator.log
echo "✓ Output files and logs cleared"
echo "Files will be created on-demand by the C++ processes"
echo ""

# Create a lock directory for coordination
lock_dir="${resultspath}/.locks"
rm -rf "${lock_dir}"
mkdir -p "${lock_dir}"

# Counter for active background jobs
active_jobs=0

# Process each latitude
for lat in $(seq 0 $((total_latitudes - 1))); do
    
    echo "[$(date +%H:%M:%S)] Starting latitude ${lat}..."
    
    # Run this latitude in background with environment variable set
    (
        # Set environment variable for the C++ program to know which latitude band this is
        export STAGEZERO_LATITUDE_BAND=${lat}
        
        mpirun -n ${mpi_processes} ${STAGEZERO_BIN} \
            --kernel=trend_model \
            --lon=${lon} \
            --lat=${lat} \
            --dts=720 \
            --startyear=${startyear} \
            --endyear=${endyear} \
            --data-path=${data_path} \
            --forcing-data-path=${forcing_data_path} \
            --resultspath=${resultspath} \
            --starting-theta=0.9 \
            --lb=0.001 \
            --ub=0.95 \
            --tolerance=7 \
            --max-mle-iterations=30 \
            --stage-zero \
            --cores=${cores} \
            --gpus=0 \
            --p=1 \
            --q=1 \
            > "${resultspath}/lat_${lat}_output.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "[$(date +%H:%M:%S)] ✓ Latitude ${lat} complete"
            touch "${lock_dir}/lat_${lat}.done"
        else
            echo "[$(date +%H:%M:%S)] ✗ Latitude ${lat} FAILED"
            touch "${lock_dir}/lat_${lat}.failed"
        fi
    ) &
    
    # Increase counter
    ((active_jobs++))
    
    # If we've reached the limit, wait for one to finish
    if [ ${active_jobs} -ge ${parallel_jobs} ]; then
        wait -n  # Wait for any one job to finish
        ((active_jobs--))
    fi
done

# Wait for all remaining jobs to finish
echo ""
echo "[$(date +%H:%M:%S)] Waiting for remaining jobs to complete..."
wait

echo ""
echo "=========================================="
echo "  Stage Zero Processing Complete"
echo "=========================================="

# Show summary
completed=$(find "${lock_dir}" -name "*.done" | wc -l)
failed=$(find "${lock_dir}" -name "*.failed" | wc -l)

echo ""
echo "Job Summary:"
echo "  ✓ Completed: ${completed}/${total_latitudes} latitudes"
if [ ${failed} -gt 0 ]; then
    echo "  ✗ Failed:    ${failed} latitude(s)"
fi

# Check output file sizes
echo ""
echo "Output Verification:"
if [ -f "${resultspath}/params.csv" ]; then
    param_lines=$(wc -l < "${resultspath}/params.csv")
    expected_param_lines=$((total_latitudes * lon))
    if [ ${param_lines} -eq ${expected_param_lines} ]; then
        echo "  ✓ params.csv:  ${param_lines}/${expected_param_lines} lines"
    else
        echo "  ⚠ params.csv:  ${param_lines}/${expected_param_lines} lines (MISMATCH)"
    fi
else
    echo "  ✗ params.csv:  FILE NOT FOUND"
fi

z_files=$(ls "${resultspath}"/z_*.csv 2>/dev/null | wc -l)
expected_z_lines=$((total_latitudes * lon))
if [ ${z_files} -gt 0 ]; then
    first_z_file=$(ls "${resultspath}"/z_*.csv | head -1)
    z_lines=$(wc -l < "${first_z_file}")
    if [ ${z_lines} -eq ${expected_z_lines} ]; then
        echo "  ✓ Z files:     ${z_files} files (${z_lines}/${expected_z_lines} lines each)"
    else
        echo "  ⚠ Z files:     ${z_files} files (${z_lines}/${expected_z_lines} lines each - MISMATCH)"
    fi
else
    echo "  ✗ Z files:     NO FILES FOUND"
fi

echo ""
echo "Results location: ${resultspath}"

if [ ${completed} -lt ${total_latitudes} ]; then
    echo ""
    echo "⚠ Failed Jobs - Check Logs:"
    for lat in $(seq 0 $((total_latitudes - 1))); do
        if [ -f "${lock_dir}/lat_${lat}.failed" ]; then
            echo "  - Latitude ${lat}: ${resultspath}/lat_${lat}_output.log"
        fi
    done
    echo "=========================================="
    echo ""
    exit 1
fi

echo "=========================================="
echo ""

# Run Climate Emulator if requested
if [ "$run_climate_emulator" = true ]; then
    echo ""
    echo "=========================================="
    echo "  Climate Emulator"
    echo "=========================================="
    echo ""
    
    if [ ! -f "$EMULATOR_BIN" ]; then
        echo "ERROR: Climate Emulator binary not found: $EMULATOR_BIN"
        exit 1
    fi
    
    emulator_log="${resultspath}/climate_emulator.log"
    
    echo "Configuration:"
    echo "  Binary:       $EMULATOR_BIN"
    echo "  Problem size: $N"
    echo "  Tile size:    $dts"
    echo "  Timeslots:    $timeslot"
    echo "  Objects:      $ObjectsNumber"
    echo "  Cores:        $cores"
    echo "  Input data:   $resultspath"
    echo "  Log file:     $emulator_log"
    echo ""
    
    echo "Running Climate Emulator (output logged to $emulator_log)..."
    
    "$EMULATOR_BIN" \
        --N="$N" \
        --dts="$dts" \
        --timeslot="$timeslot" \
        --ObjectsNumber="$ObjectsNumber" \
        --cores="$cores" \
        --add-diagonal="$add_diagonal" \
        --Accuracy="$Accuracy" \
        --band_dense="$band_dense" \
        --verbose="$verbose" \
        --gpus="$gpus" \
        --data_path="$resultspath" \
        > "${emulator_log}" 2>&1
    
    EMULATOR_EXIT=$?
    
    echo ""
    if [ $EMULATOR_EXIT -eq 0 ]; then
        echo "✓ Climate Emulator completed successfully!"
        echo "  Log: $emulator_log"
        echo "=========================================="
        echo ""
    else
        echo "✗ Climate Emulator failed with exit code: $EMULATOR_EXIT"
        echo "  Check log: $emulator_log"
        echo "=========================================="
        echo ""
        exit $EMULATOR_EXIT
    fi
else
    echo "✓ Pipeline complete!"
    echo ""
fi
