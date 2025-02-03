#!/bin/bash

# Function to display help message
function show_help() {
    echo "Usage: $0 [options]"
    echo "Run a pipeline for coverage prediction, real coverage computation, and HMM segmentation."
    echo ""
    echo "Options:"
    echo "  --keys_path PATH               Path to the keys CSV file for the model."
    echo "  --hdf5_path PATH               Path to the HDF5 file for the model."
    echo "  --experiment_config_path PATH  Path to the experiment config YAML file for the model."
    echo "  --config_path PATH             Path to the model config JSON file."
    echo "  --checkpoint_path PATH         Path to the model checkpoint file."
    echo "  --model_output_file PATH       Output file for the model predictions (default: model_predictions.tsv)."
    echo "  --batch_size SIZE              Batch size for the model inference (default: 32)."
    echo "  --labels NUM                   Number of samples used to finetune the model."
    echo "  --sample_file PATH             Path to the file with sample names and bigWig file paths."
    echo "  --regions_file PATH            Path to the file with genomic regions in CSV format."
    echo "  --real_coverage_output PATH    Output file for real coverage (default: real_coverage.csv)."
    echo "  --hmm_output_file PATH         Output file for HMM predictions (default: hmm_predictions.tsv)."
    echo "  --help                         Show this help message and exit."
    echo ""
    echo "Example:"
    echo "  $0 --keys_path keys.csv --hdf5_path data.hdf5 --experiment_config_path config.yaml \\"
    echo "     --config_path model_config.json --labels num_samples --checkpoint_path checkpoint.pt --sample_file samples.txt \\"
    echo "     --regions_file regions.csv"
    exit 0
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --keys_path) KEYS_PATH="$2"; shift ;;
        --hdf5_path) HDF5_PATH="$2"; shift ;;
        --experiment_config_path) EXPERIMENT_CONFIG_PATH="$2"; shift ;;
        --config_path) CONFIG_PATH="$2"; shift ;;
        --checkpoint_path) CHECKPOINT_PATH="$2"; shift ;;
        --model_output_file) MODEL_OUTPUT_FILE="$2"; shift ;;
        --batch_size) BATCH_SIZE="$2"; shift ;;
        --labels) LABELS="$2"; shift ;;
        --sample_file) SAMPLE_FILE="$2"; shift ;;
        --regions_file) REGIONS_FILE="$2"; shift ;;
        --real_coverage_output) REAL_COVERAGE_OUTPUT="$2"; shift ;;
        --hmm_output_file) HMM_OUTPUT_FILE="$2"; shift ;;
        --help) show_help ;;
        *) echo "Unknown parameter: $1"; show_help ;;
    esac
    shift
done

# Set default values if not provided
MODEL_OUTPUT_FILE=${MODEL_OUTPUT_FILE:-"model_predictions.tsv"}
BATCH_SIZE=${BATCH_SIZE:-32}
REAL_COVERAGE_OUTPUT=${REAL_COVERAGE_OUTPUT:-"real_coverage.csv"}
HMM_OUTPUT_FILE=${HMM_OUTPUT_FILE:-"hmm_predictions.tsv"}

# Validate required arguments
if [[ -z "$KEYS_PATH" || -z "$HDF5_PATH" || -z "$EXPERIMENT_CONFIG_PATH" || -z "$CONFIG_PATH" || -z "$CHECKPOINT_PATH" || -z "$SAMPLE_FILE" || -z "$REGIONS_FILE" || -z "$LABELS" ]]; then
    echo "Error: Missing required arguments."
    show_help
fi

# Step 1: Run the coverage prediction model
echo "Running predictions_by_the_model.py..."
python predictions_by_the_model.py \
    --keys_path "$KEYS_PATH" \
    --hdf5_path "$HDF5_PATH" \
    --experiment_config_path "$EXPERIMENT_CONFIG_PATH" \
    --config_path "$CONFIG_PATH" \
    --checkpoint_path "$CHECKPOINT_PATH" \
    --output_file "$MODEL_OUTPUT_FILE" \
    --batch_size "$BATCH_SIZE" \
    --labels "$LABELS"

# Check if the model predictions were generated
if [ ! -f "$MODEL_OUTPUT_FILE" ]; then
    echo "Error: Model predictions file not found. Exiting."
    exit 1
fi

# Step 2: Compute real coverage
echo "Running real_coverage_count.py..."
python real_coverage_count.py \
    --samples "$SAMPLE_FILE" \
    --regions "$REGIONS_FILE" \
    --output "$REAL_COVERAGE_OUTPUT"

# Check if the real coverage file was generated
if [ ! -f "$REAL_COVERAGE_OUTPUT" ]; then
    echo "Error: Real coverage file not found. Exiting."
    exit 1
fi

# Step 3: Run the HMM segmentator
echo "Running hmm_segmentator.py..."
python hmm_segmentator.py \
    --prediction_coverage "$MODEL_OUTPUT_FILE" \
    --real_coverage "$REAL_COVERAGE_OUTPUT" \
    --output_file "$HMM_OUTPUT_FILE"

# Check if the HMM predictions were generated
if [ ! -f "$HMM_OUTPUT_FILE" ]; then
    echo "Error: HMM predictions file not found. Exiting."
    exit 1
fi

echo "Pipeline completed successfully. Output saved to $HMM_OUTPUT_FILE."