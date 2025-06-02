#!/usr/bin/env bash
# process_datasets.sh
# Public script to preprocess datasets into binary form for graph processing.

set -euo pipefail

PYTHON=python3
BUILD_TYPE=release
BUILD_DIR=build/${BUILD_TYPE}
EXTRA_ARGS=""
PROCESS_JOBS=8

CUSTOM=false
GRAMI=false
PROGRESS=false

GRAMI_URLS_PREFIX=https://raw.githubusercontent.com/ehab-abdelhamid/GraMi/refs/heads/master/Datasets/
GRAMI_DATASETS=(citeseer mico)
GRAMI_URLS_SUFFIX=.lg

# Parse flags
for arg in "$@"; do
  case "$arg" in
    --progress) PROGRESS=true ;;
    --custom) CUSTOM=true ;;
    --grami) GRAMI=true ;;
    *) echo "Unknown option: $arg" >&2; exit 1 ;;
  esac
done

if $PROGRESS; then
  EXTRA_ARGS+=" --progress"
fi

# Ensure build directory exists
if [ ! -d "$BUILD_DIR" ]; then
  echo "[INFO] Generating build directory with CMake preset: ${BUILD_TYPE}"
  cmake --preset "${BUILD_TYPE}"
else
  echo "[INFO] Found existing build directory: ${BUILD_DIR}"
fi

# Build the processing binary
echo "[INFO] Building processor..."
cmake --build --preset "${BUILD_TYPE}" -j "${PROCESS_JOBS}" --target process

# Ensure required directories exist
mkdir -p ./data/info ./data/extracted ./data/processed ./data/lg ./data/custom

# Extract .gz files
for gz_file in ./data/downloaded/*.gz; do
  [ -e "$gz_file" ] || continue
  base_name=$(basename "$gz_file" .gz)
  if [ ! -f "./data/extracted/$base_name" ]; then
    echo "[INFO] Extracting $gz_file"
    gzip -dkc < "$gz_file" > "./data/extracted/$base_name"
  fi
done

# Process extracted `.txt` files into binary
for extracted in ./data/extracted/*.txt; do
  [ -e "$extracted" ] || continue
  base_name=$(basename "$extracted" .txt)
  out_file="./data/processed/$base_name.bin"
  if [ ! -f "$out_file" ]; then
    echo "[INFO] Processing extracted: $base_name"
    ./"${BUILD_DIR}/bin/process" -i "$extracted" -o "$out_file" --format citPatentsTxt --seed 0 $EXTRA_ARGS
  fi
done

# If --grami, copy datasets from user path
if $GRAMI; then
  echo "[INFO] Downloading GraMi datasets..."
  for dataset in "${GRAMI_DATASETS[@]}"; do
    lg_file="./data/lg/${dataset}${GRAMI_URLS_SUFFIX}"
    if [ ! -f "$lg_file" ]; then
      echo "[INFO] Downloading GraMi dataset: $dataset"
      curl -L "${GRAMI_URLS_PREFIX}${dataset}${GRAMI_URLS_SUFFIX}" -o "$lg_file"
    fi
  done
fi

# Process `.lg` files
for lg_file in ./data/lg/*.lg; do
  [ -e "$lg_file" ] || continue
  base_name=$(basename "$lg_file" .lg)
  out_file="./data/processed/$base_name.bin"
  if [ ! -f "$out_file" ]; then
    echo "[INFO] Processing GraMi: $base_name"
    ./"${BUILD_DIR}/bin/process" -i "$lg_file" -o "$out_file" --format lg --seed 0 $EXTRA_ARGS
  fi
done

# Process custom datasets if --custom is given
for custom_file in ./data/custom/*.txt; do
  [ -e "$custom_file" ] || continue
  base_name=$(basename "$custom_file" .txt)
  out_file="./data/processed/$base_name.bin"
  if $CUSTOM || [ ! -f "$out_file" ]; then
    echo "[INFO] Processing custom: $base_name"
    ./"${BUILD_DIR}/bin/process" -i "$custom_file" -o "$out_file" --format custom --seed 0 $EXTRA_ARGS
  fi
done

# Generate LG dataset
echo "[INFO] Generating LG dataset (Python)..."
"$PYTHON" ./python/generateLGDataset.py --verbose

echo "[INFO] Dataset processing completed successfully."
