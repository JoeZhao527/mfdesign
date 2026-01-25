#!/bin/bash
# Script to extract MFDesign dataset from compressed archives
# Usage: ./scripts/mfdesign/extract_data.sh [data_dir]

set -e

DATA_DIR="${1:-datasets/mfdesign}"

echo "=== MFDesign Data Extraction ==="
echo "Data directory: $DATA_DIR"
echo ""

cd "$DATA_DIR"

# Check if zstd is available
if ! command -v zstd &> /dev/null; then
    echo "zstd not found, trying to install..."
    # Try to install via conda
    if command -v conda &> /dev/null; then
        conda install -y zstd
    else
        echo "Error: zstd not installed and cannot auto-install"
        echo "Please install zstd: sudo apt install zstd (Ubuntu) or conda install zstd"
        exit 1
    fi
fi

# Check if pv is available for progress display
USE_PV=false
if command -v pv &> /dev/null; then
    USE_PV=true
else
    echo "Note: pv not installed, progress bar disabled. Install with: sudo apt install pv"
fi

# Helper function to extract .tar.zst with optional progress bar
extract_zst() {
    local archive="$1"
    local size=$(stat -c%s "$archive" 2>/dev/null || stat -f%z "$archive" 2>/dev/null)
    if [ "$USE_PV" = true ]; then
        pv -s "$size" "$archive" | zstd -d | tar -xf -
    else
        zstd -d -c "$archive" | tar -xf -
    fi
}

# Extract raw_data.tar.zst -> raw_data/ (原始 PDB 数据)
if [ -f "raw_data.tar.zst" ] && [ ! -d "raw_data" ]; then
    echo "Extracting raw_data.tar.zst..."
    extract_zst raw_data.tar.zst
    echo "  Done: raw_data/"
else
    if [ -d "raw_data" ]; then
        echo "raw_data/ already exists, skipping extraction"
    else
        echo "Warning: raw_data.tar.zst not found"
    fi
fi

# Extract antibody_data.tar.zst -> antibody_data/ (处理后的 JSON records)
if [ -f "antibody_data.tar.zst" ] && [ ! -d "antibody_data" ]; then
    echo "Extracting antibody_data.tar.zst..."
    extract_zst antibody_data.tar.zst
    echo "  Done: antibody_data/"
else
    if [ -d "antibody_data" ]; then
        echo "antibody_data/ already exists, skipping extraction"
    else
        echo "Warning: antibody_data.tar.zst not found"
    fi
fi

# Extract msa.tar.zst
if [ -f "msa.tar.zst" ] && [ ! -d "msa" ]; then
    echo "Extracting msa.tar.zst..."
    extract_zst msa.tar.zst
    echo "  Done: msa/"
else
    if [ -d "msa" ]; then
        echo "msa/ already exists, skipping extraction"
    else
        echo "Warning: msa.tar.zst not found"
    fi
fi

# Extract test_yaml_dir.tar.zst
if [ -f "test_yaml_dir.tar.zst" ] && [ ! -d "test_yaml_dir" ]; then
    echo "Extracting test_yaml_dir.tar.zst..."
    extract_zst test_yaml_dir.tar.zst
    echo "  Done: test_yaml_dir/"
else
    if [ -d "test_yaml_dir" ]; then
        echo "test_yaml_dir/ already exists, skipping extraction"
    else
        echo "Warning: test_yaml_dir.tar.zst not found"
    fi
fi

echo ""
echo "=== Extraction Complete ==="
echo ""
echo "Directory contents:"
ls -la

# Check extracted directories
for dir in raw_data antibody_data msa test_yaml_dir; do
    if [ -d "$dir" ]; then
        count=$(find "$dir" -type f 2>/dev/null | wc -l)
        echo "  $dir/: $count files"
    fi
done
