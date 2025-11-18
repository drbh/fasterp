#!/bin/bash
set -e

# Download real FASTQ test files for CI testing
# Files are cached by GitHub Actions to avoid repeated downloads

DATA_DIR="${1:-test_data}"
mkdir -p "$DATA_DIR"

BASE_URL="https://genedata.dholtz.com/SRX2987343"
FILES=(
    "SRR5808766_1.fastq.gz"
    "SRR5808766_2.fastq.gz"
)

echo "Downloading test data to $DATA_DIR..."

for file in "${FILES[@]}"; do
    if [ -f "$DATA_DIR/$file" ]; then
        echo "  $file already exists, skipping"
    else
        echo "  Downloading $file..."
        curl -L -o "$DATA_DIR/$file" "$BASE_URL/$file"
    fi
done

echo "Download complete. Files:"
ls -lh "$DATA_DIR"/*.fastq.gz 2>/dev/null || echo "  No .fastq.gz files found"
