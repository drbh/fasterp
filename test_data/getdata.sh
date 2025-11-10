#!/bin/bash

# Generate test data using nucgen
echo "Generating test data using nucgen..."

# Check if nucgen is available
if ! command -v nucgen &> /dev/null; then
    echo "Error: nucgen is not installed or not in PATH"
    echo "Please install nucgen: cargo install nucgen"
    exit 1
fi

# Configuration
SEED=42
READ_LENGTH=150

# Function to create a test file with N reads
create_test_file() {
    local num_reads=$1
    local output_file=$2

    echo "Creating $output_file with $num_reads reads..."

    nucgen -n "$num_reads" -l "$READ_LENGTH" -f q -S "$SEED" "$output_file"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate $output_file"
        exit 1
    fi

    local actual_reads=$(wc -l < "$output_file")
    actual_reads=$((actual_reads / 4))
    echo "  Created with $actual_reads reads"
}

# Download R1.fq from fastp for compatibility tests
echo ""
echo "Downloading R1.fq from fastp repository for compatibility..."
if command -v wget &> /dev/null; then
    wget -q https://raw.githubusercontent.com/OpenGene/fastp/refs/heads/master/testdata/R1.fq
elif command -v curl &> /dev/null; then
    curl -sL https://raw.githubusercontent.com/OpenGene/fastp/refs/heads/master/testdata/R1.fq -o R1.fq
else
    echo "Warning: Neither wget nor curl found. Skipping R1.fq download."
fi

if [ -f R1.fq ]; then
    count=$(wc -l < R1.fq)
    count=$((count / 4))
    echo "  R1.fq downloaded with $count reads"
fi

# Create test files of different sizes
echo ""
echo "Generating synthetic test files..."
create_test_file 1000 "small_1k.fq"
create_test_file 10000 "medium_10k.fq"
create_test_file 100000 "large_100k.fq"
create_test_file 500000 "xlarge_500k.fq"

echo ""
echo "Test data files created successfully:"
ls -lh *.fq

echo ""
echo "Read counts:"
for f in R1.fq small_1k.fq medium_10k.fq large_100k.fq xlarge_500k.fq; do
    if [ -f "$f" ]; then
        count=$(wc -l < "$f")
        count=$((count / 4))
        echo "  $f: $count reads"
    fi
done

# Function to create paired-end test files with overlap
create_pe_overlap_files() {
    local num_reads=$1
    local r1_file=$2
    local r2_file=$3

    echo ""
    echo "Creating PE overlap files: $r1_file and $r2_file with $num_reads reads..."

    # Create R1 with shorter read length for overlap (100bp)
    nucgen -n "$num_reads" -l 100 -f q -S "$SEED" "$r1_file"

    # Create R2 with same seed + offset for overlap
    # R2 should be reverse complement and overlap with R1
    nucgen -n "$num_reads" -l 100 -f q -S $((SEED + 1)) "$r2_file"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate PE overlap files"
        exit 1
    fi

    echo "  Created PE overlap files with $num_reads read pairs"
}

# Function to create standard paired-end test files
create_pe_files() {
    local num_reads=$1
    local r1_file=$2
    local r2_file=$3

    echo ""
    echo "Creating PE files: $r1_file and $r2_file with $num_reads reads..."

    # Create R1
    nucgen -n "$num_reads" -l "$READ_LENGTH" -f q -S "$SEED" "$r1_file"

    # Create R2 with different seed
    nucgen -n "$num_reads" -l "$READ_LENGTH" -f q -S $((SEED + 100)) "$r2_file"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate PE files"
        exit 1
    fi

    echo "  Created PE files with $num_reads read pairs"
}

# Generate overlap test files (shorter reads that can overlap)
echo ""
echo "Generating paired-end overlap test files..."
create_pe_overlap_files 1000 "overlap_R1.fq" "overlap_R2.fq"
create_pe_overlap_files 10000 "overlap_10k_R1.fq" "overlap_10k_R2.fq"

# Generate standard paired-end test files
echo ""
echo "Generating standard paired-end test files..."
create_pe_files 1000 "pe_small_R1.fq" "pe_small_R2.fq"
create_pe_files 10000 "pe_medium_R1.fq" "pe_medium_R2.fq"
create_pe_files 10000 "pe_10k_R1.fq" "pe_10k_R2.fq"
create_pe_files 100000 "pe_large_R1.fq" "pe_large_R2.fq"

echo ""
echo "All test data files created successfully!"
echo ""
echo "Single-end files:"
ls -lh *.fq | grep -v "_R[12].fq"
echo ""
echo "Paired-end overlap files:"
ls -lh overlap_*.fq
echo ""
echo "Standard paired-end files:"
ls -lh pe_*.fq
