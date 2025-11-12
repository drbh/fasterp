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

# Function to create merge test files (R2 is reverse complement of R1)
create_merge_test_files() {
    local num_reads=$1
    local r1_file=$2
    local r2_file=$3
    local read_length=${4:-60}  # Default 60bp for perfect overlaps

    echo ""
    echo "Creating merge test files: $r1_file and $r2_file with $num_reads reads..."

    # Create R1 first
    > "$r1_file"
    > "$r2_file"

    # Function to reverse complement a DNA sequence
    reverse_complement() {
        echo "$1" | rev | tr 'ACGTacgt' 'TGCAtgca'
    }

    for i in $(seq 0 $((num_reads - 1))); do
        # Generate random sequence for R1
        local seq=""
        local bases="ACGT"
        for ((j=0; j<read_length; j++)); do
            seq="${seq}${bases:$((RANDOM % 4)):1}"
        done

        # Create quality string (all high quality)
        local qual=$(printf 'I%.0s' $(seq 1 $read_length))

        # Write R1
        echo "@read${i}" >> "$r1_file"
        echo "$seq" >> "$r1_file"
        echo "+" >> "$r1_file"
        echo "$qual" >> "$r1_file"

        # Write R2 as reverse complement of R1
        local r2_seq=$(reverse_complement "$seq")
        echo "@read${i}" >> "$r2_file"
        echo "$r2_seq" >> "$r2_file"
        echo "+" >> "$r2_file"
        echo "$qual" >> "$r2_file"
    done

    echo "  Created merge test files with $num_reads read pairs ($read_length bp each)"
}

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

# Generate merge test files (R2 is reverse complement of R1 for perfect merging)
echo ""
echo "Generating merge test files for read merging tests..."
create_merge_test_files 3 "merge_test_R1.fq" "merge_test_R2.fq" 60
create_merge_test_files 1000 "merge_test_10k_R1.fq" "merge_test_10k_R2.fq" 100

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

# Function to create adapter trimming test file
create_adapter_test_file() {
    echo ""
    echo "Generating adapter_trim_test.fastq with embedded adapters..."

    local output="adapter_trim_test.fastq"
    local adapter="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    local total_reads=1000
    local reads_with_adapters=50

    # Function to generate random DNA sequence
    random_dna() {
        local length=$1
        local seq=""
        local bases="ACGT"
        for ((i=0; i<length; i++)); do
            seq="${seq}${bases:$((RANDOM % 4)):1}"
        done
        echo "$seq"
    }

    # Clear output file
    > "$output"

    for i in $(seq 1 $total_reads); do
        # Read header
        echo "@read_${i}" >> "$output"

        # Decide if this read should have an adapter
        if [ $i -le $reads_with_adapters ]; then
            # Generate read with adapter at various positions
            insert_pos=$((50 + RANDOM % 80))  # Random position between 50-130
            prefix=$(random_dna $insert_pos)

            # Calculate how much adapter to include (partial or full)
            remaining=$((150 - insert_pos))
            adapter_len=${#adapter}
            if [ $remaining -lt $adapter_len ]; then
                # Partial adapter at the end
                adapter_part="${adapter:0:$remaining}"
                seq="${prefix}${adapter_part}"
            else
                # Full adapter + some random bases
                suffix_len=$((150 - insert_pos - adapter_len))
                suffix=$(random_dna $suffix_len)
                seq="${prefix}${adapter}${suffix}"
            fi
        else
            # Regular read without adapter
            seq=$(random_dna 150)
        fi

        echo "$seq" >> "$output"
        echo "+" >> "$output"

        # Quality string (all high quality)
        qual=$(printf 'I%.0s' $(seq 1 ${#seq}))
        echo "$qual" >> "$output"
    done

    echo "  Created $output with $total_reads reads ($reads_with_adapters with adapters)"
}

# Create adapter trimming test file
create_adapter_test_file

echo ""
echo "All test data files created successfully!"
echo ""
echo "Single-end files:"
ls -lh *.fq | grep -v "_R[12].fq"
echo ""
echo "Merge test files:"
ls -lh merge_test_*.fq
echo ""
echo "Paired-end overlap files:"
ls -lh overlap_*.fq
echo ""
echo "Standard paired-end files:"
ls -lh pe_*.fq
echo ""
echo "Adapter test file:"
ls -lh adapter_trim_test.fastq
