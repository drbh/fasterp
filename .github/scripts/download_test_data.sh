#!/bin/bash
set -euo pipefail

# Download real FASTQ files from the genedata.dholtz.com mirror and check them
# against the SHA-256 sums below.
#
#   bash .github/scripts/download_test_data.sh [dir] [set ...]
#
# With no sets, only the SRR5808766 pair used in CI is fetched. The other sets
# are benchmark data, described in book/src/benchmarks.md:
#
#   srr1766443-2m       human WGS, 2x148 bp, first 2M pairs, gzip -4       0.5 GB
#   srr1766443-2m-bgzf  the same records compressed with bgzip -l 4        0.5 GB
#   err3534523-2m       NovaSeq 6000 WGS, 2x150 bp, first 2M pairs         0.3 GB
#   err3534523          the full ERR3534523 run, 33.1M pairs, gzip -4      5.5 GB
#   srr891269-2m        GM12878 ATAC-seq, 2x50 bp, first 2M pairs          0.2 GB
#   all                 every set above
#
# Files already present with the right checksum are skipped.

DATA_DIR="${1:-test_data}"
shift || true
SETS=("${@:-ci}")
BASE_URL="https://genedata.dholtz.com"

files_for() {
    case "$1" in
        ci)
            echo "dceb9047c28c81d5a8eef45d1618bc4e0e10e936212e7b67099c63302045cfe8 SRX2987343/SRR5808766_1.fastq.gz"
            echo "fd59538e0470a1910fadf633399340456f1ca5baab736bb6d0b6ae0ed48866b6 SRX2987343/SRR5808766_2.fastq.gz"
            ;;
        srr1766443-2m)
            echo "443aca7a7ac499d39cabef7a74f4b30107fa11ab32962ecfc2fe8d16515b77d8 SRX847862/SRR1766443_2M_1.fastq.gz"
            echo "1fc21353380cf250c45bb77f92b096484b56091cf70cfe8cc06695bd059c7528 SRX847862/SRR1766443_2M_2.fastq.gz"
            ;;
        srr1766443-2m-bgzf)
            echo "78146eee1f03dc7ffb371c7d5c8a54e87d4486e6ac590d1fc843215f8301db16 SRX847862/SRR1766443_2M_bgzf_1.fastq.gz"
            echo "57be405d9dd82a1f0e52e73fbba0fdea089878ef39795c66c1daba06b736951c SRX847862/SRR1766443_2M_bgzf_2.fastq.gz"
            ;;
        err3534523-2m)
            echo "e4b84c6e0b2729970d6915036c8d6543a58d57ca1a15010a469e6df2d054322c ERX3269667/ERR3534523_2M_1.fastq.gz"
            echo "28ba0e8594de9ffa0cd1f18860e783901ceac22f9be14a86e1a2cc4c8e011c86 ERX3269667/ERR3534523_2M_2.fastq.gz"
            ;;
        err3534523)
            echo "048001e206687d1049d3285b8e7fa6d0b9971943e20550c8f25c419440413151 ERX3269667/ERR3534523_1.fastq.gz"
            echo "f016a2cad6f57b9c26bd8c53f0e1b7824b25252d645e9c329a3f6dd83512d9fe ERX3269667/ERR3534523_2.fastq.gz"
            ;;
        srr891269-2m)
            echo "10691cc21fe6ea6a754b589824ccf9108daeea573bc3e8fd2b9873bffab45338 SRX298001/SRR891269_2M_1.fastq.gz"
            echo "99327587f5ba557960c0dfaf0bb89b95ad81e2184d75841f56c591dcaaf6d7bd SRX298001/SRR891269_2M_2.fastq.gz"
            ;;
        all)
            for set in ci srr1766443-2m srr1766443-2m-bgzf err3534523-2m err3534523 srr891269-2m; do
                files_for "$set"
            done
            ;;
        *)
            echo "Unknown set: $1" >&2
            return 1
            ;;
    esac
}

sha256() {
    if command -v sha256sum >/dev/null; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

mkdir -p "$DATA_DIR"
echo "Downloading test data to $DATA_DIR..."

for set in "${SETS[@]}"; do
    files_for "$set" | while read -r sum path; do
        file="$DATA_DIR/$(basename "$path")"
        if [ -f "$file" ] && [ "$(sha256 "$file")" = "$sum" ]; then
            echo "  $(basename "$path") already exists, skipping"
            continue
        fi
        echo "  Downloading $path..."
        curl -fL --retry 3 -o "$file.part" "$BASE_URL/$path"
        if [ "$(sha256 "$file.part")" != "$sum" ]; then
            echo "Checksum mismatch for $path" >&2
            rm -f "$file.part"
            exit 1
        fi
        mv "$file.part" "$file"
    done
done

echo "Download complete. Files:"
ls -lh "$DATA_DIR"/*.fastq.gz 2>/dev/null || echo "  No .fastq.gz files found"
