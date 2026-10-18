#!/usr/bin/env bash
# Build the exact compatibility reference using installed ISA-L, libdeflate
# and Highway >= 1.1.0 development libraries. Never reuse an existing checkout.
set -euo pipefail

destination=${1:?Usage: bash scripts/build_fastp_reference.sh NEW_DIRECTORY}
reference=ee125a52fddae9d8189e13cd8a33ad04988a1351
if [[ -e "$destination" ]]; then
    echo "Destination already exists: $destination" >&2
    exit 1
fi
mkdir -p "$(dirname "$destination")"
git init --quiet "$destination"
git -C "$destination" remote add origin https://github.com/OpenGene/fastp.git
git -C "$destination" fetch --quiet --depth 1 origin "$reference"
git -C "$destination" checkout --quiet --detach FETCH_HEAD
test "$(git -C "$destination" rev-parse HEAD)" = "$reference"
make -C "$destination" -j"${FASTP_BUILD_JOBS:-2}" \
    CXX="${CXX:-c++}" \
    INCLUDE_DIRS="${FASTP_INCLUDE_DIRS:-}" \
    LIBRARY_DIRS="${FASTP_LIBRARY_DIRS:-}"
"$destination/fastp" --version
