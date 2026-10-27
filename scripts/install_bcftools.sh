#!/usr/bin/env bash
# Install bcftools pinned to commit 385820f with native --check-ref C uppercase support
set -euo pipefail

TARGET_PREFIX="${1:-}"

if [ -z "${TARGET_PREFIX}" ]; then
    if [ -n "${CONDA_PREFIX:-}" ]; then
        TARGET_PREFIX="${CONDA_PREFIX}"
    elif [ -n "${VIRTUAL_ENV:-}" ]; then
        TARGET_PREFIX="${VIRTUAL_ENV}"
    elif [ -w "/usr/local" ]; then
        TARGET_PREFIX="/usr/local"
    else
        TARGET_PREFIX="$(dirname "$(dirname "$(which bcftools)")")"
    fi
fi

echo "=== Target installation prefix: ${TARGET_PREFIX} ==="
mkdir -p "${TARGET_PREFIX}/bin" "${TARGET_PREFIX}/libexec/bcftools"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "--> Cloning htslib..."
git clone --depth 1 --recurse-submodules https://github.com/samtools/htslib.git "${TMP_DIR}/htslib"

echo "--> Cloning bcftools (commit 385820f)..."
git clone https://github.com/samtools/bcftools.git "${TMP_DIR}/bcftools"
cd "${TMP_DIR}/bcftools"
git checkout 385820f

echo "--> Building bcftools and plugins with prefix=${TARGET_PREFIX}..."
make prefix="${TARGET_PREFIX}" HTSDIR="${TMP_DIR}/htslib" -j"$(nproc)"
make prefix="${TARGET_PREFIX}" install

echo "=== Verifying installation ==="
"${TARGET_PREFIX}/bin/bcftools" --version
echo "=== Verifying plugins ==="
"${TARGET_PREFIX}/bin/bcftools" plugin -l | grep -E "fill-tags|setGT|fixploidy"
echo "Successfully installed bcftools and plugins with --check-ref C support."
