#!/bin/bash
set -e # Exit immediately if any command or unit test fails
cd ..

expected="tests/Snakefile_results"
default="results/tests/Snakefile/default"
stream="results/tests/Snakefile/stream"


rm -r "$default"
snakemake -c 8 --drop-metadata --configfile tests/config.unit_test.yaml \
  --config results_dir="$default"

# Run tests, compare default to expected
DIR_1="$expected" DIR_2="$default" python -m unittest tests.test_Snakefile

rm -r "$stream"
snakemake -c 8 --drop-metadata --configfile tests/config.unit_test.yaml \
  --config  use_stream=true results_dir="$stream"

# Run tests, compare stream to expected (keep steam and default in sync)
DIR_1="$expected" DIR_2="$stream" python -m unittest tests.test_Snakefile
