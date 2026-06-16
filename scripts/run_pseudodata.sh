#!/bin/bash

set -e

# Pseudo-data production script for the FCC-ee generator-level study.
# This script generates an independent pseudo-data sample at alphaS = 0.118.
# The final thesis production used 1,000,000 events.

NEVENTS=1000000

mkdir -p data logs

echo "Running pseudo-data production"
echo "alphaS  = 0.118"
echo "NEVENTS = ${NEVENTS}"
echo "seed    = 999"
echo "Output  = data/pseudo_data_alphaS_0118.root"

./generate_fccee_alphas_root 0.118 ${NEVENTS} 999 data/pseudo_data_alphaS_0118.root > logs/pseudo_data_alphaS_0118.log 2>&1

echo "Pseudo-data production finished."
echo "Log written to logs/pseudo_data_alphaS_0118.log"
