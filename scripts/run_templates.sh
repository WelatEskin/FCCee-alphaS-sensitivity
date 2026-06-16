#!/bin/bash

set -e

NEVENTS=1000000

mkdir -p data logs figures

echo "Starting FCC-ee alphaS template production"
echo "Events per sample: ${NEVENTS}"

./generate_fccee_alphas_root 0.110 ${NEVENTS} 101 data/template_alphaS_0110.root > logs/template_alphaS_0110.log 2>&1
./generate_fccee_alphas_root 0.114 ${NEVENTS} 102 data/template_alphaS_0114.root > logs/template_alphaS_0114.log 2>&1
./generate_fccee_alphas_root 0.118 ${NEVENTS} 103 data/template_alphaS_0118.root > logs/template_alphaS_0118.log 2>&1
./generate_fccee_alphas_root 0.122 ${NEVENTS} 104 data/template_alphaS_0122.root > logs/template_alphaS_0122.log 2>&1
./generate_fccee_alphas_root 0.126 ${NEVENTS} 105 data/template_alphaS_0126.root > logs/template_alphaS_0126.log 2>&1
./generate_fccee_alphas_root 0.130 ${NEVENTS} 106 data/template_alphaS_0130.root > logs/template_alphaS_0130.log 2>&1

./generate_fccee_alphas_root 0.118 ${NEVENTS} 999 data/pseudo_data_alphaS_0118.root > logs/pseudo_data_alphaS_0118.log 2>&1

echo "All template and pseudo-data samples completed successfully."