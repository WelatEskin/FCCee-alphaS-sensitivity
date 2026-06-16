#!/bin/bash

set -e

mkdir -p data logs figures

echo "Running validation sample..."

./generate_fccee_alphas_root 0.118 10000 12345 data/validation_alphaS_0118.root \
  > logs/validation_alphaS_0118.log 2>&1

echo "Validation sample completed."
echo "Output: data/validation_alphaS_0118.root"