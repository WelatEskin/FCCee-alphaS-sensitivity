#!/bin/bash

set -e

echo "Building FCC-ee alphaS event generator..."

g++ src/generate_fccee_alphas_root.cc -o generate_fccee_alphas_root \
  $(pythia8-config --cxxflags --libs) \
  $(fastjet-config --cxxflags --libs) \
  $(root-config --cflags --libs)

echo "Build completed successfully."
