#!/usr/bin/env bash
set -euo pipefail

bash EOS/download_Neos2D.sh 2

./MUSIChydro tests/test_1+1D_with_Akihiko/music_input_test_baryon_diffusion

rm -fr EOS/neos_2

if python3 tests/test_1+1D_with_Akihiko/TestOutputFiles.py; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
