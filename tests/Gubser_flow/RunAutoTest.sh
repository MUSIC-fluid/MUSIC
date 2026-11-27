#!/usr/bin/env bash
set -euo pipefail

./MUSIChydro tests/Gubser_flow/music_input_Gubser

if python3 tests/Gubser_flow/TestOutputFiles.py; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
