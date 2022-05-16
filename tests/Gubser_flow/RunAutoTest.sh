#!/usr/bin/env bash

./MUSIChydro tests/Gubser_flow/music_input_Gubser

python3 tests/Gubser_flow/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
