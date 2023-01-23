#!/usr/bin/env bash

./MUSIChydro tests/testIPGlasma2D/music_input_mode_2

python3 tests/testIPGlasma2D/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
