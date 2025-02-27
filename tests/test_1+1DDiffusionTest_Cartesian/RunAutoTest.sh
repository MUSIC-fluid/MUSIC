#!/usr/bin/env bash

./MUSIChydro tests/test_1+1DDiffusionTest_Cartesian/music_input_1+1DDiffusionTest_x

python3 tests/test_1+1DDiffusionTest_Cartesian/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi

./MUSIChydro tests/test_1+1DDiffusionTest_Cartesian/music_input_1+1DDiffusionTest_y

python3 tests/test_1+1DDiffusionTest_Cartesian/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi

./MUSIChydro tests/test_1+1DDiffusionTest_Cartesian/music_input_1+1DDiffusionTest_z

python3 tests/test_1+1DDiffusionTest_Cartesian/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
