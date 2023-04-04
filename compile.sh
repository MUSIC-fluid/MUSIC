#!/usr/bin/env bash

mkdir -p build
cd build
rm -fr *
#CXX=g++-12 cmake .. -DCMAKE_BUILD_TYPE=Release
CXX=g++-12 cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j4
make install
cd ..
rm -fr build
