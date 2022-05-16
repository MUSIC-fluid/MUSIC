#!/usr/bin/env bash

mkdir -p build
cd build
rm -fr *
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
make install
cd ..
rm -fr build
