#!/usr/bin/env bash


find src -iname '*.h' -o -iname '*.cpp' | xargs clang-format -i -style=file
