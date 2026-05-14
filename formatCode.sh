#!/usr/bin/env bash

if command -v clang-format &> /dev/null
then
    echo "formatting code ..."
    find src -iname '*.h' -o -iname '*.cpp' | xargs clang-format -i -style=file
else
    echo "clang-format not found, skip formatting code"
fi
