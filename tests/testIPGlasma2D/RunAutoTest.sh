#!/usr/bin/env bash

cd EOS; bash download_hotQCD.sh; cd ..

(cd example_inputfiles/2D_IPGlasma/IPGlasma_2D_testEvent/input; bash download_testIPGevent.sh)

./MUSIChydro tests/testIPGlasma2D/music_input_mode_2

rm -fr EOS/hotQCD

python3 tests/testIPGlasma2D/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
