#!/usr/bin/env bash
# run from the repository root
set -euo pipefail

# download the hotQCD table only if it is missing, and remove afterwards
# only what this test downloaded
EOS_TABLE=EOS/hotQCD/hrg_hotqcd_eos_binary.dat
if [ ! -f "${EOS_TABLE}" ]; then
    trap 'rm -f "${EOS_TABLE}"; rmdir EOS/hotQCD 2>/dev/null || true' EXIT
    bash EOS/download_hotQCD.sh
fi

bash example_inputfiles/IPGlasma_2D/input/download_testIPGevent.sh

./MUSIChydro tests/testIPGlasma2D/music_input_mode_2

if python3 tests/testIPGlasma2D/TestOutputFiles.py; then
    echo "All tests passed! :)"
    rm -fr ./*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
