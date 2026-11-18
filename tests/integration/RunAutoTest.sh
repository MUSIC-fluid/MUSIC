#!/usr/bin/env bash
# High-level integration tests. Run from the repository root after building
# MUSIChydro (make install). The memory_history test also needs the driver
# from a -Dunittest=ON build in build_unittest/; it is skipped otherwise.
# Extra arguments are passed on, e.g. a list of tests or --workdir DIR.

driver=()
if [ -x build_unittest/src/integration_hydro_memory.e ]; then
    driver=(--memory-driver build_unittest/src/integration_hydro_memory.e)
fi

python3 tests/integration/music_integration_tests.py --music ./MUSIChydro \
    "${driver[@]}" "$@"
