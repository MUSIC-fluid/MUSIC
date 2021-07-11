#!/usr/bin/env bash

REFERENCE_DIR="tests/Gubser_flow/SimulationReferences"

./MUSIChydro tests/Gubser_flow/music_input_Gubser

DIFF_1=$(diff ./Gubser_flow_check_tau_1.dat $REFERENCE_DIR/Gubser_flow_check_tau_1.dat)
DIFF_2=$(diff ./Gubser_flow_check_tau_1.2.dat $REFERENCE_DIR/Gubser_flow_check_tau_1.2.dat)
DIFF_3=$(diff ./Gubser_flow_check_tau_1.5.dat $REFERENCE_DIR/Gubser_flow_check_tau_1.5.dat)
DIFF_4=$(diff ./Gubser_flow_check_tau_2.dat $REFERENCE_DIR/Gubser_flow_check_tau_2.dat)

N_PASSED=0
N_TESTS=0

if [ "${DIFF_1}" == "" ]
then
    N_PASSED=$((${N_PASSED}+1))
else
    echo "Test Gubser_flow_check_tau_1.dat failed"
fi
N_TESTS=$((${N_TESTS}+1))

if [ "${DIFF_2}" == "" ]
then
    N_PASSED=$((${N_PASSED}+1))
else
    echo "Test Gubser_flow_check_tau_1.2.dat failed"
fi
N_TESTS=$((${N_TESTS}+1))

if [ "${DIFF_3}" == "" ]
then
    N_PASSED=$((${N_PASSED}+1))
else
    echo "Test Gubser_flow_check_tau_1.5.dat failed"
fi
N_TESTS=$((${N_TESTS}+1))

if [ "${DIFF_4}" == "" ]
then
    N_PASSED=$((${N_PASSED}+1))
else
    echo "Test Gubser_flow_check_tau_2.dat failed"
fi
N_TESTS=$((${N_TESTS}+1))

N_FAILED=$(($N_TESTS-$N_PASSED))
if [[ $N_FAILED -eq 0 ]]
then
    echo "All $N_TESTS tests passed! :)"
    rm -fr ./*.dat
else
    echo ""
    echo "Tests FAILED :("
    echo "$N_FAILED/$N_TESTS tests FAILED"
    exit 1
fi
echo ""
