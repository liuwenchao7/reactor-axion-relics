#!/bin/bash

mode='S1'
files=100
events=100000
parallel=50
clean=''
justprint=''

while getopts m:f:e:j:cn flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
        f) files=${OPTARG};;
        e) events=${OPTARG};;
        j) parallel=${OPTARG};;
        c) clean='clean';;
        n) justprint='-n';;
    esac
done

export FILEN=$files
export EVENTN=$events
export PMODE=$mode

if [[ -z "${RELICSSIM}" ]]; then
    echo "Error: Must provide directory to RelicsSim"
    exit 1
fi

if [[ -z "${ANALYSIS}" ]]; then
    echo "Error: Must provide directory to RelicsAnalysis"
    exit 1
fi

if [[ -z "${FOLDER}" ]]; then
    echo "Error: Must provide directory to simulation results"
    exit 1
fi

make -j $parallel $clean $justprint
