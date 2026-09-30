#!/bin/bash

if [[ -z "${GEO}" ]]; then
    echo "Using default geometry setting"
else
    echo "Using customized geometry setting $GEO"
fi

export RELICSSIM=$HOME/RelicsSim
export ANALYSIS=$HOME/RelicsAnalysis/SimAnalysis/optical
# export PYTHONPATH="$RELICSSIM/scripts:$PYTHONPATH"

signal='S1'
parallel=50
files=100
events=100
clean=''
justprint=''

while getopts f:e:s:j:cn flag
do
    case "${flag}" in
        j) parallel=${OPTARG};;
        f) files=${OPTARG};;
        e) events=${OPTARG};;
        s) signal=${OPTARG};;
        c) clean='-c';;
        n) justprint='-n';;
    esac
done

# optical simulation
export FOLDER=$HOME/results/optical
export SUFFIX="_$signal"
./run.sh -m optical -f $files -j $parallel -e $events $clean $justprint
