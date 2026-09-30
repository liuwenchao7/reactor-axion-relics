#!/bin/bash

mode='muon'
isotope='geantino'
component='steel'
reactor='ON'
topside='TOP'
files=100
events=100000
parallel=50
clean=''
justprint=''

while getopts m:i:t:r:s:f:e:j:cn flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
        i) isotope=${OPTARG};;
        t) component=${OPTARG};;
        r) reactor=${OPTARG};;
        s) topside=${OPTARG};;
        f) files=${OPTARG};;
        e) events=${OPTARG};;
        j) parallel=${OPTARG};;
        c) clean='-c';;
        n) justprint='-n';;
    esac
done

export RELICSSIM=$HOME/RelicsSim
export ANALYSIS=$HOME/relicsanalysis/muonneutron
export FOLDER=$HOME/results/$mode

# PE thickness scan

increment=10
innPE=($(seq -f '%02g' 0 $increment 60))
outPE=($(seq -f '%02g' 0 $increment 60 | sort -r))

for i in ${!innPE[@]}; do
    export GEO="| jq '.inn_pe=${innPE[i]}' | jq '.out_pe=${outPE[i]}'";
    export SUFFIX="_${innPE[i]}_${outPE[i]}"
    ./run.sh -m $mode -r $reactor -i $isotope -t $component -s $topside -f $files -j $parallel -e $events $clean $justprint
done
