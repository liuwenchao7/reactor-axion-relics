#!/bin/bash

if [[ -z "${GEO}" ]]; then
    echo "Using default geometry setting"
else
    echo "Using customized geometry setting $GEO"
fi

export RELICSSIM=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm     # Need update
export ANALYSIS=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/Multi/energy  # Need update
# export PYTHONPATH="$RELICSSIM/scripts:$PYTHONPATH"

parallel=50
files=1000
events=10000
clean=''
justprint=''

while getopts f:e:j:cn flag
do
    case "${flag}" in
        j) parallel=${OPTARG};;
        f) files=${OPTARG};;
        e) events=${OPTARG};;
        c) clean='-c';;
        n) justprint='-n';;
    esac
done

norm='norm.json'
jq -n ".N=$(expr $files \* $events)" > $norm


# # Neutron when reactor ON
export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/gamma_OFF_direct/gammaOFF_Y_direct_30B
reactors="OFF"
topsides="SIDE_Y_OFF"
for reactor in $reactors; do
    for topside in $topsides; do
        export SUFFIX="_gamma_${reactor}_${topside}"
        jq ".$SUFFIX=10000" $norm | sponge $norm
        ./run.sh -m gamma -r $reactor -s $topside -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 1000) $clean $justprint
    done
done

# reactors="OFF"
# topsides="TOP_Z_OFF"
# for reactor in $reactors; do
#     for topside in $topsides; do
#         export SUFFIX="_gamma_${reactor}_${topside}"
#         jq ".$SUFFIX=10000" $norm | sponge $norm
#         ./run.sh -m gamma -r $reactor -s $topside -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 1000) $clean $justprint
#     done
# done