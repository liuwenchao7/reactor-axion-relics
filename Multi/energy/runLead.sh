#!/bin/bash

if [[ -z "${GEO}" ]]; then
    echo "Using default geometry setting"
else
    echo "Using customized geometry setting $GEO"
fi

export RELICSSIM=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm     # Need update
export ANALYSIS=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/Multi/energy   # Need update
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

# # Muon induced neutron and gamma
# export FOLDER=$HOME/RelicsSim_v/RelicsSim_v5/results/muon   # Need update
# export SUFFIX="_muon"
# jq ".$SUFFIX=500" $norm | sponge $norm
# # ./run.sh -m muon -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 50) $clean $justprint

# # CRN induced neutron and gamma
# export FOLDER=$HOME/RelicsSim_v/RelicsSim_v5/results/CRN   # Need update
# export SUFFIX="_CRN"
# jq ".$SUFFIX=500" $norm | sponge $norm
# # ./run.sh -m CRN -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 50) $clean $justprint


# # Material background from metal, PE, and lXe
# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material/PS_50M_K40_1B
# isotopes="U238 Th232"
# components="PS"
# # copper particles number is about 1M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=100000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 10000) $clean $justprint
#     done
# done

# isotopes="K40"
# components="PS"
# # copper particles number is about 1M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=100000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100000) $clean $justprint
#     done
# done

export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material/Lead_50B
isotopes="Pb210"
components="lead"
# copper particles number is about 1M
for isotope in $isotopes; do
    for component in $components; do
        export SUFFIX="_${isotope}_${component}"
        jq ".$SUFFIX=100000" $norm | sponge $norm
        ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 1) -j $parallel -e $(expr $events \* 1000000) $clean $justprint
    done
done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/steel
# isotopes="U238 Th232 Co60 K40 Cs137 U235"
# components="steel"
# # steel particles number is about 10M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=10000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 1000) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/Teflon
# isotopes="U238 Th232 Co60 K40 Cs137 U235"
# components="Teflon"
# # particles number is about 1M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=1000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/Flange
# isotopes="U238 Th232 Co60 K40 Cs137 U235"
# components="Flange"
# # particles number is about 1M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=1000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/PMTbase
# isotopes="U238 Th232 Co60 K40 Cs137 U235"
# components="PMTbase"
# # particles number is about 1M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=1000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/PMTwindow
# isotopes="U238 Th232 Co60 K40 Cs137 U235"
# components="PMTwindow"
# # particles number is about 10M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=10000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 1000) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/PE_100M_K40_1B
# isotopes="K40"
# components="PE"
# # particles number is about 100M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=100000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 10000) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/lead_300M_self_shielding
# isotopes="Pb210"
# components="lead"
# # particles number is about 300M
# for isotope in $isotopes; do
#     for component in $components; do
#         export SUFFIX="_${isotope}_${component}"
#         jq ".$SUFFIX=100000" $norm | sponge $norm
#         ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 10000) $clean $justprint
#         # ./run.sh -m material -i $isotope -t $component -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100) $clean $justprint
#     done
# done

# export FOLDER=/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/Material_10M_steel_window_1M_others/lXe
# # particles number is about 1M
# export SUFFIX="_Kr85_lXe"
# jq ".$SUFFIX=1000" $norm | sponge $norm
# ./run.sh -m material -i Kr85 -t lXe -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100) $clean $justprint

# export SUFFIX="_Rn222_lXe"
# jq ".$SUFFIX=1000" $norm | sponge $norm
# ./run.sh -m material -i Rn222 -t lXe -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 100) $clean $justprint


python3 genNormalize.py --InputFile $norm --paramsFile "$RELICSSIM/config/geo_params.json" --GeoFile "$RELICSSIM/config/cevns.json" --OutputFile normalization.json
#python3 genNormalize.py --InputFile $norm --GeoFile $RELICSSIM/config/cevns.json --OutputFile normalization.json

# if [ -z $justprint ]
# then
#     python3 genTemplate.py --InputFolder $HOME/results --InputNorm normalization.json --OutputFile templateSim.npz
# fi



# ignore this particles

# # Neutron when reactor ON and OFF
# export FOLDER=$HOME/results/neutron
# reactors="ON OFF"
# topsides="TOP SIDE"
# for reactor in $reactors; do
#     for topside in $topsides; do
#         export SUFFIX="_neutron_${reactor}_${topside}"
#         jq ".$SUFFIX=100" $norm | sponge $norm
#         ./run.sh -m neutron -r $reactor -s $topside -f $(expr $files \* 10) -j $parallel -e $(expr $events \* 10) $clean $justprint
#     done
# done

# # Gamma when reactor ON
# export FOLDER=$HOME/results/gamma
# reactors="ON"
# topsides="SIDE"
# for reactor in $reactors; do
#     for topside in $topsides; do
#         export SUFFIX="_gamma_${reactor}_${topside}"
#         jq ".$SUFFIX=1000" $norm | sponge $norm
#         ./run.sh -m gamma -r $reactor -s $topside -f $(expr $files \* 1) -j $parallel -e $(expr $events \* 10) $clean $justprint
#     done
# done