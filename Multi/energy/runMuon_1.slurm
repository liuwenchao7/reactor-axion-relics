#!/bin/bash
#SBATCH -J runMuon           # 设置作业名称
#SBATCH -o /public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/muon_30M_6m_out/runMuon_1.out       # 标准输出文件
#SBATCH -e /public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/solid_v2_0_4_2_PE_shielding_bot_side_Lead_5cm_thickness_260622/side_distance_side_15cm_top_bottom_15cm/result/muon_30M_6m_out/runMuon_1.err       # 错误输出文件
#SBATCH -n 1                 # 总任务数（由 runMuon.sh 自行控制并行任务）
#SBATCH -c 15                # 每个任务分配的 CPU 核心数（根据服务器总核数调整）
#SBATCH --mem=85000          # 每个任务的内存分配（单位 MB）
#SBATCH -p compute           # 指定分区
#SBATCH -w node02            # 指定在 node02 上运行

# 激活 conda 环境
source /public/apps/anaconda3/etc/profile.d/conda.sh
conda activate relics3

# 运行脚本
# BASE_RUN=0 ./runMuon.sh -j 25 -f 5 -e 2
./runMuon.sh -j 15 -f 3 -e 100
# ./runMuon.sh -j 50 -f 1 -e 6