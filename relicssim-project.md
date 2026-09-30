# RELICS project memory (2026-08-18 更新)

## 项目概述
RELICS (REactor neutrino LIquid xenon Coherent elastic Scattering experiment) 实验的 Geant4 模拟。
位置: `/home/liuwch/R10_v2_clean`
框架: BambooMC (https://github.com/pandax-experiments/BambooMC), Geant4 11.1.1

最终目标：基于已经计算出的流强和反应截面，模拟Relics中MeV量级的反应堆轴子的单步/双步散射信号，得到散射角的分布。加上探测器效应，得到探测的散射角的分布。

实现计划：1.模拟反应堆光子的信号(done)
         2.进行散射角分布等分析(done)
         3.模拟反应堆轴子的信号+analysis(on-going)
         4.add detector response(not started)

**背景**：R10 六面铅层夹心 + 反应堆γ能谱（1.014–7.957 MeV，93 点，`data/gammaON_merged.txt`）。
目标：得到 γ 在 LXe TPC 中**Compton 散射角 θ = ∠(p_in, p_out)**  distribution

## 关键命令

### 生成 XML 配置
```bash
cd /home/liuwch/R10_v2_clean
python3 scripts/relics_xml.py config/cevns_clean.json \
  -o config/cevns_clean.xml --gen gamma --check_overlap
```

### 运行 Source1 / Source2 的 1e6 clean 模拟
```bash
cd /home/liuwch/R10_v2_clean
GAMMAHIST=data/gammaON_merged.txt ./build/BambooMC \
  -c config/cevns_clean.xml -m macro/gammaBeam_source1.mac \
  -n 1000000 -o result/gamma_source1_1e6_clean.root
GAMMAHIST=data/gammaON_merged.txt ./build/BambooMC \
  -c config/cevns_clean.xml -m macro/gammaBeam_source2.mac \
  -n 1000000 -o result/gamma_source2_1e6_clean.root
```

### Sorting (ROOT → HDF5)
```bash
cd /home/liuwch/R10_v2_clean
python3 Multi/energy/Sorting.py \
  --InputFile result/gamma_source1_1e6_clean.root \
  --OutputFile result/gamma_source1_1e6_clean.h5
python3 Multi/energy/Sorting.py \
  --InputFile result/gamma_source2_1e6_clean.root \
  --OutputFile result/gamma_source2_1e6_clean.h5
```

## 模拟配置

### source config（macro）
SimpleGPSGenerator (gamma)
- 类型: planar circle source 
- 半径: 8cm
- 方向: -X
- 能量：`GAMMAHIST=data/gammaON_merged.txt` 

- **source1**： (30,0,-7.55)cm，锥 0–18°（覆盖 TPC r=8.25cm），`gammaBeam_source1.mac`
- **source2**： (200,0,-7.55)cm，锥 0–0.1°（近平行束），`gammaBeam_source2.mac`


### 当前 RELICS clean 探测器几何
- 生效配置：`config/cevns_clean.json` → `config/cevns_clean.xml`；当前模拟不使用完整几何 `config/cevns.xml`。
- World：半长 1000 cm；`World.cc` 中材料为 `G4_AIR`。
- `XenonDetector`：LXe 母体，外半径 102.579 mm、高 178.85 mm、中心 z = -89.425 mm；内部 `XenonVeto1` 为半径 82.5–102.579 mm 的 LXe veto 壳，因此中心区 r < 82.5 mm 作为 TPC 分析区域。
- **保留组件**：XenonDetector、XenonVeto1、BotPMTTeflon、XeSideTeflon1、Shapingrings、TopKee2、Cathode、BottomScreening、Gate 及三组对应 rings。
- **clean 配置未包含**：OuterContainer、Vacuum、Flange、InnerContainer、OuterGasXe、旧 XenonVeto、PMT 阵列、Anode、水屏蔽、混凝土和铅层。
- 内部材料仍可能在 γ 进入 TPC 前散射：source1 的记录占比诊断为 Shapingrings 14.2%、XenonVeto1 19.1%、TopKee2 2.4%、XeSideTeflon1 3.3%、World 3.8%；所以“首次散射必在 TPC”不成立。

## 当前数据产品与 HDF5 结构

### 1e6 clean 样本（`result/`）
- Source1：`gamma_source1_1e6_clean.root`（约 1.9 GB）、`gamma_source1_1e6_clean.h5`（约 1.3 GB）。
- Source2：`gamma_source2_1e6_clean.root`（约 2.6 GB）、`gamma_source2_1e6_clean.h5`（约 1.8 GB）。
- ROOT `mcTree` 是 gammaComp 体积、PreStep/PostStep 动量和 TrackID 的权威来源；这些字段没有完整写入 HDF5。当前 $\theta_{reco}$ 的顶点与排序改用 HDF5 reconstructed clusters，ROOT gammaComp 仅提供匹配事件的 $\theta_{true}$。

### `events` (`event_dtype`)
- `runId, eventId, nHits, nClusters, nPrimaries`。
- TPC/veto 的 `sum/num/max/alt` ER/NR 汇总；`max/alt` 还含位置、时间、parent/type/process 等信息。

### `primaries` (`primary_dtype`)
- `runId, eventId, nPrimaries, primaryEnergy, primaryPx/Py/Pz, primaryX/Y/Z, primaryId, primaryType`。

### `clusters` (`cluster_dtype`)
- `runId, eventId, energy, nr, er`。
- `er_*/nr_*`：`x, y, z, t, parentId, type, parentType, creatorProcess, depositionProcess, volume`。

### `gammaComp` (`gamma_comp_dtype`)
- 一条 HDF5 row 对应一条 Sorting 导出的 γ Compton 记录：`runId, eventId, nGammaComp, x/y/z, px/py/pz, t`。
- HDF5 `gammaComp` **不含** `gammaCompVol`、`gammaCompPrePx/Py/Pz`、`gammaCompTrackId`；按体积过滤和计算真实偏转角必须读取 ROOT。
- `energy_dtype` 只用于 Sorting 中间处理，不单独写入 HDF5。
- dtype 定义：`Multi/energy/dtypes.py`；ROOT→HDF5：`Multi/energy/Sorting.py`。

## 当前分析 notebook

### `analysis/gamma_beam_analysis_R10.ipynb`
- 当前为 23 cells，同时分析 Source1/Source2 的 1e6 clean ROOT 和 HDF5：
  - HDF5：`../result/gamma_source{1,2}_1e6_clean.h5`，用于 events/primaries/clusters、初级能谱、沉积能谱、簇位置和 θ–E 配对。
  - ROOT：`../result/gamma_source{1,2}_1e6_clean.root`，用于 `gammaCompVol`、PreStep/PostStep 动量和逐次 TPC Compton 真值角。
- `compute_inTPC_theta()`：按 `gammaCompVol == "XenonDetector"` 过滤，构造 θ_all/θ_first 并区分单重、多重。
- `pair_theta_energy()`：用 `events.eventId` 完成 ROOT 树行号→G4 eventId 映射，并把 cluster 时间从秒换算为 ns 后配对 E_dep。
- `compute_theta_reco_true()` / `analysis/python_code/plot_theta_reco_vs_true.py`：不使用 TrackID；从 HDF5 选择 TPC ER reconstructed clusters，按 `er_t` 排序并取最早两个 cluster 的位置构造 $\theta_{reco}$。ROOT 不再提供重建顶点，只按 `eventId` 提供首个 TPC Compton 的 $\theta_{true}$。
- 图件统一保存在 `analysis/plot/`；θ–E_dep 和 θ–E_init 系列分别位于其同名子目录。

## 关键文件路径

### 当前 clean 模拟直接使用
- 几何参数及配置：`config/geo_params.json`、`config/cevns_clean.json`、`config/cevns_clean.xml`。
- XML 生成：`scripts/relics_xml.py`。
- Source macros：`macro/gammaBeam_source1.mac`、`macro/gammaBeam_source2.mac`；旧束流配置保存在 `macro/legacy/` 和 `macro/iso_angle/`，不用于当前结果。
- 输入能谱：`data/gammaON_merged.txt`；原始/对照谱为 `data/gammaON.txt`、`data/gammaOFF.txt`。
- 模拟结果：`result/gamma_source{1,2}_1e6_clean.{root,h5}`；对应运行日志位于 `result/`。
- 主分析：`analysis/gamma_beam_analysis_R10.ipynb`。
- 角分布辅助脚本：`analysis/python_code/plot_theta_first_cdf.py`、`analysis/python_code/plot_theta_reco_vs_true.py`。
- ROOT→HDF5：`Multi/energy/Sorting.py`；dtype：`Multi/energy/dtypes.py`。

### 几何输入与位置文件（完整保留）
- 容器 RZ 轮廓：`data/OuterContainerRZ.txt`、`data/InnerContainerRZ.txt`、`data/VacuumRZ.txt`、`data/OuterGasXeRZ.txt`、`data/XenonVetoRZ.txt`。
- LXe 容器轮廓：`data/Lxe_InnerContainer.txt`、`data/Lxe_OutContainer.txt`。
- PMT 位置：`data/topPMTs.txt`、`data/botPMTs.txt`。
- PMT Teflon 开孔/氙区位置：`data/XeInTopPMTTeflon.txt`、`data/XeInBotPMTTeflon.txt`；当前 clean 几何中的 `BotPMTTeflon` 直接读取后者。
- 支撑与法兰位置：`data/shapingrings.txt`、`data/keel_GXe.txt`、`data/keel_LXe.txt`、`data/Flange35.txt`、`data/Flange63.txt`；当前 clean 几何中的 Shapingrings/TopKee2 分别使用 `shapingrings.txt`/`keel_LXe.txt`。

### C++ 实现
- gamma Compton 记录：`src/PandaXSteppingAction.cc`、`include/PandaXDataManager.hh`、`src/PandaXDataManager.cc`。
- 探测器实现：`user/detectors/relics/*.cc/.hh`，包括 World、Cuboid、Cylinder、Prism、Polyhedra、Array、HolesBoard、MeshGrid、PMTsR8520、Monoblock、Material、GetMass。
- 生成器/物理：`user/mc/relics/`。
- 框架核心：`src/SimpleGPSGenerator.cc`、`src/PandaXAnalysis.cc`。

## clean XML 生成方式
```bash
cd /home/liuwch/R10_v2_clean
python3 scripts/relics_xml.py config/cevns_clean.json \
  -o config/cevns_clean.xml --gen gamma --check_overlap
```
- `--gen gamma` 配置 `SimpleGPSGenerator`；束流位置、方向和 Arb 能谱由运行时 macro/GPS 命令设置。
- `--check_overlap` 打开几何重叠检查。
- `XenonDetector` 和 `XenonVeto1` 在 clean JSON 中显式带 `sendname`；不能再沿用旧版“仅 XenonDetector 有 sensitive detector”的描述。

## 单位体系（当前数据）
| 物理量 | ROOT/HDF5 中的单位 | 备注 |
|--------|-------------------|------|
| 能量、primary momentum、gammaComp momentum | keV、keV/c | notebook 中 `M_E_C2 = 511.0` keV |
| 位置 | mm | ROOT gammaComp、HDF5 cluster/gammaComp 均按 mm |
| ROOT `gammaCompT`、HDF5 `gammaComp.t` | ns | Sorting 对 gammaComp 时间直接抄写 |
| HDF5 cluster `er_t/nr_t` | s | 与 gammaComp 配对时需乘 `1e9` 转成 ns |
- `gammaON_merged.txt` 的能量横轴为 MeV，notebook 读取后乘 1000 转为 keV。

### LaTeX 渲染
- `analysis/relics_new.mplstyle` 开启了 `text.usetex: True` + siunitx
- **禁止**在 matplotlib 标签中使用 Unicode 特殊字符（如 `≥`），会导致 LaTeX 编译失败
- 替代方案：使用 ASCII (`>=`) 或 LaTeX 命令 (`$\geq$`)

### Notebook 常见陷阱
- `awkward` 数组与 numpy 操作需谨慎，`.to_numpy()` 前确保数据结构一致
- HDF5 clusters 表中无 `x/y/z` 直接字段，位置用 `er_x/er_y/er_z`
- 空数组的 `.min()/.max()` 会抛 ValueError

### experience

1. **GPS 束流方向：极轴 = `rot2 × rot1`（not rot1×rot2！）**
   - `rot1=(0,0,1) rot2=(0,-1,0)` → 极轴 = (0,-1,0)×(0,0,1) = **(-1,0,0) = 朝向探测器**
   - **1e6 → 714,736 事件（71.5% 命中率）**
   - 验证方法：检查 `primaryPx` 均值必须 < 0

2. **World 体积 = `G4_AIR`（`user/detectors/relics/World.cc` 硬编码）！**
   - 注释掉 Air 体积 ≠ 移除空气！束流仍穿过 10m 空气盒（World half=1000cm）
   - 诊断显示 source1 有 57,682 个 gammaCompton 记录在 World 中（3.8%）
   - 要真正"clean"：World 材料改为真空（`G4_Galactic`）

3. **γ 的 Compton 步能量沉积 = 0** → 灵敏探测器永远看不到它 → 必须用 stepping action 在 γ 自己的散射步上捕获 PostStep 动量

4. **`-n` 参数**：Geant4 11 MT 只打印一次 "Processing event 0..."，实际处理全部事件，勿误判

5. **Arb 能谱文件带 header 行会崩 GPS**（`G4SPSEneDistribution::ArbInterpolate`）→ 必须无 header

6. **`/gps/ang/type beam` 会忽略 Arb 能谱**（固定 1000 keV）→ 用 `iso` + 0.0017 rad (0.1°) 锥

7. **120 GeV 离群点**：每个 1e6 文件恰好 1 个（GPS Arb 插值边缘伪影，0.0001%，可忽略）

### C++ 修改（gamma 出射动量记录）
- `src/PandaXSteppingAction.cc`：`def->GetParticleName()=="gamma"` 且 PostStep 过程=="compt" → `addGammaCompRecord(pos, mom(out), t, trackId, pmom(in), volName)`
- 目前是 **`if(true)` 全体积记录**（诊断模式），未按计划限制在 XenonDetector/XenonVeto1
- `include/PandaXDataManager.hh` + `src/PandaXDataManager.cc`：`nGammaComp, gammaCompX/Y/Z, gammaCompPx/Py/Pz, gammaCompT, gammaCompTrackId, gammaCompPrePx/Py/Pz, gammaCompVol` 分支；gammaComp 容器在 **fillEvent 末尾**清空（不能在 resetData，因 stepping 在事件中填充）
- `Multi/energy/Sorting.py`：写 HDF5 数据集 `gammaComp`（dtype 见 `Multi/energy/dtypes.py`，字段 runId/eventId/nGammaComp/x/y/z/px/py/pz/t）

### 运行命令
```bash
cd /home/liuwch/R10_v2_clean
GAMMAHIST=data/gammaON_merged.txt ./build/BambooMC -c config/cevns_clean.xml -m macro/gammaBeam_source1.mac -n 1000000 -o result/gamma_source1_1e6_clean.root
# Sorting (ROOT→HDF5)
cd Multi/energy && python3 Sorting.py --InputFile ../../result/gamma_source1_1e6_clean.root --OutputFile ../../result/gamma_source1_1e6_clean.h5
```

### 当前结果（1e6 事件）
| | source1 | source2 |
|---|---|---|
| 事件数 | 714,736 (71.5%) | 930,959 (93%) |
| ROOT 大小 | 1.9GB | 2.5GB |
| primaryPx 均值 | -1403 | -1439 |
| gammaComp 记录/事件 | ~2.15 | — |
| 散射在 TPC (r<82.5) | 80.9% | — |

### 分析 notebook（`analysis/gamma_beam_analysis_R10.ipynb`）
**偏转角**：γ 在 LXe TPC（volume==XenonDetector）内每一次 Compton 相互作用的偏转角
> θ = ∠(p_in, p_out)，p_in=该步入射动量（PreStep），p_out=出射动量（PostStep）
> 无论 γ 之前在别处散射过 n 次（第 n+1 次在 TPC 内也算）；**单重 + 多重都要**

- `compute_inTPC_theta(root_path)`：从 **ROOT 直接读**（h5 的 gammaComp 无体积、无 PrePx/Py/Pz，**不够用**！）→ 过滤 volume=="XenonDetector" → 每条记录算 θ → 每事件取时间最早（θ_first）+ 全部（θ_all）→ 单重/多重切分

### 结果（1e6）
| | Source 1 | Source 2 |
|---|---|---|
| 有≥1次TPC内Compton | 442,778 (61.9%) | 649,632 (69.8%) |
| 单重/多重 | 197,359/245,419 | 277,382/372,250 |

- **单重(≈68°) vs 多重(≈48°) 显著不同**（物理：小角前向散射 γ 保留能量继续穿越→多重；大角→转向/逃逸→单重）
- **Source1 与 Source2 分布几乎重合**（θ 相对各 γ 自身入射方向定义，与源几何无关）→ 自洽性 ✓

**🚨 eventId vs 树行号 bug 已修复（重要）**
- h5 的 `eventId` 是 G4 事件号（0..999998，如 `[0,1,2,4,8...]`），**≠ 树行号**（0..N-1）！
- 修复：`row_to_evtid = data['events']['eventId']`（events 表是树序）把行号映射到 G4 eventId
- h5 `er_volume` 是 bytes 对象数组，比较必须 `== b'XenonDetector'`（str 全 False）

**add θ_reco（当前 reconstructed-cluster 版本）**
- 对每个事件选择 `energy > 0` 且 `er_volume == b'XenonDetector'` 的 HDF5 ER clusters；要求至少两个有效 cluster，按 `er_t` 排序，束流轴投影仅作时间完全相同时的稳定 tie-break。
- 以最早两个 cluster 的重建位置为 $r_1,r_2$，使用已知束流轴 $\hat d_{beam}=(-1,0,0)$，定义 $\theta_{reco}=\angle(\hat d_{beam},r_2-r_1)$。不再使用“名义源中心 → $r_1$”作为入射方向。
- ROOT gammaComp 只用于按 `eventId` 匹配 $\theta_{true}=\angle(p_{in,1},p_{out,1})$；重建选择、顶点和入射方向均不读取 TrackID 或 gammaComp X/Y/Z。
- 单 cluster 事件无法用两点连线重建散射角；当前 `er_t` 排序是 cluster-level baseline，尚不等于真实探测器中已可靠解决的交互顺序。

**reconstructed-cluster 样本术语（2026-08-25 明确）**
- **multi-scatter sample**：TPC 内有效 ER cluster 数量 $N_{cluster}\ge2$；$\theta_{reco}$ 使用最早两个 clusters。
- **double-scatter sample**：$N_{cluster}=2$；它是 multi-scatter sample 的严格子集。
- ROOT 中恰好 1 条或至少 2 条 TPC gammaComp 记录的分组只称为 **truth-record multiplicity diagnostic**，不再用它定义 reconstructed single/multi/double sample。
- reconstructed $\theta_{reco}$ PDF/CDF 直接使用全部合格 HDF5 cluster 事件，不要求存在 ROOT $\theta_{true}$ 匹配；只有 $\theta_{reco}-\theta_{true}$ 分辨率检查才要求 ROOT truth match。

**pre-cluster double-first 对照定义**
- `theta_doublefirst_inTPC_pdf_cdf.png` 恢复改用 cluster 之前的 truth 定义：ROOT 中恰好 2 条 `gammaCompVol == XenonDetector` 的 Compton 真值记录，按 `gammaCompT` 取第一条，画 $\theta_{first}=\angle(p_{in},p_{out})$。
- Source1/2 的事件数分别为 135,256/202,184，median $\theta_{first}=53.982°/53.697°$。
- 当前 $N_{cluster}=2$ 的 reconstructed double 图另存为 `theta_double_cluster_reco_inTPC_pdf_cdf.png`，避免与 pre-cluster truth 图混淆。

**新旧 multi-first PDF 形状差异**
- 旧图选择 ROOT 中 ≥2 条 TPC Compton 真值记录，并画首次 Compton 的真实动量偏转角；“发生后续 Compton”会偏好较小的首次散射角，因此峰后呈明显的先快后慢下降。
- 新图选择 $N_{cluster}\ge2$，画固定束流轴与两个 cluster 质心连线的几何角。Source2 的 486,914 个 cluster-multi 事件中，只有 74.0% 真正有 ≥2 条 TPC Compton 真值记录；25.1% 只有 1 条，1.0% 一条也没有。后两类的额外 clusters 可来自电子轨迹分裂、光电吸收、次级粒子或其他相邻沉积，不再服从“首次 Compton 后继续散射”的条件分布。
- 新 multi 还混合了 224,591 个 $N_{cluster}=2$ 事件（46.1%，median $\theta_{reco}=73.502°$）和 ≥3-cluster 事件；高角度、宽分布的 double-cluster 成分抬高了中大角尾部。再叠加 cluster 能量加权质心和 `er_t` 排序误差，峰后的曲率被填平，看起来更接近线性下降。

**Compton-like cluster multi 对照**
- 额外图 `theta_multifirst_compton_cluster_inTPC_pdf_cdf.png` 先过滤 `er_creatorProcess == b'compt'`，要求事件内至少两个这样的 TPC ER clusters，再按 `er_t` 取最早两个构造 $\theta_{reco}$。
- Source1：N=212,884，median $\theta_{reco}=47.785°$；Source2：N=324,746，median $47.280°$。
- 过滤后 PDF 从全部-cluster multi 的近线性下降恢复为更接近旧 truth multi-first 的“峰后先快后慢”，说明 photoelectric/其他来源 clusters 确实显著抬高了中大角尾部。
- 这里的 `er_creatorProcess` 只是 cluster 内最大能量 ER hit 的 creator-process 标签，所以应称 Compton-like，而不能解释为纯 Compton cluster truth。

**plot**
- deflection_vs_energy_vs_einit.png uses 2 panels（count density / sixty thousand points E_init colored）
- 三张 θ-E 图统一范围 0–2500 keV
- energy_deposition use log scale
- primary_energy_spectrum compared by bins, verify GPS sampling is correct

**single scatter samples have Compton+photoelectric absorption**
- "单重"（1 次 TPC 内 Compton）事件里 **~60% 有 ≥2 个 cluster**（s1: 2簇44.8%+≥3簇15.5%；s2: 46.6%+14.9%）
- 若想要"纯单次逃逸"拓扑，应过滤到 **恰好 1 个 cluster** 的事件

**🔍 measurement selection bias**
- 当前 measurable angle 样本要求至少两个 HDF5 TPC ER clusters，并且存在可匹配的首个 TPC Compton 真值角；它不再等同于“ROOT 中至少两次 TPC Compton”的 multi-scatter 样本。

### θ–E plot

**composition of θ_all（ROOT truth-record multiplicity diagnostic）**：θ_all 包含每一次 TPC gammaComp 真值相互作用
- 源1：856,234 条；其中来自 ≥2 条 TPC gammaComp 记录事件的占 77.0%（658,875），1 条记录事件占 23.0%（197,359）。
- 源2：1,282,956 条；其中来自 ≥2 条记录事件的占 78.4%（1,005,574），1 条记录事件占 21.6%。

| 样本 | E_dep vertical axis | E_init vertical axis |
|------|-----------|------------|
| θ_all | `deflection_vs_energy_vs_einit_{s}.png` | `deflection_vs_einit_{s}.png` |
| θ_multi | `..._multiscatter_{s}.png` | `deflection_vs_einit_multiscatter_{s}.png` |
| θ_multi_first | `..._multifirst_{s}.png` | `deflection_vs_einit_multifirst_{s}.png` |

- **历史变量名的过滤逻辑**：`multi = res['cnt'][res['ev_all']] >= 2` 实际表示 ROOT 中 ≥2 条 TPC gammaComp 真值记录；`first` 是按 `(event, time)` 取第一条。它不是当前 cluster-defined multi-scatter。
- 该 truth-record `multi+first` 记录数：源1 = 245,419、源2 = 372,250。
- E_dep vertical axis bin：θ 2°/格 × E_dep 25 keV/格；E_init vertical axis bin：θ 2°/格 × E_init 100 keV/格

**E_dep 定义**
- **E_dep 轴 = local energy deposit**（DBSCAN 簇总能量 `clusters['energy']`），not γ energy loss（γ 损失 = |p_in|−|p_out|，可由 `gammaCompPrePx/..` 与 `gammaCompPx/..` 算出，单位 keV/c）。
- **`pair_theta_energy` 配对单位**：`clusters['er_t']` use second（Sorting.py 除以 `scale_to_ns` 后）, `gammaCompT` use ns（~6.4）
- **unify**：`c_t = cl['er_t'][cm] * 1e9` unify to ns。簇≈γ损失(±5 keV) **74–84%**；`multi_first` E_dep 中位数 237 → **533 keV**（≈γ损失 512）。剩余 ~20% 簇>γ损失 = **合并簇**（散射后 γ 在 ~3mm/10ns 内被再次散射/光电吸收并入同簇，物理合理）。

**E_init 图水平亮线解释**
- θ–E_init 图里的水平亮线（~1.4/1.8/2.7/6.1 MeV）= 输入反应堆 γ 谱的**强谱线**（`gammaON_merged.txt` 峰：1.33/1.47、1.77、2.75、6.13 MeV）投影。θ 的 0–180° 全覆盖是**所有能量**共有的 Klein–Nishina 散射行为；只是峰位统计量大 → 形成横跨全 θ 的亮线，连续区计数低仅为微弱带。

### 无 TrackID 的 reconstructed-cluster $\theta_{reco}$–$\theta_{true}$ 检查

**目标与事件选择**
- 实现位于 `analysis/python_code/plot_theta_reco_vs_true.py`，notebook 为 `analysis/gamma_beam_analysis_R10.ipynb`；重建不读取、不使用 `gammaCompTrackId`。
- 从 HDF5 `clusters` 中选择 `energy > 0`、`er_volume == b'XenonDetector'` 且 `er_x/y/z/t` 有限的 ER clusters。按 `(eventId, er_t, beam-axis projection)` 排序，仅保留至少两个 cluster 的事件；最早两个 cluster 定义 $r_1,r_2$。
- 已知入射方向为固定束流轴 $\hat d_{beam}=(-1,0,0)$，定义 $\theta_{reco}=\angle(\hat d_{beam},r_2-r_1)$。ROOT 只过滤 `gammaCompVol == "XenonDetector"`、按 `gammaCompT` 找到该事件首次相互作用，并提供 $\theta_{true}=\angle(p_{in,1},p_{out,1})$；gammaComp X/Y/Z 不参与重建。
- HDF5 的 `events.eventId` 用于将 ROOT 树行号映射到 G4 eventId。每个匹配事件只构造一个 $\Delta\theta=\theta_{reco}-\theta_{true}$。

**全样本结果**

| | Source 1 | Source 2 |
|---|---:|---:|
| 入选匹配事件数（TPC 内 ≥2 个 ER cluster） | 320,114 | 481,952 |
| mean signed $\Delta\theta$ | 2.217° | 2.213° |
| median signed $\Delta\theta$ | 1.017° | 0.323° |
| mean absolute $\Delta\theta$ | 14.614° | 10.034° |
| median absolute $\Delta\theta$ | 8.423° | 1.038° |
| 68% containment of $|\Delta\theta|$ | 12.734° | 3.974° |
| 90% containment of $|\Delta\theta|$ | — | 34.378° |
| corr($\theta_{true},\theta_{reco}$) | 0.8197 | 0.8445 |

Source2 的束流仅有 ±0.1° 发散，固定 $-x$ 轴适合作为已知入射方向；改动后其 signed median 从 0.836° 降为 0.323°，median absolute 从 1.781° 降为 1.038°。Source1 有 ±18° 发散，固定轴只代表名义轴而不是逐事件方向，因此其分辨率仍明显较差。两者都有长尾，所以 mean absolute/std 不能单独代表中心分辨率，优先同时报告 signed median、median absolute 与 68% containment。

**Source2 中心峰 Gaussian 检查（0.25° bin）**
- Gaussian $\mu_{fit}=0.306°$、$\sigma_{fit}=0.512°$，Gaussian FWHM = $1.206°$，但 $\chi^2/\mathrm{ndf}=2320.6$，说明“很窄且非对称的中心峰 + 长尾”远非单 Gaussian。
- reconstructed-cluster 版本中 $|\Delta\theta|\le10^{-8}°$ 的比例为 0；旧 MC 顶点版本的数值 delta 尖峰已经消失。
- 推荐主结果使用 median $\Delta\theta=0.323°$、median $|\Delta\theta|=1.038°$ 和 68% containment $=3.974°$；Gaussian 宽度只作为中心核的辅助描述。

**cluster 排序方法诊断**
- 当前 baseline 按 `er_t` 排序；只有 `er_t` 完全相同时，才按束流轴投影选择更上游的 cluster 作为稳定 tie-break。旧“名义源距离排序”比较属于已废弃的名义源方向版本，不再用于当前重建。
- 这不证明真实数据中可以无歧义恢复先后顺序；HDF5 `er_t` 是模拟重建 cluster 的时间字段，后续仍需研究实际时间分辨率、cluster 合并和顺序判别效率。

**旧 MC 真值顶点方法的 TrackID 事后诊断（历史记录，不参与当前 cluster 重建）**
- Source2 最早两条 TPC Compton 记录属于同一 TrackID：**352,829 / 372,250 = 94.7828%**；其中 **337,958 / 352,829 = 95.7852%** 满足 $|\Delta\theta|\le10^{-8}°$。
- 属于不同 TrackID：**19,421 / 372,250 = 5.2172%**；其中满足 $|\Delta\theta|\le10^{-8}°$ 的事件为 **0**。
- 因而旧方法的全部 337,958 个近零事件都来自同 TrackID；这解释了旧真值顶点图的零点尖峰。该标签分解不用于当前 reconstructed-cluster 样本的筛选或性能定义。

**当前图件**
- 基础分布（`analysis/plot/`）：`primary_energy_spectrum.png`、`energy_deposition_spectrum.png`、`hit_position_distribution.png`。
- ROOT truth-angle diagnostic：`theta_first_inTPC.png`、`theta_first_inTPC_cdf.png`。
- reconstructed-cluster $\theta_{reco}$ PDF/CDF：`theta_multifirst_inTPC_pdf_cdf.png`（$N_{cluster}\ge2$）和 `theta_double_cluster_reco_inTPC_pdf_cdf.png`（$N_{cluster}=2$）。multi 文件名中的 `first` 为兼容旧命名，图中实际量是 $\theta_{reco}$。
- Compton-like cluster multi：`theta_multifirst_compton_cluster_inTPC_pdf_cdf.png`，选择 `er_creatorProcess == b'compt'` 且 $N_{cluster}^{compt}\ge2$。
- pre-cluster ROOT truth 对照：`theta_doublefirst_inTPC_pdf_cdf.png`，表示恰好 2 条 TPC Compton 真值记录的首次真值动量角。三张图均由 `analysis/python_code/plot_theta_first_cdf.py` 生成。
- cluster PDF/CDF 统计：Source1 multi/double = 323,639/153,215，median $\theta_{reco}=59.611°/73.578°$；Source2 = 486,914/224,591，median $59.151°/73.502°$。
- reconstructed-cluster θ_reco 对比：`theta_reco_vs_true.png`、`theta_reco_minus_true_hist.png`（线性百分数）、`theta_reco_minus_true_hist_log.png`（对数百分数）、`theta_reco_minus_true_gaussian_fwhm.png`（仅 Source2 Gaussian 中心核诊断）。
- 兼容文件名：`theta_geo_vs_out.png` 与当前 `theta_reco_vs_true.png` 同步；`θ_reco_vs_true.png` 为旧文件，不作为当前主图名。
- θ–E_dep（`analysis/plot/θ_vs_E_dep/`）：`deflection_vs_energy.png`；all 为 `deflection_vs_energy_vs_einit_{s1,s2}.png`，另有 `_multiscatter_{s1,s2}`、`_multifirst_{s1,s2}`，共 6 张 Source 分图。
- θ–E_init（`analysis/plot/θ_vs_E_init/`）：`deflection_vs_einit_{s1,s2}.png`，以及 `_multiscatter_{s1,s2}`、`_multifirst_{s1,s2}`，共 6 张。
- 场景草图（`analysis/sketch/`）：`sketch_scenario_v2.png`、`sketch_scenario_v2_source2.png`。
