# 1 First step,  Determine your file path, get the geometry.
Generates the geometry and the corresponding particle emission locations.

```python
python3 ${path_to_RelicsSim}/scripts/relics_xml_first.py ${path_to_RelicsSim}/config/geo_params.json -o ${path_to_RelicsSim}/config/cevns.xml [ --gen muon ] [ --run_number 1234 ] [ --enable_track ] [ --force_sd ] [ --check_overlap ] [ --save_txt ] [ --optical ] [ --get_geo_mass ]
```

Attention! the option in [] have different function.Please add as needed.
-- gen muon, means  Muon Generator. You can change it to CRN or material.
--run_number 1234, means simulating random numbers, which can be ignored in general.
--enable_track, Whether to enable throughput recording.
--force_sd, Record the volume of energy deposited.
--check_overlap, Check if the geometry overlaps.
--save_txt, Whether to save information such as PMT arrays.
--optical, Used for optical simulation.
--get_geo_mass, Record the mass within the geometry.

So Generally run！!In ${path_to_RelicsSim}, run code.
```python
python3 ${path_to_RelicsSim}/scripts/relics_xml_first.py ${path_to_RelicsSim}/config/geo_params.json -o ${path_to_RelicsSim}/config/cevns.xml  --gen muon  --check_overlap 
```

After this step, the geometry generated according to the {path_to_RelicsSim}/config/geo_params.json parameters is obtained.
Also updated the data in {path_to_RelicsSim}/data and the position of particle emission in {path_to_RelicsSim}/macro.

# 2 Second step, Compile executable file
Under the path of RelicsSim
```bash
mkdir build
cd build
cmake -DENABLE_DETECTOR_SETS=relics -DENABLE_USER_MC=pandax,relics -DCMAKE_BUILD_TYPE=Debug ..
make
```

# 3 Third step, Simulation
Make sure which backgroun you want to simulate.
## Muon
For Muon, you need make sure you generator in first step is right. 
```bash
${path_to_RelicsSim}/build/BambooMC -c ${path_to_RelicsSim}/config/cevns.xml [ -n 100 ] [ -o muon.root ] [ -i ]
```
-n, means number;
-o, means output file;
-i, means visualization; Ignore! Unless you want to look at the trajectory or geometry.

## Material background

```bash
ISOTOPEFILE=${path_to_RelicsSim}/macro/isotopes/Kr85.mac ${path_to_RelicsSim}/build/BambooMC -c ${path_to_RelicsSim}/config/cevns.xml -m ${path_to_RelicsSim}/macro/components/lXe.mac [ -n 100 ] [ -o material.root ] [ -i ]
```

## CRN
For CRN, you need make sure you generator in first step is right. 
```bash
${path_to_RelicsSim}/build/BambooMC -c ${path_to_RelicsSim}/config/cevns.xml [ -n 100 ] [ -o muon.root ] [ -i ]
```

# 4 Four step, Sorting and Normalization
In ${path_to_RelicsSim}/Multi/energy, have a file names Sorting.py
```python
python3 ${path_to_RelicsSim}/Multi/energy/Sorting.py --InputFile ${path_to_RelicsSim}/build/muon.root --OutputFile ${path_to_RelicsSim}/build/muon.h5
```

# Multi-threaded simulation
In ${path_to_RelicsSim}, run code.
```python
python3 ./scripts/relics_xml_first.py ./config/geo_params.json -o ./config/cevns.xml  --gen muon  --check_overlap 
```

```bash
mkdir build
cd build
cmake -DENABLE_DETECTOR_SETS=relics -DENABLE_USER_MC=pandax,relics -DCMAKE_BUILD_TYPE=Debug ..
make
```

Attention! Note that you give .sh permission to run the file
```bash
cd ..
cd Multi
cd energy
./runall.sh -j 60 -f 1000 -e 10000 > runall.file
```
./runall.sh, you can choose which you want to run, muon or material 
-j, means Number of threads
-f, file number 1000, suggestion
-e, event number 10000, suggestion

