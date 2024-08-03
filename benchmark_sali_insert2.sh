#!/bin/bash

# chmod +x benchmark_sali_insert2.sh

# nohup ./benchmark_sali_insert2.sh &
# ./benchmark_sali_insert2.sh

# ps aux | head -n 1 && ps aux | grep ./benchmark_sali_insert2
# ps aux | head -n 1 && ps aux | grep ./sali_poi3

# actual paras
#datasets=("genome")
datasets=("fb" "covid" "osm" "genome")

#datasets=("uniform")
poisons=("1" "0")
poison_sizes=("20000000")
#poison_sizes=("60000000" "80000000")
insert_props=("1")

# testing paras
# datasets=("fb_200M_uint64")
# poisons=("1")
# poison_sizes=("10000")

# inserts=("1")
# insert_props=("0.2")

# format
#argv[1] data
#argv[2] size
#argv[3] poison or not argv[3][0]
#argv[4] poison size
#argv[5] insert prop

# datasets=("fb_200M_uint64" "osm_cellids_200M_uint64")

data_folder="../../../ext/Data"

# Compile the C++ program
g++ -fopenmp -I /opt/intel/oneapi/tbb/2021.12/include sali_poi3.cpp -std=c++17 -o sali_poi3 -L /opt/intel/oneapi/tbb/2021.12/lib -ltbb -march=native -mpopcnt

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
    exit 1
fi

# Iterate over the datasets and run the program
for dataset in "${datasets[@]}"; do
    for insert_prop in "${insert_props[@]}"; do
        for poison_size in "${poison_sizes[@]}"; do
            for poison in "${poisons[@]}"; do
                #echo "Dataset : Poison: Poison_size: Insert_prop"
                echo "Running with dataset: $dataset: poison $poison: poison size $poison_size: insert prop $insert_prop"
                #echo "=============="
                LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi3 "$dataset" 200000000 "$poison" "$poison_size" "$insert_prop"
            done
        done
    done 
done