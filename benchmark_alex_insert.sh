#!/bin/bash

# chmod +x benchmark_alex_insert.sh

# nohup ./benchmark_alex_insert.sh &
# ./benchmark_alex_insert.sh
# ps aux | head -n 1 && ps aux | grep ./benchmark_alex_insert
# ps aux | head -n 1 && ps aux | grep ./alex_poi2

# actual paras
datasets=("fb" "covid" "osm" "genome")
# datasets=("fb")
# datasets=("test" "test")
# datasets=("genome_25M" "genome_12.5M" "genome_50M" "genome_100M" )
poisons=("1" "0")
# poison_sizes=("0.05" "0.1" "0.2" "0.4" "0.8")
poison_sizes=("0.1")
insert_props=("1")
cutoffs=(-30 -100 -50 -50)

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
g++ alex_poi2.cpp -std=c++17 -o alex_poi2 -march=native -mpopcnt

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
    exit 1
fi
index=0
# Iterate over the datasets and run the program
for dataset in "${datasets[@]}"; do
    cutoff=${cutoffs[$index]}
    for insert_prop in "${insert_props[@]}"; do
        for poison_size in "${poison_sizes[@]}"; do
            for poison in "${poisons[@]}"; do
                #echo "Dataset : Poison: Poison_size: Insert_prop"
                echo "Running with dataset: $dataset: poison $poison: poison size $poison_size: insert prop $insert_prop"
                #echo "=============="
                ./alex_poi2 "$dataset" 200000000 "$poison" "$poison_size" "$insert_prop" "$cutoff"
            done
        done
    done 
    index=$((index + 1))
done