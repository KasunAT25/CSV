#!/bin/bash

# chmod +x benchmark_alex.sh

# nohup ./benchmark_alex.sh &
# ./benchmark_alex.sh
# ps aux | head -n 1 && ps aux | grep ./benchmark_alex
# ps aux | head -n 1 && ps aux | grep ./alex_poi

# actual paras
# datasets=("fb" "covid" "osm" "genome")
datasets=("osm" "genome")

# datasets=("genome_25M" "genome_12.5M" "genome_50M" "genome_100M" )
poisons=("1" "0")
# poison_sizes=("0.05" "0.1" "0.2" "0.4" "0.8")
poison_sizes=("0.8")
insert_props=("0")
dataset_int=(100000000 50000000 25000000 12500000)

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
g++ alex_poi.cpp -std=c++17 -o alex_poi -march=native -mpopcnt

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
                ./alex_poi "$dataset" 200000000 "$poison" "$poison_size" "$insert_prop"
            done
        done
    done 
done