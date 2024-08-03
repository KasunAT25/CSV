#!/bin/bash

# chmod +x benchmark_lipp_insert.sh
# nohup ./benchmark_lipp_insert.sh &
# ps aux | head -n 1 && ps aux | grep ./benchmark_lipp
# ps aux | head -n 1 && ps aux | grep ./lipp_poi
# actual paras
datasets=("fb" "covid" "osm" "genome")
# datasets=("genome")
#datasets=("uniform")
poisons=("1" "0")
poison_sizes=("20000000")
insert_props=("1")

# testing paras

# datasets=("fb_200M_uint64")
# poisons=("1")
# poison_sizes=("10000")
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
g++ lipp_poi2.cpp -std=c++17 -o lipp_poi2 -march=native -mpopcnt

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
                # ./lipp_poi "$data_folder/$dataset" 100000 "$poison" "$poison_size" "$insert_prop"
                ./lipp_poi2 "$dataset" 200000000 "$poison" "$poison_size" "$insert_prop"
                 # nohup ./lipp_poi "$dataset" 200000000 "$poison" "$poison_size" "$insert_prop" &
            done
        done
    done 
done