#!/bin/bash

# chmod +x benchmark_alex_insert.sh

# nohup ./benchmark_alex_insert.sh &
# ./benchmark_alex_insert.sh

datasets=("fb" "covid" "osm" "genome")
# datasets=("test")

smooths=("1" "0")
smooth_sizes=("0.1")
insert_props=("1")
cutoffs=(-30 -100 -50 -50)

# Compile the C++ program
# g++ alex_csv.cpp -std=c++17 -o alex_csv -march=native -mpopcnt

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
        for smooth_size in "${smooth_sizes[@]}"; do
            for smooth in "${smooths[@]}"; do
                echo "Running with dataset: $dataset: poison $smooth: poison size $smooth_size: insert prop $insert_prop"
                #echo "=============="
                ./alex_csv "$dataset" 200000000 "$smooth" "$smooth_size" "$insert_prop" "$cutoff"
            done
        done
    done 
    index=$((index + 1))
done