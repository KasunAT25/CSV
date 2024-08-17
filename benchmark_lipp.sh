#!/bin/bash

# chmod +x benchmark_lipp.sh
# nohup ./benchmark_lipp.sh &

datasets=("fb" "covid" "osm" "genome")
# datasets=("test")

smooths=("1" "0")
smooth_sizes=("0.05" "0.1" "0.2" "0.4" "0.8")

insert_props=("0")

# Compile the C++ program
# g++ lipp_csv.cpp -std=c++17 -o lipp_csv -march=native -mpopcnt

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
        for smooth_size in "${smooth_sizes[@]}"; do
            for smooth in "${smooths[@]}"; do
                echo "Running with dataset: $dataset: smooth $smooth: smooth size $smooth_size: insert prop $insert_prop"
                #echo "=============="
                ./lipp_csv "$dataset" 200000000 "$smooth" "$smooth_size" "$insert_prop"
            done
        done
    done 
done