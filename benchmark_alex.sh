#!/bin/bash

# chmod +x benchmark_alex.sh



datasets=("fb" "covid" "osm" "genome")

smooths=("1" "0")
# smooth_sizes=("0.1")
smooth_sizes=("0.05" "0.1" "0.2" "0.4" "0.8")
insert_props=("0")
cutoffs=(-50 -50 -50 -50)

# Compile the C++ program
# g++ alex_csv.cpp -std=c++17 -o alex_csv -march=native -mpopcnt
g++ -pthread alex_csv_par.cpp -std=c++17 -o alex_csv_par -march=native -mpopcnt

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
                echo "Running with dataset: $dataset: smooth $smooth: smooth size $smooth_size: insert prop $insert_prop"
                #echo "=============="
                ./alex_csv_par "$dataset" 200000000 "$smooth" "$smooth_size" "$insert_prop" "$cutoff"
            done
        done
    done 
    index=$((index + 1))
done