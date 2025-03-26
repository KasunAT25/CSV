#pragma once

#include <iostream>
#include "./competitors/LIPP/src/core/lipp.h"

// void flushCache() {
//     constexpr size_t cacheSize = 32 * 1024 * 1024; // Adjust size as needed
//     static std::vector<char> dummyData(cacheSize, 0);
//     for (volatile char &val : dummyData) {
//         val = 0;
//     }
// }

long benchmark_lookup_lipp(std::vector<KEY_TYPE>& lookups, LIPP<KEY_TYPE, PAYLOAD_TYPE> & index){
   
    int size = lookups.size();
    auto start = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < size; i++){
         
         auto payload = index.at(lookups[i]);
         //assert(index.exists(lookup));
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    return duration.count();
}

std::vector<long> perform_benchmark(LIPP<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE>& data, std::vector<KEY_TYPE>& lookups, int num){
    std::vector<long> measurements;
    
    for (int i = 0; i < num; i++){
        flushCache();
        long time = benchmark_lookup_lipp(lookups, index);
        measurements.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());

    return measurements;
}


void benchmark_lipp_real(LIPP<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile,
bool poison, bool insert, KEY_TYPE orignal_P, double insert_threshold){


    int reps = 25;
    int size_data = lookups.size();

    std::vector<long> measurements = perform_benchmark(index, data, lookups, reps);

    double mean = 0;
    double median = measurements[measurements.size()/2];
    for(int i = 0; i < measurements.size(); i++){
        mean += measurements[i];
    }
    mean /= measurements.size();

    mean /= size_data;
    median /= size_data;

    long log_error = -1;
    long d_log_error = -1;
    long mse_error = -1;

     std::string model_results = "LIPP;" + data_name + ";" + std::to_string(data.size()) + ";"+ std::to_string(lookups.size()) + ";" +
          std::to_string(mean) + ";" +std::to_string(median) +";"+ std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) ;


    std::string outputfile = outfile+"_"+data_name+"_ind_"+std::to_string(poison) +"_" + std::to_string(insert) + "_" +
          std::to_string(orignal_P) +"_" + std::to_string(insert_threshold) +".csv";


    saveToCSV(model_results,outfile);
   
}