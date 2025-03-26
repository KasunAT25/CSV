#pragma once

#include <iostream>
#include "./competitors/ALEX/src/core/alex.h"

typedef alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> data_node_type_alex;


void flushCache() {
    constexpr size_t cacheSize = 32 * 1024 * 1024; // Adjust size as needed
    static std::vector<char> dummyData(cacheSize, 0);
    for (volatile char &val : dummyData) {
        val = 0;
    }
}
/*
    Calculates the error for a single element for a certain linear function
*/
template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
inline double calculate_error_single_element_alex(std::vector<KEY_TYPE>& data, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, int i){
    double x = data[i];
    double y = i;
    auto it = index.lower_bound(x);
    return 0;
    // Cannot calculate error since lower_bound always return correct position
    //return apply_errorfn<E, ROUND, BOUNDED>(low_pos, y, data.size()-1);
}

/*
    Calculates the total error for a linear function
*/
template<ERROR_TYPE E, bool CORRECT = true, bool BOUNDED = true>
long double calculate_error_alex(std::vector<KEY_TYPE>& data, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index) {
    long double total_error = 0;
    for (long i = 0; i < data.size(); i++) {
        total_error += calculate_error_single_element_alex<E, CORRECT, BOUNDED>(data, index, i);
    }
    return total_error;
}

/*
    Benchmark ALEX by measuring lookup times
*/

long benchmark_lookup_alex(std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE>& lookups, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index){
   
    int size = lookups.size();
    auto start = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < size; i++){
        
        data_node_type_alex* leaf = index.get_leaf(lookups[i]);
        //std::cout << leaf << std::endl;
        int location_key_real = leaf->find_key(lookups[i]);
        
    }
    
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    return duration.count();
}

long benchmark_lookup_alex_query(std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE>& lookups, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index){

    int size = lookups.size();

    auto start = std::chrono::high_resolution_clock::now();
    
    for(int i = 0; i < size; i++){
         data_node_type_alex* leaf = index.get_leaf(lookups[i]);
    }   

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    return duration.count();
}

std::vector<long> perform_benchmark(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE>& data, std::vector<KEY_TYPE>& lookups, int num){
    std::vector<long> measurements;
   
    for (int i = 0; i < num; i++){
        flushCache();
        long time = benchmark_lookup_alex(data, lookups, index);
        measurements.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());

    return measurements;
}

std::pair<std::vector<long>, std::vector<long>> perform_benchmark_seperate(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE>& data, std::vector<KEY_TYPE>& lookups, int num){
    std::vector<long> measurements;
    std::vector<long> measurements_q;

    for (int i = 0; i < num; i++){
        flushCache();
        long time = benchmark_lookup_alex(data, lookups, index);
        measurements.push_back(time);
    }
    //get just the query times
    for (int i = 0; i < num; i++){
        flushCache();
        long time = benchmark_lookup_alex_query(data, lookups, index);
        measurements_q.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());
    std::sort(measurements_q.begin(), measurements_q.end());

    return std::make_pair(measurements, measurements_q);
}


void benchmark_alex_real(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile){

    
    int size_data = lookups.size();

    std::vector<long> measurements = perform_benchmark(index, data, lookups, 30);

    double mean = 0;
    double median = measurements[measurements.size()/2];
    for(int i = 0; i < measurements.size(); i++){
        mean += measurements[i];
    }
    mean /= measurements.size();

    mean /= size_data;
    median /= size_data;

    long log_error = calculate_error_alex<LogNorm, true, false>(data, index);
    long d_log_error = calculate_error_alex<DiscreteLogNorm, true, false>(data, index);
    long mse_error = calculate_error_alex<L2Norm, true, false>(data, index);


    std::cout << regression_name << " data_name:" << data_name << "; poisoning_threshold: " << poisoning_threshold << "; data size:" << data.size() << " lookups size:" << lookups.size() << " mean lookup ns:" << mean << " median lookup ns:" << median << " log error:" << log_error << " discrete log error:" << d_log_error << " mse error:" << mse_error  << std::endl;


    std::string outputfile = outfile+"_ALEX_detailed.csv";

    std::cout << "MEAN: "<< mean << std::endl;


}

void benchmark_alex_real_seperate(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile,
bool poison, bool insert, std::string orignal_P, double insert_threshold){

    int reps = 25;
    int size_data = lookups.size();

    auto start = std::chrono::high_resolution_clock::now();

    auto stop = std::chrono::high_resolution_clock::now();
    long build_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

    std::pair<std::vector<long>, std::vector<long>> measurements_all = perform_benchmark_seperate(index, data, lookups, reps);
    std::vector<long> measurements = measurements_all.first;
    std::vector<long> measurements_q = measurements_all.second;

    double mean = 0;
    double median = measurements[measurements.size()/2];
    for(int i = 0; i < measurements.size(); i++){
        mean += measurements[i];
    }
    mean /= measurements.size();

    mean /= size_data;
    median /= size_data;


    std::string outputfile = outfile+"_"+data_name+"_ind_"+std::to_string(poison) +"_" + std::to_string(insert) + "_" +
          (orignal_P) +"_" + std::to_string(insert_threshold) +".csv";


    double mean_q = 0;
    double median_q = measurements_q[measurements_q.size()/2];
    for(int i = 0; i < measurements_q.size(); i++){
        mean_q += measurements_q[i];
    }
    mean_q /= measurements_q.size();

    mean_q /= size_data;
    median_q /= size_data;

    std::cout << "Total Mean: "<< mean << std::endl;

    std::cout << "Query Mean: "<< mean_q << std::endl;

    std::cout << "Search Mean: "<< (mean - mean_q) << std::endl;


    std::string model_results = "ALEX;" + data_name + ";" + std::to_string(data.size()) + ";"+ std::to_string(lookups.size()) + ";" +
          std::to_string(mean) + ";" +std::to_string(median) +";"+ std::to_string(mean_q) + ";" +std::to_string(median_q) +";"+
          std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          (orignal_P)  + ";"  + std::to_string(insert_threshold) ;

    saveToCSV(model_results,outfile);
   
}