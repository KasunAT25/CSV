#pragma once

#include <iostream>
#include <omp.h>
#include "./competitors/SALI/src/core/sali.h"

//#define KEY_TYPE double
//#define PAYLOAD_TYPE double
//typedef alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> data_node_type_alex;



/*
    Calculates the error for a single element for a certain linear function
*/
// template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
// inline double calculate_error_single_element_lipp(std::vector<KEY_TYPE>& data, LIPP<KEY_TYPE, PAYLOAD_TYPE> & index, int i){
//     double x = data[i];
//     double y = i;
//     auto it = index.lower_bound(x);
//     return 0;
//     // Cannot calculate error since lower_bound always return correct position
//     //return apply_errorfn<E, ROUND, BOUNDED>(low_pos, y, data.size()-1);
// }

// /*
//     Calculates the total error for a linear function
// */
// template<ERROR_TYPE E, bool CORRECT = true, bool BOUNDED = true>
// long double calculate_error_lipp(std::vector<KEY_TYPE>& data, LIPP<KEY_TYPE, PAYLOAD_TYPE> & index) {
//     long double total_error = 0;
//     for (long i = 0; i < data.size(); i++) {
//         total_error += calculate_error_single_element_alex<E, CORRECT, BOUNDED>(data, index, i);
//     }
//     return total_error;
// }

/*
    Benchmark ALEX by measuring lookup times
*/
// Add the time codes here
// namespace {
//     static int size_data = 0; // Static variable limited to this translation unit
//     static std::vector<double> times_total;
//     static std::vector<double> times_query;
//     static std::vector<double> times_search;
//     static std::vector<KEY_TYPE> keys;
//     static std::vector<short> level;
//     static std::vector<int> diffs;
// }
long benchmark_lookup_sali(KEY_TYPE lookup, sali::SALI<KEY_TYPE, PAYLOAD_TYPE> & index, int num, int idx){
    //Add the vectors
    //std::cout << "In " << num << std::endl;
    //int size = lookups.size();
    //times_total.assign(size,0);
    

    // times_query.assign(size,0);
    // times_search.assign(size,0);
    // std::fill(times_total.begin(), times_total.end(), 0);
    // std::fill(times_query.begin(), times_query.end(), 0);
    // std::fill(times_search.begin(), times_search.end(), 0);

    //int data_size = data.size()-1;
    auto start = std::chrono::high_resolution_clock::now();
    PAYLOAD_TYPE  p;
    for(int i = 0; i < num; i++){
         
         auto payload = index.at(lookup, p);
         //assert(index.exists(lookup));

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    times_total[idx] = duration.count()/(1.0*num);
    
    keys[idx] = lookup;

    // std::transform(times_search.begin(), times_search.end(), times_query.begin(), times_total.begin(),
    //                    [](int a, int b) { return a + b; });

    return duration.count();
}

std::vector<long> perform_benchmark(sali::SALI<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE>& data, std::vector<KEY_TYPE>& lookups, int num){
    std::vector<long> measurements;
    std::fill(times_total.begin(), times_total.end(), 0);
    std::fill(times_query.begin(), times_query.end(), 0);
    std::fill(times_search.begin(), times_search.end(), 0);
    //std::fill(level.begin(), level.end(), 0);

    //std::cout << "Out 1" << std::endl;

     for (int i = 0; i < lookups.size(); i++){
        long time = benchmark_lookup_sali(lookups[i], index, num, i);
        measurements.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());

    // //Taking the average of the times (dividing by replicates)
    // std::transform(times_query.begin(), times_query.end(), times_query.begin(),
    //                [num](double x) { return x / num; });

    // std::transform(times_total.begin(), times_total.end(), times_total.begin(),
    //                [num](double x) { return x / num; });

    // std::transform(times_search.begin(), times_search.end(), times_search.begin(),
    //                [num](double x) { return x / num; });

    return measurements;
}

// void benchmark_lipp(std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile){

//     // Combine pre-defined keys with randomly generated payloads (can not pass a simple vector to ALEX)
//     auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()];
//     std::mt19937_64 gen_payload(std::random_device{}());
//     for (int i = 0; i < data.size(); i++) {
//         values[i].first = data[i];
//         values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
//     }
//     auto start = std::chrono::high_resolution_clock::now();

//     // Create ALEX and bulk load
//     alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
//     index.bulk_load(values, data.size());

//     auto stop = std::chrono::high_resolution_clock::now();
//     long build_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

//     std::vector<long> measurements = perform_benchmark(index, data, lookups, 30);

//     double mean = 0;
//     double median = measurements[measurements.size()/2];
//     for(int i = 0; i < measurements.size(); i++){
//         mean += measurements[i];
//     }
//     mean /= measurements.size();

//     mean /= lookups.size();
//     median /= lookups.size();

//     long log_error = -1;
//     long d_log_error = -1;
//     long mse_error = -1;

//     double average_search = static_cast<double>(std::accumulate(times_search.begin(), times_search.end(), 0)) / times_search.size();
//     double average_query = static_cast<double>(std::accumulate(times_query.begin(), times_query.end(), 0)) / times_query.size();
//     double average_tot = static_cast<double>(std::accumulate(times_total.begin(), times_total.end(), 0)) / times_total.size();

//     std::ofstream file;
//     file.open(outfile+".txt", std::ios_base::app);
//     //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
//     file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time <<  ";"  << average_tot << ";" << average_query << ";" << average_search <<std::endl;

//     file.close();

//     std::cout << regression_name << " data_name:" << data_name << "; poisoning_threshold: " << poisoning_threshold << "; data size:" << data.size() << " lookups size:" << lookups.size() << " mean lookup ns:" << mean << " median lookup ns:" << median << " log error:" << log_error << " discrete log error:" << d_log_error << " mse error:" << mse_error << " build time:" << build_time << std::endl;

//     std::cout << regression_name << " data_name:" << data_name << "; Time Total: " << average_tot <<"; Time Query: " << average_query << "; Time Search:" << average_search <<  std::endl;

//     std::string outputfile = outfile+"_ALEX_detailed.csv";

//     //save_to_csv(times_total,times_query,times_search,headers,outputfile);


//     delete[] values;
// }

void benchmark_sali_real(sali::SALI<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile,
bool poison, bool insert, KEY_TYPE orignal_P, double insert_threshold){

    // Combine pre-defined keys with randomly generated payloads (can not pass a simple vector to ALEX)
    // auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()];
    // std::mt19937_64 gen_payload(std::random_device{}());
    // for (int i = 0; i < data.size(); i++) {
    //     values[i].first = data[i];
    //     values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    // }
    auto start = std::chrono::high_resolution_clock::now();

    // // Create ALEX and bulk load
    // alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    // index.bulk_load(values, data.size());

    auto stop = std::chrono::high_resolution_clock::now();
    long build_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

    int reps = 100;
    size_data = lookups.size();
    times_total.resize(size_data); // Resize the vector
    std::fill(times_total.begin(), times_total.end(), 0.0);

    times_query.resize(size_data); // Resize the vector
    std::fill(times_query.begin(), times_query.end(), 0.0);

    times_search.resize(size_data); // Resize the vector
    std::fill(times_search.begin(), times_search.end(), 0.0);

    keys.resize(size_data); // Resize the vector
    std::fill(keys.begin(), keys.end(), 0.0);

    level.resize(size_data); // Resize the vector
    std::fill(level.begin(), level.end(), 0.0);

    diffs.resize(size_data); // Resize the vector
    std::fill(diffs.begin(), diffs.end(), 0.0);

    std::vector<long> measurements = perform_benchmark(index, data, lookups, reps);

    double mean = 0;
    double median = measurements[measurements.size()/2];
    for(int i = 0; i < measurements.size(); i++){
        mean += measurements[i];
    }
    mean /= measurements.size();

    mean /= reps;
    median /= reps;

    long log_error = -1;
    long d_log_error = -1;
    long mse_error = -1;

    double average_search = static_cast<double>(std::accumulate(times_search.begin(), times_search.end(), 0)) / lookups.size();
    double average_query = static_cast<double>(std::accumulate(times_query.begin(), times_query.end(), 0)) / lookups.size();
    double average_tot = static_cast<double>(std::accumulate(times_total.begin(), times_total.end(), 0)) / lookups.size();

    // std::ofstream file;
    // file.open(outfile+".txt", std::ios_base::app);
    // //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
    // // file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time <<  ";"  << average_tot << ";" << average_query << ";" << average_search <<std::endl;
    // file << regression_name << ";" << data_name << ";"  << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" <<std::endl;

    // file.close();

     std::string model_results = "SALI;" + data_name + ";" + std::to_string(data.size()) + ";"+ std::to_string(lookups.size()) + ";" +
          std::to_string(mean) + ";" +std::to_string(median) +";"+ std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) ;

    std::cout << regression_name << " data_name:" << data_name << "; poisoning_threshold: " << poisoning_threshold << "; data size:" << data.size() << " lookups size:" << lookups.size() << " mean lookup ns:" << mean << " median lookup ns:" << median << " log error:" << log_error << " discrete log error:" << d_log_error << " mse error:" << mse_error << " build time:" << build_time << std::endl;

    std::cout << regression_name << " data_name:" << data_name << "; Time Total: " << average_tot <<"; Time Query: " << average_query << "; Time Search:" << average_search <<  std::endl;

    std::cout << "MEAN: "<< mean << std::endl;

   std::string outputfile = outfile+"_"+data_name+"_ind_"+std::to_string(poison) +"_" + std::to_string(insert) + "_" +
          std::to_string(orignal_P) +"_" + std::to_string(insert_threshold) +".csv";

    std::vector<std::string> headers2 = {"Key", "Time Total","Query Time", "Search Time", "Level"};

     saveToCSV(model_results,outfile);

    //save_level_to_csv(keys,times_total,times_query,times_search,level,headers2,outputfile);
    //save_to_csv(keys,times_total,times_query,times_search,headers,outputfile);

    std::vector<std::string> headers = {"Key", "Time Total"};
    //save_level_to_csv(keys,times_total,times_query,times_search,level,headers2,outputfile);
    // save_to_csv2(keys,times_total,headers,outputfile);


    //delete[] values;
}