#pragma once

#include <iostream>
#include "./competitors/ALEX/src/core/alex.h"

//#define KEY_TYPE double
//#define PAYLOAD_TYPE double
typedef alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> data_node_type_alex;

//int size = 1000000;
namespace {
    static int size_data = 0; // Static variable limited to this translation unit
    static std::vector<double> times_total;
    static std::vector<double> times_query;
    static std::vector<double> times_search;
    static std::vector<KEY_TYPE> keys;
    static std::vector<short> level;
    static std::vector<int> diffs;
}
std::vector<std::string> headers = {"Key","Time Total","Query Time", "Search Time"};
// std::vector<double> times_total(size_data,0);
// std::vector<double> times_query(size_data,0);
// std::vector<double> times_search(size_data,0);

// std::vector<KEY_TYPE> keys(size_data,0);
// std::vector<short> level(size_data,0);


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
// Add the time codes here

long benchmark_lookup_alex(std::vector<KEY_TYPE> & data, KEY_TYPE lookup, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, int num, int idx){
    //Add the vectors
    //std::cout << "In " << num << std::endl;
    //int size = lookups.size();
    //times_total.assign(size,0);
    

    // times_query.assign(size,0);
    // times_search.assign(size,0);
    // std::fill(times_total.begin(), times_total.end(), 0);
    // std::fill(times_query.begin(), times_query.end(), 0);
    // std::fill(times_search.begin(), times_search.end(), 0);

    int data_size = data.size()-1;
    auto start = std::chrono::high_resolution_clock::now();
    //int diffs = 0;
    for(int i = 0; i < num; i++){
        //DoNotOptimize(index.lower_bound(lookups[i]));
        //auto start_q = std::chrono::high_resolution_clock::now();

        //auto payload = index.get_payload(lookups[i]);
        // auto payload = index.find(lookups[i]);

         data_node_type_alex* leaf = index.get_leaf(lookup);

         //auto end_q = std::chrono::high_resolution_clock::now();

        //auto start_search = std::chrono::high_resolution_clock::now();

        //  std::cout << "Leaf :"<< leaf << std::endl;
        //   std::cout << "Lookup :"<< lookups[i] << std::endl;

         int location_key_real = leaf->find_key(lookup);
        // int location_key = leaf->predict_position(lookups[i]);
        // diffs+= pow(location_key_real-location_key,2);
        //auto end_search = std::chrono::high_resolution_clock::now();


        
            
        


        //int predicted_pos = index.predict(lookups[i]);
         //auto end_q = std::chrono::high_resolution_clock::now();
         //times_query[i] += std::chrono::duration_cast<std::chrono::nanoseconds>(end_q - start_q).count();
         //times_search[i] += std::chrono::duration_cast<std::chrono::nanoseconds>(end_search - start_search).count();

    //      if(num ==1){
    //         //level[i] = leaf->level_;
    //         //keys[i] = lookups[i];
    //         //level[i] = 1;
    //         diffs+= leaf->number_of_searches_;
    //        //std::cout << "Key :"<< lookups[i] << " Found :" << location_key_real <<  " | ";
    //        //std::cout << "Key :"<< lookups[i]  <<" Diff :" << location_key_real-location_key<<  " | ";
    //    }
        if(i==0){
            level[idx] = leaf->level_;
            diffs[idx] = leaf->number_of_searches_;
        }

        // if(lookups[i] !=leaf->key_slots_[location_key]){
        //     std::cout << "THIS DOES NOT WORK" << std::endl;
        //     break;
        // }
        // if(lookup !=leaf->key_slots_[location_key_real]){
        //     std::cout << "Error "<< " | ";
        // }
        //std::cout << lookups[i] <<" : " << leaf->key_slots_[location_key_real] << " , "<< std::endl;
    }
    
    auto stop = std::chrono::high_resolution_clock::now();

    // if(num == 1){
    //     std::cout << "Diffs :"<<  diffs << std::endl;
    // }
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    times_total[idx] = duration.count()/(1.0*num);
    
    keys[idx] = lookup;

    
            // level[idx] = 1;
            // diffs[idx] = 1;
    

    // std::transform(times_search.begin(), times_search.end(), times_query.begin(), times_total.begin(),
    //                    [](int a, int b) { return a + b; });

    return duration.count();
}

long benchmark_lookup_alex_query(std::vector<KEY_TYPE> & data, KEY_TYPE lookup, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, int num, int idx){
    //Add the vectors
    //std::cout << "In " << num << std::endl;
    //int size = lookups.size();
    //times_total.assign(size,0);
    

    // times_query.assign(size,0);
    // times_search.assign(size,0);
    // std::fill(times_total.begin(), times_total.end(), 0);
    // std::fill(times_query.begin(), times_query.end(), 0);
    // std::fill(times_search.begin(), times_search.end(), 0);

    int data_size = data.size()-1;
    int diffs = 0;
    auto start = std::chrono::high_resolution_clock::now();
    
    for(int i = 0; i < num; i++){

         data_node_type_alex* leaf = index.get_leaf(lookup);

    }   

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    // std::transform(times_search.begin(), times_search.end(), times_query.begin(), times_total.begin(),
    //                    [](int a, int b) { return a + b; });

    times_query[idx] = duration.count()/(1.0*num);

    return duration.count();
}

std::vector<long> perform_benchmark(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE>& data, std::vector<KEY_TYPE>& lookups, int num){
    std::vector<long> measurements;
    std::fill(times_total.begin(), times_total.end(), 0);
    std::fill(times_query.begin(), times_query.end(), 0);
    std::fill(times_search.begin(), times_search.end(), 0);
    std::fill(level.begin(), level.end(), 0);

    //std::cout << "Out 1" << std::endl;

    for (int i = 0; i < lookups.size(); i++){
        long time = benchmark_lookup_alex(data, lookups[i], index,num, i);
        measurements.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());



    //Taking the average of the times (dividing by replicates)
    // std::transform(times_query.begin(), times_query.end(), times_query.begin(),
    //                [num](double x) { return x / num; });

    // std::transform(times_total.begin(), times_total.end(), times_total.begin(),
    //                [num](double x) { return x / num; });

    // std::transform(times_search.begin(), times_search.end(), times_search.begin(),
    //                [num](double x) { return x / num; });

    return measurements;
}

std::pair<std::vector<long>, std::vector<long>> perform_benchmark_seperate(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE>& data, std::vector<KEY_TYPE>& lookups, int num){
    std::vector<long> measurements;
    std::vector<long> measurements_q;

    std::fill(times_total.begin(), times_total.end(), 0);
    std::fill(times_query.begin(), times_query.end(), 0);
    std::fill(times_search.begin(), times_search.end(), 0);
    std::fill(level.begin(), level.end(), 0);

    //std::cout << "Out 1" << std::endl;
    //get the total time
    for (int i = 0; i < lookups.size(); i++){
        long time = benchmark_lookup_alex(data, lookups[i], index,num, i);
        measurements.push_back(time);
    }
    //get just the query times
    for (int i = 0; i < lookups.size(); i++){
        long time = benchmark_lookup_alex_query(data, lookups[i], index,num, i);
        measurements_q.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());
    std::sort(measurements_q.begin(), measurements_q.end());

    std::transform(times_total.begin(), times_total.end(), times_query.begin(), times_search.begin(),
                       [](double a, double b) { return a - b; });

    //Taking the average of the times (dividing by replicates)
    // std::transform(times_query.begin(), times_query.end(), times_query.begin(),
    //                [num](double x) { return x / num; });

    // std::transform(times_total.begin(), times_total.end(), times_total.begin(),
    //                [num](double x) { return x / num; });

    // std::transform(times_search.begin(), times_search.end(), times_search.begin(),
    //                [num](double x) { return x / num; });

    return std::make_pair(measurements, measurements_q);
}

// void benchmark_alex(std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile){

//     // Combine pre-defined keys with randomly generated payloads (can not pass a simple vector to ALEX)
//     auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()];
//     std::mt19937_64 gen_payload(std::random_device{}());
//     //Try sending the values themselves for testing next
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

//     long log_error = calculate_error_alex<LogNorm, true, false>(data, index);
//     long d_log_error = calculate_error_alex<DiscreteLogNorm, true, false>(data, index);
//     long mse_error = calculate_error_alex<L2Norm, true, false>(data, index);

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

void benchmark_alex_real(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile){

    // Combine pre-defined keys with randomly generated payloads (can not pass a simple vector to ALEX)
    // auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()];
    // std::mt19937_64 gen_payload(std::random_device{}());
    // for (int i = 0; i < data.size(); i++) {
    //     values[i].first = data[i];
    //     values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    // }
    size_data = lookups.size();

    auto start = std::chrono::high_resolution_clock::now();

    // // Create ALEX and bulk load
    // alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    // index.bulk_load(values, data.size());

    auto stop = std::chrono::high_resolution_clock::now();
    long build_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

    std::vector<long> measurements = perform_benchmark(index, data, lookups, 30);

    double mean = 0;
    double median = measurements[measurements.size()/2];
    for(int i = 0; i < measurements.size(); i++){
        mean += measurements[i];
    }
    mean /= measurements.size();

    mean /= lookups.size();
    median /= lookups.size();

    long log_error = calculate_error_alex<LogNorm, true, false>(data, index);
    long d_log_error = calculate_error_alex<DiscreteLogNorm, true, false>(data, index);
    long mse_error = calculate_error_alex<L2Norm, true, false>(data, index);

    double average_search = static_cast<double>(std::accumulate(times_search.begin(), times_search.end(), 0)) / lookups.size();
    double average_query = static_cast<double>(std::accumulate(times_query.begin(), times_query.end(), 0)) / lookups.size();
    double average_tot = static_cast<double>(std::accumulate(times_total.begin(), times_total.end(), 0)) / lookups.size();

    // std::ofstream file;
    // file.open(outfile+".txt", std::ios_base::app);
    // //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
    // file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time <<  ";"  << average_tot << ";" << average_query << ";" << average_search <<std::endl;

    // file.close();

    std::cout << regression_name << " data_name:" << data_name << "; poisoning_threshold: " << poisoning_threshold << "; data size:" << data.size() << " lookups size:" << lookups.size() << " mean lookup ns:" << mean << " median lookup ns:" << median << " log error:" << log_error << " discrete log error:" << d_log_error << " mse error:" << mse_error << " build time:" << build_time << std::endl;

    std::cout << regression_name << " data_name:" << data_name << "; Time Total: " << average_tot <<"; Time Query: " << average_query << "; Time Search:" << average_search <<  std::endl;

    std::string outputfile = outfile+"_ALEX_detailed.csv";

    //std:: cout << "==============" <<std::endl;
    std::cout << "MEAN: "<< mean << std::endl;

    std::vector<std::string> headers2 = {"Key", "Time Total","Query Time", "Search Time", "Level"};

    //save_level_to_csv(keys,times_total,times_query,times_search,level,headers2,outputfile);
    //save_to_csv(keys,times_total,times_query,times_search,headers,outputfile);


    //delete[] values;
}

void benchmark_alex_real_seperate(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> & index, std::vector<KEY_TYPE> & data, std::vector<KEY_TYPE> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile,
bool poison, bool insert, std::string orignal_P, double insert_threshold){

    // Combine pre-defined keys with randomly generated payloads (can not pass a simple vector to ALEX)
    // auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()];
    // std::mt19937_64 gen_payload(std::random_device{}());
    // for (int i = 0; i < data.size(); i++) {
    //     values[i].first = data[i];
    //     values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    // }
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


    // times_query(size_data,0);
    // times_search(size_data,0);

    // keys(size_data,0);
    // level(size_data,0);

    //std::cout << size_data << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // // Create ALEX and bulk load
    // alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    // index.bulk_load(values, data.size());

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

    mean /= reps;
    median /= reps;

    long log_error = calculate_error_alex<LogNorm, true, false>(data, index);
    long d_log_error = calculate_error_alex<DiscreteLogNorm, true, false>(data, index);
    long mse_error = calculate_error_alex<L2Norm, true, false>(data, index);

    // double average_search = static_cast<double>(std::accumulate(times_search.begin(), times_search.end(), 0)) / lookups.size();
    // double average_query = static_cast<double>(std::accumulate(times_query.begin(), times_query.end(), 0)) / lookups.size();
    // double average_tot = static_cast<double>(std::accumulate(times_total.begin(), times_total.end(), 0)) / lookups.size();

    // std::ofstream file;
    // file.open(outfile+".txt", std::ios_base::app);
    // //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
    // file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time <<  ";"  << average_tot << ";" << average_query << ";" << average_search <<std::endl;

    // file.close();

    //std::cout << regression_name << " data_name:" << data_name << "; poisoning_threshold: " << poisoning_threshold << "; data size:" << data.size() << " lookups size:" << lookups.size() << " mean lookup ns:" << mean << " median lookup ns:" << median << " log error:" << log_error << " discrete log error:" << d_log_error << " mse error:" << mse_error << " build time:" << build_time << std::endl;

    //std::cout << regression_name << " data_name:" << data_name << "; Time Total: " << average_tot <<"; Time Query: " << average_query << "; Time Search:" << average_search <<  std::endl;

    std::string outputfile = outfile+"_"+data_name+"_ind_"+std::to_string(poison) +"_" + std::to_string(insert) + "_" +
          (orignal_P) +"_" + std::to_string(insert_threshold) +".csv";

    //std:: cout << "==============" <<std::endl;
    

    double mean_q = 0;
    double median_q = measurements_q[measurements_q.size()/2];
    for(int i = 0; i < measurements_q.size(); i++){
        mean_q += measurements_q[i];
    }
    mean_q /= measurements_q.size();

    mean_q /= reps;
    median_q /= reps;

    std::cout << "Total Mean: "<< mean << std::endl;

    std::cout << "Query Mean: "<< mean_q << std::endl;

    std::cout << "Search Mean: "<< (mean - mean_q) << std::endl;

    std::vector<std::string> headers2 = {"Key", "Time Total","Query Time", "Search Time","Searches", "Level"};

    std::ofstream file;
    file.open(outfile+".txt", std::ios_base::app);
    //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
    // file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time <<  ";"  << average_tot << ";" << average_query << ";" << average_search <<std::endl;
    file << regression_name << ";" << data_name << ";"  << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << mean_q << ";" << median_q << ";" <<std::endl;

    file.close();

    std::string model_results = "ALEX;" + data_name + ";" + std::to_string(data.size()) + ";"+ std::to_string(lookups.size()) + ";" +
          std::to_string(mean) + ";" +std::to_string(median) +";"+ std::to_string(mean_q) + ";" +std::to_string(median_q) +";"+
          std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          (orignal_P)  + ";"  + std::to_string(insert_threshold) ;

    //save_level_to_csv(keys,times_total,times_query,times_search,diffs,level,headers2,outputfile);
    //save_to_csv(keys,times_total,times_query,times_search,headers,outputfile);
    saveToCSV(model_results,outfile);

    std::vector<std::string> headers = {"Key", "Time Total"};
    //save_to_csv2(keys,times_total,headers,outputfile);

    //delete[] values;
}