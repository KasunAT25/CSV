#pragma once

#include <chrono>
#include <random>
#include <iostream>
#include <algorithm>
#include <string>
#include "../exponential_search.h"
#include "../error_function.h"
#include "../linear_regression.h"
#include "./io_handler.h"

template <class Tp>
inline void DoNotOptimize(Tp const& value) {
    asm volatile("" : : "r,m"(value) : "memory");
}
// void save_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename);

/*
    Benchmark Linear Regressions by measuring lookuptimes
*/
// Add the time codes here
int size = 1000000;
std::vector<std::string> headers = {"Time Total","Query Time", "Search Time"};
std::vector<double> times_total(size,0);
std::vector<double> times_query(size,0);
std::vector<double> times_search(size,0);


long benchmark_lookup_linear_regression(std::vector<double> & data, std::vector<double> lookups, linear_regression lr){ 
    //Add the vectors
    // int size = lookups.size();
    // // std::vector<double> times_total(size,0);
    // // std::vector<double> times_query(size,0);
    // // std::vector<double> times_search(size,0);
    // times_total.resize(size,0);
    // times_query.resize(size,0);
    // times_search.resize(size,0);
    //std::cout << "Benchmark lookup" << std::endl;

    int data_size = data.size()-1;
    auto start = std::chrono::high_resolution_clock::now();
    
    for(int i = 0; i < lookups.size(); i++){
        //Start time
        auto start_q = std::chrono::high_resolution_clock::now();
        long pos = lr.predict(lookups[i]);
        pos = std::max<long>(std::min<long>(pos, data_size), 0);
        auto end_q = std::chrono::high_resolution_clock::now();
        auto start_search = std::chrono::high_resolution_clock::now();
        DoNotOptimize(exponential_search_lower_bound_linear_head(pos, lookups[i], data));
        //End time
        auto end_search = std::chrono::high_resolution_clock::now();
        
        //std::cout << "Adding" << std::endl;
        times_query[i] += std::chrono::duration_cast<std::chrono::nanoseconds>(end_q - start_q).count();
        times_search[i] += std::chrono::duration_cast<std::chrono::nanoseconds>(end_search - start_search).count();
        //std::cout << "Added" << std::endl;
        
    }
    //times_search[i] += times_total[i]-times_query[i];

    
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    std::transform(times_search.begin(), times_search.end(), times_query.begin(), times_total.begin(),
                       [](int a, int b) { return a + b; });
   
    

    return duration.count();
}

typedef linear_regression (*build_function)(std::vector<double> &);
std::vector<long> perform_benchmark(linear_regression lr, std::vector<double>& data, std::vector<double>& lookups, int num){
    std::vector<long> measurements;
    int size = lookups.size();
    // std::vector<double> times_total(size,0);
    // std::vector<double> times_query(size,0);
    // std::vector<double> times_search(size,0);
    // times_total.resize(size,0);
    // times_query.resize(size,0);
    // times_search.resize(size,0);
   // std::cout << "Performing" << std::endl;
    std::fill(times_total.begin(), times_total.end(), 0);
    std::fill(times_query.begin(), times_query.end(), 0);
    std::fill(times_search.begin(), times_search.end(), 0);
    //std::cout << "Done" << std::endl;

    for (int i = 0; i < num; i++){
        long time = benchmark_lookup_linear_regression(data, lookups, lr);
        measurements.push_back(time);
    }

    std::sort(measurements.begin(), measurements.end());

    //Taking the average of the times (dividing by replicates)
    std::transform(times_query.begin(), times_query.end(), times_query.begin(),
                   [num](double x) { return x / num; });

    std::transform(times_total.begin(), times_total.end(), times_total.begin(),
                   [num](double x) { return x / num; });

    std::transform(times_search.begin(), times_search.end(), times_search.begin(),
                   [num](double x) { return x / num; });


    return measurements;
}


template<build_function f>
void benchmark_regression(std::vector<double> & data, std::vector<double> & lookups, std::string regression_name, std::string data_name, double poisoning_threshold, std::string outfile){
    //std::cout << "Starting" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    linear_regression lr = f(data);
    auto stop = std::chrono::high_resolution_clock::now();
    long build_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

    //std::cout << "Benchmarking" << std::endl;
    std::vector<long> measurements = perform_benchmark(lr, data, lookups, 30);

    double mean = 0;
    double median = measurements[measurements.size()/2];
    for(int i = 0; i < measurements.size(); i++){
        mean += measurements[i];
    }
    mean /= measurements.size();

    mean /= lookups.size();
    median /= lookups.size();

    long log_error = calculate_error<LogNorm, true, false>(data, lr);
    long d_log_error = calculate_error<DiscreteLogNorm, true, false>(data, lr);
    long mse_error = calculate_error<L2Norm, true, false>(data, lr);

    double average_search = static_cast<double>(std::accumulate(times_search.begin(), times_search.end(), 0)) / times_search.size();
    double average_query = static_cast<double>(std::accumulate(times_query.begin(), times_query.end(), 0)) / times_query.size();
    double average_tot = static_cast<double>(std::accumulate(times_total.begin(), times_total.end(), 0)) / times_total.size();

    

    std::ofstream file;
    //file.open(outfile+".csv", std::ios_base::app);
    file.open(outfile+".txt", std::ios_base::app);
    //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
    file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time <<  ";"  << average_tot << ";" << average_query << ";" << average_search <<std::endl;

    file.close();
    
    std::cout << regression_name << " data_name:" << data_name << "; poisoning_threshold: " << poisoning_threshold << "; data size:" << data.size() << " lookups size:" << lookups.size() << " mean lookup ns:" << mean << " median lookup ns:" << median << " log error:" << log_error << " discrete log error:" << d_log_error << " mse error:" << mse_error << " build time:" << build_time << std::endl;
    std::cout << regression_name << " data_name:" << data_name << "; Time Total: " << average_tot <<"; Time Query: " << average_query << "; Time Search:" << average_search <<  std::endl;

    //Write to file
    std::string outputfile = outfile+"_SLR_detailed.csv";
    save_to_csv(times_total,times_query,times_search,headers,outputfile);

    

}

// void save_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename) {
    
//     std::ofstream outputFile(filename);

//     if (!outputFile) {
//         std::cerr << "Error opening file." << std::endl;
//         return;
//     }

//     for (int i = 0; i < headers.size(); i++) {
//         outputFile << headers[i];
//         if (i < headers.size() - 1) {
//             outputFile << ",";
//         }
//     }
//     outputFile << "\n";

//     // Write the vectors to the CSV file
//     for (size_t i = 0; i < time_total.size(); i++) {
//         if (i < time_total.size()) {
//             outputFile << time_total[i];
//         }
//         outputFile << ",";

//         if (i < time_query.size()) {
//             outputFile << time_query[i];
//         }
//         outputFile << ",";

//         if (i < time_search.size()) {
//             outputFile << time_search[i];
//         }
//         outputFile << "\n";
//     }

//     outputFile.close();
// }