#pragma once

#include <algorithm>
#include <chrono>
#include <iostream>
#include <functional>
#include <fstream>
#include <thread>
#include <vector>
#include <cassert>
#include <cmath>
#include <math.h>
#include <cstring>
#include <climits>
#include <random>
#include <string>


#include <algorithm>
#include <sstream>

#include "competitors/ALEX/src/benchmark/zipf.h"
//#define KEY_TYPE_io uint64_t



void save_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename);
void save_to_csv2(const std::vector<KEY_TYPE>& keys,const std::vector<double>& time_total,const std::vector<std::string>& headers, const std::string& filename);
void save_segs_to_csv(const std::vector<std::vector<uint64_t>>& data, const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename);
void save_level_to_csv(const std::vector<KEY_TYPE>& keys,const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<short>& level,const std::vector<std::string>& headers, const std::string& filename);

template<typename T>
static std::vector<T> load_data(char * filename_cs) {
    std::string filename(filename_cs);
    std::vector<T> data;

    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    // Read size.
    uint64_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    data.resize(size);
    // Read values.
    in.read(reinterpret_cast<char*>(data.data()), size*sizeof(T));
    in.close();
    
    return data;
}


std::vector<u_int64_t> cache_data;
/*
    Parses first 2 arguments to create data
*/
std::vector<u_int64_t> parse_arguments(int argc, char *argv[]){
    if(argc < 3){
        std::cout << "Specify Dataset and Size" << std::endl;
        std::cout << "Available option are:" << std::endl;
        std::cout << "path_to_dataset" << std::endl;
        std::cout << "rand" << std::endl;
        std::cout << "normal" << std::endl;
        std::cout << "lognormal" << std::endl;
        std::cout << "randint" << std::endl;
        std::cout << "debug" << std::endl;
        std::cout << "exp" << std::endl;
        std::cout << "outlier" << std::endl;
        std::cout << "poisoning" << std::endl;
        exit(0);
    }
    
    long num_data = strtol(argv[2], nullptr, 0);
    std::vector<u_int64_t> data;

    std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(12346);
    
    if(strcmp(argv[1],"rand") == 0){
        std::uniform_real_distribution<double> dis(0, 2147483647);
        //std::uniform_real_distribution<double> dis(0, 50000);
        // for (long i = 0; i < num_data; i++){
        //     data.push_back(dis(gen));
        // }
        while (data.size() < num_data) {
            //double randomNum = dis(gen);
            uint64_t randomNum = static_cast<uint64_t>(std::round(dis(gen)));
            if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
                data.push_back(randomNum);
            }
    }
    }
    else if(strcmp(argv[1],"normal") == 0){
        // std::normal_distribution<double> dis(6000,6000*0.5);
       
        // for (long i = 0; i < num_data; i++){
        //     data.push_back(dis(gen));
        // }
        // std::normal_distribution<double> dis(600000,600000*1.5);
        // while (data.size() < num_data) {
        //     int randomNum = static_cast<uint64_t>(std::round(std::abs(dis(gen))));
        //     if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
        //         data.push_back(randomNum);
        //     }
        // }
         std::vector<uint64_t> all_values(2147483647);
         std::cout<< "Data vector" << std::endl;
         std::iota(all_values.begin(), all_values.end(), 0);
         std::cout<< "Data created" << std::endl;

        // Shuffle the array
        std::shuffle(all_values.begin(), all_values.end(), gen);
        std::cout<< "Data shuffled" << std::endl;
        for (uint64_t i = 0; i < num_data; ++i) {
            data.push_back(all_values[i]);
         }

    }
     else if(strcmp(argv[1],"normald") == 0){
        std::normal_distribution<double> dis(6000,6000*0.5);
        // for (long i = 0; i < num_data; i++){
        //     data.push_back(dis(gen));
        // }
        while (data.size() < num_data) {
            // double randomNum = std::abs(dis(gen));
            
            // //if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
            //     data.push_back(randomNum);
            //}
            double value = std::abs(dis(gen));
            value = static_cast<double>(static_cast<long long>(value * 10000)) / 10000.0;

            if (std::find(data.begin(), data.end(), value) == data.end()) {
                data.push_back(value);
            }

    }
    }
    else if(strcmp(argv[1],"lognormal") == 0){
        std::lognormal_distribution<double> dis(0.0, 4.0);
        // for (long i = 0; i < num_data; i++){
        //     data.push_back(dis(gen));
        // }
        while (data.size() < num_data) {
            double randomNum = dis(gen);
            if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
                data.push_back(randomNum);
            }
    }
    }
    else if(strcmp(argv[1],"randint") == 0){
        // for (long i = 0; i < num_data; i++){
        //     long val = rand();
        //     data.push_back(val);
        // }
        while (data.size() < num_data) {
            double randomNum = rand();
            if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
                data.push_back(randomNum);
            }
    }

    }
    else if(strcmp(argv[1],"debug") == 0){
        for (long i = 0; i < num_data; i++){
            long val = rand()%30;
            data.push_back(val);
        }
    }
    else if(strcmp(argv[1],"randfixed") == 0){
        for (long i = 0; i < num_data; i++){
            std::uniform_real_distribution<double> dis(0, 10000);
            data.push_back(dis(gen));
        }
    }
    else if(strcmp(argv[1],"exp") == 0){
        for (long i = 0; i < num_data; i++){
            long val = std::pow(2,i);
            data.push_back(val);
        }

        
    }
    else if(strcmp(argv[1],"manual") == 0){
        for (long i = 2; i < argc; i++){
            data.push_back(strtol(argv[i], nullptr, 0));
        }
        std::cout << data.size() << std::endl;
        std::sort(data.begin(), data.end());
        return data;
    }
    else if(strcmp(argv[1],"outlier") == 0){
        int num_outliers = 1;
        for (long i = 0; i < num_data - num_outliers; i++){
            std::uniform_real_distribution<double> dis(0, 100000);
            data.push_back(dis(gen));
        }
        for (long i = 0; i < num_outliers; i++){
            std::uniform_real_distribution<double> dis(200000000000000, 300000000000000);
            data.push_back(dis(gen));
        }
    }
    else if(strcmp(argv[1],"poisoning") == 0){
            std::uniform_int_distribution<int> dis(1, 100000);
            for (long i = 0; i < num_data; i++){
                data.push_back(dis(gen));
            }
        }
    else if(strcmp(argv[1],"test") == 0){
        //     std::uniform_int_distribution<int> dis(1, num_data+25);
        //     // for (long i = 0; i < num_data; i++){
        //     //     data.push_back(dis(gen));
        //     // }

        //     while (data.size() < num_data) {
        //         int randomNum = dis(gen);
        //         if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
        //             data.push_back(randomNum);
        //         }
        // }
        data.push_back(3);
        data.push_back(5);
        data.push_back(6);
        //data.push_back(6);
        data.push_back(7);
        data.push_back(12);
        data.push_back(24);
        data.push_back(33);
        data.push_back(35);
        data.push_back(36);
        data.push_back(38);
        }
    else if(strcmp(argv[1],"test1") == 0){
        //     std::uniform_int_distribution<int> dis(1, num_data+25);
        //     // for (long i = 0; i < num_data; i++){
        //     //     data.push_back(dis(gen));
        //     // }

        //     while (data.size() < num_data) {
        //         int randomNum = dis(gen);
        //         if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
        //             data.push_back(randomNum);
        //         }
        // }
        data.push_back(608033450040);
        data.push_back(608296912800);
        data.push_back(612449859856);
        data.push_back(624354272984);
        data.push_back(629948030656);
        data.push_back(630758380512);
        data.push_back(630991398176);
        data.push_back(631319143936);
        data.push_back(631665836968);
        data.push_back(631830789600);

        data.push_back(632387432296);
        data.push_back(633041194440);
        data.push_back(633202968824);
        data.push_back(642455171888);
        data.push_back(647995867792);
        data.push_back(656571223560);
        data.push_back(656748982424);
        data.push_back(667255469032);
        
        }
    else if(strcmp(argv[1],"rest") == 0){
            std::uniform_int_distribution<int> dis(1, num_data+num_data*0.5);
            // for (long i = 0; i < num_data; i++){
            //     data.push_back(dis(gen));
            // }

            while (data.size() < num_data) {
                int randomNum = dis(gen);
                if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
                    data.push_back(randomNum);
                }
        }
        }
    // else if(){

    // }
    else{
        if (cache_data.size() == 0 || true){
            std::cout << argv[1] << "\n";
            cache_data = load_data<u_int64_t>(argv[1]);
        }

        data.reserve(num_data);
        std::uniform_int_distribution<long> dis(0, std::max<long>((long)cache_data.size()-1-num_data, 0));


       
        // int start = dis(gen);
        //int start = 1000000;
        
        //int start = 5000000;
        //int start = 50000000;
        
        //I used 1000 for all testing
        int start = 0;
        // int start = 1000;
        //int start = 100000000;
        //int start = 98000000;
        for(long i = 0; cache_data[i] < LONG_MAX && i < num_data && i < cache_data.size(); i++){
            data.push_back((double)cache_data[start+i]);
        }

        //To get randomly
        // std::uniform_int_distribution<size_t> rand_dis(0, cache_data.size() - 1);
        // //unsigned seed = 123;
        
        // unsigned seed = 987;
        // //std::random_device rd;
        // std::mt19937 generator(seed);
        // std::vector<int> data_test;

        // //get the random indexes to search
        // for(int i =0; i < num_data; i++){
        //     int idx = rand_dis(generator);
        //     data_test.push_back(idx);
        // }
        // for(long i = 0; cache_data[i] < LONG_MAX && i < num_data && i < cache_data.size(); i++){
        //     data.push_back((double)cache_data[data_test[i]]);
        // }
    }
    std::sort(data.begin(), data.end());
    //data.erase( unique( data.begin(), data.end() ), data.end() );
    //std::sort(data.begin(), data.end());
    return data;
}
void save_to_csv(const std::vector<KEY_TYPE>& keys,const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename) {
    
    std::ofstream outputFile(filename);

    if (!outputFile) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    for (int i = 0; i < headers.size(); i++) {
        outputFile << headers[i];
        if (i < headers.size() - 1) {
            outputFile << ",";
        }
    }
    outputFile << "\n";

    // Write the vectors to the CSV file
    for (size_t i = 0; i < time_total.size(); i++) {
        if (i < keys.size()) {
            outputFile << keys[i];
        }
        outputFile << ",";

        if (i < time_total.size()) {
            outputFile << time_total[i];
        }
        outputFile << ",";

        if (i < time_query.size()) {
            outputFile << time_query[i];
        }
        outputFile << ",";

        if (i < time_search.size()) {
            outputFile << time_search[i];
        }
        outputFile << "\n";
    }

    outputFile.close();
}

void save_to_csv2(const std::vector<KEY_TYPE>& keys,const std::vector<double>& time_total,const std::vector<std::string>& headers, const std::string& filename) {
    
    std::ofstream outputFile(filename);

    if (!outputFile) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    for (int i = 0; i < headers.size(); i++) {
        outputFile << headers[i];
        if (i < headers.size() - 1) {
            outputFile << ",";
        }
    }
    outputFile << "\n";

    // Write the vectors to the CSV file
    for (size_t i = 0; i < time_total.size(); i++) {
        if (i < keys.size()) {
            outputFile << keys[i];
        }
        outputFile << ",";

        if (i < time_total.size()) {
            outputFile << time_total[i];
        }
        outputFile << "\n";
    }

    outputFile.close();
}

void save_level_to_csv(const std::vector<KEY_TYPE>& keys,const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,
 const std::vector<short>& level,const std::vector<std::string>& headers, const std::string& filename) {
    
    std::ofstream outputFile(filename);

    if (!outputFile) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    for (int i = 0; i < headers.size(); i++) {
        outputFile << headers[i];
        if (i < headers.size() - 1) {
            outputFile << ",";
        }
    }
    outputFile << "\n";

    // Write the vectors to the CSV file
    for (size_t i = 0; i < time_total.size(); i++) {
        if (i < keys.size()) {
            outputFile << keys[i];
        }
        outputFile << ",";

        if (i < time_total.size()) {
            outputFile << time_total[i];
        }
        outputFile << ",";

        if (i < time_query.size()) {
            outputFile << time_query[i];
        }
        outputFile << ",";

        if (i < time_search.size()) {
            outputFile << time_search[i];
        }
         outputFile << ",";

        // if (i < diffs.size()) {
        //     outputFile << diffs[i];
        // }
        //  outputFile << ",";

        if (i < level.size()) {
            outputFile << level[i];
        }
       
        outputFile << "\n";
    }

    outputFile.close();
}

void save_segs_to_csv(const std::vector<std::vector<uint64_t>>& data, const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename) {
    
    std::ofstream outputFile(filename);

    if (!outputFile) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    for (int i = 0; i < headers.size(); i++) {
        outputFile << headers[i];
        if (i < headers.size() - 1) {
            outputFile << ",";
        }
    }
    outputFile << "\n";

    // Write the vectors to the CSV file
    int index =0;
    for(int j=0; j< data.size();j++){

        
        for (size_t i = 0; i < data[j].size(); i++) {
            if (i < data[j].size()) {
                outputFile << data[j][i];
            }
            outputFile << ",";

            outputFile << j;
            outputFile << ",";
            if (index < time_total.size()) {
                outputFile << time_total[index];
            }
            outputFile << ",";

            if (index < time_query.size()) {
                outputFile << time_query[index];
            }
            outputFile << ",";

            if (index < time_search.size()) {
                outputFile << time_search[index];
            }
            outputFile << "\n";
            index += 1;
        }
    }

    outputFile.close();
}


void save_data(std::vector<uint64_t>& vector, const std::string& filename){
    std::ofstream outputFile(filename); // Open a file for writing

    if (outputFile.is_open()) {
        for (KEY_TYPE element : vector) {
            outputFile << element << "\n"; // Write each element to the file with a newline
        }
        outputFile.close(); // Close the file
        std::cout << "Vector elements saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file for writing!" << std::endl;
    }
}

void saveIndexesToFile(const std::vector<int>& indexes, const std::string& filename) {

    // std::ofstream outputFile(filename, std::ios::binary);
    // if (outputFile.is_open()) {
    //     for (int index : indexes) {
    //         outputFile << index << std::endl;
    //     }
    //     outputFile.close();
    //     std::cout << "Indexes saved to " << filename << std::endl;
    // } else {
    //     std::cerr << "Unable to open file: " << filename << std::endl;
    // }
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
    std::cout << "Indexes saved to " << filename << std::endl;
    size_t size = indexes.size();
    outfile.write(reinterpret_cast<const char*>(&size), sizeof(size_t)); // Write the size of the vector

    for (const auto& elem : indexes) {
        outfile.write(reinterpret_cast<const char*>(&elem), sizeof(int)); // Write each element
    }

    outfile.close();
}

// Function to read indexes from a text file into a vector
std::vector<int> readIndexesFromFile(const std::string& filename) {
    // std::vector<int> indexes;
    // std::ifstream inputFile(filename, std::ios::binary);
    // if (inputFile.is_open()) {
    //     int index;
    //     while (inputFile >> index) {
    //         indexes.push_back(index);
    //     }
    //     inputFile.close();
    //     std::cout << "Indexes read from " << filename << std::endl;
    // } else {
    //     std::cerr << "Unable to open file: " << filename << std::endl;
    // }
    // return indexes;
    std::ifstream infile(filename, std::ios::binary);
    std::vector<int> vec;
    std::cout << "Indexes read from " << filename << std::endl;
    if (!infile) {
        std::cerr << "Error: Unable to open file for reading." << std::endl;
        return vec;
    }

    size_t size;
    infile.read(reinterpret_cast<char*>(&size), sizeof(size_t)); // Read the size of the vector

    vec.resize(size);
    for (auto& elem : vec) {
        infile.read(reinterpret_cast<char*>(&elem), sizeof(int)); // Read each element
    }

    infile.close();
    return vec;
}

void save_data_bin(const std::vector<KEY_TYPE>& indexes, const std::string& filename) {

    // std::ofstream outputFile(filename, std::ios::binary);
    // if (outputFile.is_open()) {
    //     for (int index : indexes) {
    //         outputFile << index << std::endl;
    //     }
    //     outputFile.close();
    //     std::cout << "Indexes saved to " << filename << std::endl;
    // } else {
    //     std::cerr << "Unable to open file: " << filename << std::endl;
    // }
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
    //std::cout << "dar read from " << filename << std::endl;
    size_t size = indexes.size();
    outfile.write(reinterpret_cast<const char*>(&size), sizeof(size_t)); // Write the size of the vector

    for (const auto& elem : indexes) {
        outfile.write(reinterpret_cast<const char*>(&elem), sizeof(KEY_TYPE)); // Write each element
    }

    outfile.close();
}

// Function to read indexes from a text file into a vector
std::vector<KEY_TYPE> read_data_bin(const std::string& filename) {
    // std::vector<int> indexes;
    // std::ifstream inputFile(filename, std::ios::binary);
    // if (inputFile.is_open()) {
    //     int index;
    //     while (inputFile >> index) {
    //         indexes.push_back(index);
    //     }
    //     inputFile.close();
    //     std::cout << "Indexes read from " << filename << std::endl;
    // } else {
    //     std::cerr << "Unable to open file: " << filename << std::endl;
    // }
    // return indexes;
    std::ifstream infile(filename, std::ios::binary);
    std::vector<KEY_TYPE> vec;
    
    if (!infile) {
        std::cerr << "Error: Unable to open file for reading." << std::endl;
        return vec;
    }

    size_t size;
    infile.read(reinterpret_cast<char*>(&size), sizeof(size_t)); // Read the size of the vector

    vec.resize(size);
    for (auto& elem : vec) {
        infile.read(reinterpret_cast<char*>(&elem), sizeof(KEY_TYPE)); // Read each element
    }

    infile.close();
    std::cout << "data read from " << filename << std::endl;

    return vec;
}

std::vector<KEY_TYPE> get_search_keys(std::vector<KEY_TYPE> array, int num_keys, int num_searches, unsigned seed) {
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<int> dis(0, num_keys - 1);
  std::vector<KEY_TYPE> keys(num_searches);
  for (int i = 0; i < num_searches; i++) {
    int pos = dis(gen);
    keys[i] = array[pos];
  }
  return keys;
}

std::vector<KEY_TYPE> get_search_keys_zipf(std::vector<KEY_TYPE> array, int num_keys, int num_searches, unsigned seed) {
 // auto* keys = new KEY_TYPE[num_searches];
  std::vector<KEY_TYPE> keys(num_searches);
  ScrambledZipfianGenerator zipf_gen(num_keys,seed);
  for (int i = 0; i < num_searches; i++) {
    int pos = zipf_gen.nextValue();
    keys[i] = array[pos];
  }
  return keys;
}

void saveToCSV(const std::string& data, const std::string& filename) {
   std::ofstream outFile(filename, std::ios::app); // Open file in append mode

    if (!outFile) {
        std::cerr << "Error: Unable to open or create file " << filename << std::endl;
        return;
    }

    // Parse the input data string
    std::stringstream ss(data);
    std::string token;
    while (std::getline(ss, token, ';')) {
        outFile << token << ",";
    }
    outFile << std::endl; // End the line
    outFile.close();
    std::cout << "Data appended to " << filename << " successfully." << std::endl;
}