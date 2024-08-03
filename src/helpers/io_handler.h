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

void save_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename);
void save_level_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<short>& level,const std::vector<std::string>& headers, const std::string& filename);

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
std::vector<double> parse_arguments(int argc, char *argv[]){
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
    std::vector<double> data;
    //std::vector<uint64_t> data;

    std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(123457);
    
    if(strcmp(argv[1],"rand") == 0){
        // std::uniform_real_distribution<double> dis(0, 2147483647);
        std::uniform_real_distribution<double> dis(0, 50000);
        // for (long i = 0; i < num_data; i++){
        //     data.push_back(dis(gen));
        // }
        while (data.size() < num_data) {
            //double randomNum = dis(gen);
            int randomNum = static_cast<int>(std::round(dis(gen)));
            if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
                data.push_back(randomNum);
            }
    }
    }
    else if(strcmp(argv[1],"normal") == 0){
        std::normal_distribution<double> dis(6000,6000*0.5);
        // for (long i = 0; i < num_data; i++){
        //     data.push_back(dis(gen));
        // }
        while (data.size() < num_data) {
            int randomNum = static_cast<int>(std::round(dis(gen)));
            if (std::find(data.begin(), data.end(), randomNum) == data.end()) {
                data.push_back(randomNum);
            }
    }
    }
     else if(strcmp(argv[1],"normald") == 0){
        std::normal_distribution<double> dis(50000,50000*0.5);
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
        data.push_back(7);
        data.push_back(12);
        data.push_back(24);
        data.push_back(33);
        data.push_back(35);
        data.push_back(36);
        data.push_back(38);
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
    else if(strcmp(argv[1],"read") == 0){
            std::cout << "Reading File " << std::endl;

    std::ifstream file("/home/ubuntu/codes/poisoning/LogarithmicErrorRegression-master/poisoned_data2.csv");

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open the file." << std::endl;
        //return 1;
    }

    std::vector<double> keyset2;

    std::string header;
    std::getline(file, header);
   

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream linestream(line);
        std::string cell;
        
        // Split the line using commas as delimiters
        std::getline(linestream, cell, ',');
        double value1 = std::stod(cell);
        data.push_back((value1));

    }

    // Close the file
    file.close();
    std::cout << "File Read" << std::endl;
        }
    else{
        if (cache_data.size() == 0){
            std::cout << argv[1] << "\n";
            cache_data = load_data<uint64_t>(argv[1]);
        }
        

        data.reserve(num_data);
        std::uniform_int_distribution<long> dis(0, std::max<long>((long)cache_data.size()-1-num_data, 0));
        int start = dis(gen);
        for(long i = 0; cache_data[i] < LONG_MAX && i < num_data && i < cache_data.size(); i++){
            data.push_back((double)cache_data[start+i]);
        }
    }
    std::sort(data.begin(), data.end());
    //data.erase( unique( data.begin(), data.end() ), data.end() );
    //std::sort(data.begin(), data.end());
    return data;
}
void save_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename) {
    
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

void save_level_to_csv(const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<short>& level,const std::vector<std::string>& headers, const std::string& filename) {
    
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
        if (i < level.size()) {
            outputFile << level[i];
        }
       
        outputFile << "\n";
    }

    outputFile.close();
}

void save_segs_to_csv(const std::vector<std::vector<double>>& data, const std::vector<double>& time_total,const std::vector<double>& time_query,const std::vector<double>& time_search,const std::vector<std::string>& headers, const std::string& filename) {
    
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

void save_data(std::vector<double>& vector, const std::string& filename){
    std::ofstream outputFile(filename); // Open a file for writing

    if (outputFile.is_open()) {
        for (const auto& element : vector) {
            outputFile << element << "\n"; // Write each element to the file with a newline
        }
        outputFile.close(); // Close the file
        std::cout << "Vector elements saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file for writing!" << std::endl;
    }
}