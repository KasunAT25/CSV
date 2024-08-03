
#define KEY_TYPE uint64_t
#define PAYLOAD_TYPE double

#include <iostream>
#include <vector>
#include <algorithm>
#include "src/helpers/io_handler_real.h"
// g++ data_prep.cpp -std=c++17 -o data_prep
// ./data_prep
// ./data_prep ../../../mnt2/Data/genome 200000000

// ./data_prep ../../../ext/Data/fb_200M_uint64 1000000
// nohup ./data_prep normal 200000000 &
// ps aux | head -n 1 && ps aux | grep ./data_prep

// Function to read CSV file and save its contents into a vector of uint64_t
std::vector<uint64_t> readCSV(const std::string& filename);
void csv_to_bin(std::string data_intput, std::string data_output);
void remove_dup_data(std::string data_output,int argc, char *argv[]);
void split_bulk_insert(std::string data_input, double insert_prop, int inserts);
void first_x(std::string data_input,std::string data_out,  int n);

int main(int argc, char *argv[]) {

    // std::string data_output = "../../../mnt2/Data/genome.bin";
    // remove_dup_data(data_output,argc, argv);

    // std::string data_input = "genome";
    // split_bulk_insert(data_input, 0.5, 5);
    std::string data_input = "genome";
    std::string data_out = "test_genome";

    first_x(data_input,data_out, 1000000);

    // std::vector<int> bulk_load_indexes;
    

    // std::string folder = "../../../mnt2/Data/";
    // std::string bulk_index_output = folder+ "Splits/"+ data_input +"_bulk.bin";

    // bulk_load_indexes = readIndexesFromFile(bulk_index_output);

    // // for(int i =0; i <bulk_load_indexes.size(); i++ ){
    // //     std::cout << bulk_load_indexes[i] << " ";
    // // }

    // std::cout << "=========== "<< std::endl;
    // for(int j = 0; j < 5; j++){
    //     std::string insert_index_output = folder+ "Splits/"+ data_input +"_insert"+"_"+std::to_string(j)+".bin";
    //     std::vector<int> insert_indexes = readIndexesFromFile(insert_index_output);
    //     for(int i =0; i <insert_indexes.size(); i++ ){
    //         std::cout << insert_indexes[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    
    //insert_indexes = readIndexesFromFile(insert_index_output);


    //CSV to data bin
    //================
    // std::string data_intput = "../../../ext/Data/order_price_full.csv";
    // std::string data_output = "../../../ext/Data/skew.bin";
    // std::vector<KEY_TYPE> legitimate_data = readCSV(data_intput);
    // std::cout<< "Data size before" << legitimate_data.size() << std::endl;
    // std::sort(legitimate_data.begin(), legitimate_data.end());

    // // for(int i = 0; i < 100 ; i++){
    // //     std::cout<<  legitimate_data[i] << std::endl;
    // // }

    // std::cout<< "Data sorted" << std::endl;
    // legitimate_data.erase(std::unique(legitimate_data.begin(), legitimate_data.end()), legitimate_data.end());
    // save_data_bin(legitimate_data,data_output);
    // std::cout<< "Data size after" << legitimate_data.size() << std::endl;

    //Remove dup data
    //================
    // std::string data_output = "../../../ext/Data/skew.bin";

    // std::vector<KEY_TYPE> legitimate_data = parse_arguments(argc, argv);
    // std::cout<< "Data read" << std::endl;
    // std::sort(legitimate_data.begin(), legitimate_data.end());

    // std::cout<< "Data sorted" << std::endl;
    // legitimate_data.erase(std::unique(legitimate_data.begin(), legitimate_data.end()), legitimate_data.end());
    // save_data_bin(legitimate_data,data_output);

    // std::cout<< "Data saved to " << data_output << std::endl;

    

    return 0;
}

std::vector<uint64_t> readCSV(const std::string& filename) {
    std::vector<uint64_t> data; // Vector to store CSV data as uint64_t

    std::ifstream file(filename); // Open CSV file
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return data; // Return empty vector if file cannot be opened
    }

    // Skip the first line (header)
    std::string header;
    std::getline(file, header);

    std::string line;
    while (std::getline(file, line)) { // Read each line from the file
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) { // Split line into cells using comma as delimiter
            // Convert string to uint64_t and add to vector
            uint64_t value;
            try {
                value = std::stoull(cell);
            } catch (...) {
                std::cerr << "Error: Unable to convert \"" << cell << "\" to uint64_t" << std::endl;
                continue; // Skip invalid values
            }
            data.push_back(value);
        }
    }

    file.close(); // Close file
    return data;
}

void csv_to_bin(std::string data_intput, std::string data_output) {
    //std::string data_intput = "../../../ext/Data/order_price_full.csv";
    //std::string data_output = "../../../ext/Data/skew.bin";
    std::vector<KEY_TYPE> legitimate_data = readCSV(data_intput);
    std::cout<< "Data size before" << legitimate_data.size() << std::endl;
    std::sort(legitimate_data.begin(), legitimate_data.end());

    // for(int i = 0; i < 100 ; i++){
    //     std::cout<<  legitimate_data[i] << std::endl;
    // }

    std::cout<< "Data sorted" << std::endl;
    legitimate_data.erase(std::unique(legitimate_data.begin(), legitimate_data.end()), legitimate_data.end());
    save_data_bin(legitimate_data,data_output);
    std::cout<< "Data size after" << legitimate_data.size() << std::endl;
}

void remove_dup_data(std::string data_output,int argc, char *argv[]) {
   //std::string data_output = "../../../ext/Data/skew.bin";

    std::vector<KEY_TYPE> legitimate_data = parse_arguments(argc, argv);
    std::cout<< "Data read" << std::endl;
    std::sort(legitimate_data.begin(), legitimate_data.end());

    std::cout<< "Data sorted" << std::endl;
    legitimate_data.erase(std::unique(legitimate_data.begin(), legitimate_data.end()), legitimate_data.end());
    save_data_bin(legitimate_data,data_output);

    std::cout<< "Data saved to " << data_output << std::endl;
}

void split_bulk_insert(std::string data_input, double insert_prop, int inserts) {
    std::string folder = "../../../mnt2/Data/";
    std::string data_input2 = folder+data_input +".bin";
    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_input2);
    int values_size = legitimate_data.size();
    std::vector<int> all_indexes(values_size);
    std::vector<int> bulk_load_indexes;
    std::vector<int> insert_indexes;

    std::iota(all_indexes.begin(), all_indexes.end(), 0);

    unsigned seed = 123;

    std::random_device rd;
    std::mt19937 gen(seed);

    std::shuffle(all_indexes.begin(), all_indexes.end(), gen);
    std::cout << "Shuffuled." << values_size << std::endl;

    int numberOfIndexes = values_size*(1-insert_prop);
    bulk_load_indexes.resize(numberOfIndexes);

    std::vector<std::vector<int>> parts;
    insert_indexes.resize(values_size - numberOfIndexes);
    
    //Saving part
    bulk_load_indexes.assign(all_indexes.begin(), all_indexes.begin() + numberOfIndexes);
    std::sort(bulk_load_indexes.begin(), bulk_load_indexes.end());
    std::cout << "Shuffuled2." << values_size << std::endl;
    insert_indexes.assign(all_indexes.begin() + numberOfIndexes, all_indexes.end());
    //std::sort(insert_indexes.begin(), insert_indexes.end());

    size_t partSize = insert_indexes.size() / inserts;
    size_t remainder = insert_indexes.size() % inserts;

    std::string bulk_index_output = folder+ "Splits/"+ data_input +"_bulk.bin";

    saveIndexesToFile(bulk_load_indexes, bulk_index_output);
    

    auto it = insert_indexes.begin();

    for (int i = 0; i < inserts; ++i) {
        size_t currentPartSize = partSize + (remainder > 0 ? 1 : 0);
        parts.push_back(std::vector<int>(it, it + currentPartSize));
        it += currentPartSize;
        if (remainder > 0) {
            --remainder;
        }
    }
    for (int i = 0; i < inserts; ++i) {
        std::vector<int> part = parts[i];
        std::string insert_index_output = folder+ "Splits/"+ data_input +"_insert"+"_"+std::to_string(i)+".bin";
        saveIndexesToFile(part, insert_index_output);
    }
    

}

void first_x(std::string data_input,std::string data_out,  int n) {
    std::string folder = "../../../mnt2/Data/";
    std::string data_input2 = folder+data_input +".bin";
    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_input2);

    std::vector<KEY_TYPE> copiedVector;
    copiedVector.reserve(n); // Reserve space for n elements (optional but can improve performance)

    std::copy_n(legitimate_data.begin(), n, std::back_inserter(copiedVector));

    std::string data_output = folder+data_out+".bin";
    save_data_bin(copiedVector,data_output);


}