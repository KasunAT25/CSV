
#define KEY_TYPE uint64_t
#define PAYLOAD_TYPE double


#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>

#include "src/log_regression.h"
#include "src/competitor_regression.h"
#include "src/irls.h"


#include "src/helpers/io_handler_real.h"
#include "src/smooth_simple_v2.h"


#include "src/helpers/alex_benchmark.h"
#include "src/helpers/sali_benchmark.h"

#include "src/fast_brute_force_real.h"

#include "src/theil_sen.h"

#include <unordered_set>
#include <unistd.h>


#include <thread>
#include <mutex>
#include <unordered_map>

typedef sali::SALI<KEY_TYPE, PAYLOAD_TYPE, true>::Node Node;

bool compare_level_des(const Node* a, const Node* b);
void scan_node(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload);
void scan_node_all(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload);

void get_nodes_by_level( std::vector<Node*> *nodes, std::vector<Node*> *nodes_by_level, int level);
void get_poisoned_keys(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi, std::vector<PAYLOAD_TYPE> payload, KEY_TYPE* keys, PAYLOAD_TYPE* payloads, int size);
void get_nodes_by_level_with_child( std::vector<Node*> *nodes, std::vector<Node*> *nodes_by_level, int level);
void shift_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
void shift_back_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
void get_nodes_children(Node* node, std::vector<Node*> *nodes_by_level);
void count_data_node(Node* node, int * count);

void get_data_node_max_height(Node* node, int * height);
void get_data_node_data_by_level(Node* node, std::vector<std::vector<KEY_TYPE>*> & node_data);
void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1,  std::vector<KEY_TYPE> vec2);

std::pair<uint64_t, PAYLOAD_TYPE> * create_values(std::vector<KEY_TYPE> data,int * size);

//To create the values from smooth data
std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_poisoned_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, int size);


bool check_duplicates(std::vector<KEY_TYPE> data);

int csv_node(std::vector<Node*> nodes_by_level, KEY_TYPE P_i, size_t start, size_t end);
void get_node_data_by_level(Node* node, std::set<KEY_TYPE>& node_data, int level);
void get_node_data_below_level(Node* node, std::set<KEY_TYPE>& node_data, int level);


void clearMemoryCache() {
      if (system("sync && sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'") != 0) {
        std::cerr << "Failed to clear memory cache." << std::endl;
    }
}

std::mutex resultsMutex;
std::unordered_map<int, std::tuple<
                    KEY_TYPE, 
                    PAYLOAD_TYPE,
                    KEY_TYPE*,
                    PAYLOAD_TYPE*,
                    int,
                    KEY_TYPE*,
                    int>> values_for_inner;

sali::SALI<KEY_TYPE, PAYLOAD_TYPE> index2;

int main(int argc, char *argv[]){
    // clearMemoryCache();
    // flushCache();

    std::string dataset_name = argv[1];
    std::string poi_size_str = argv[4];
    std::string insert_prop_str = argv[5];
    
    bool benchmark = 1;
    //std::string data_folder = "data/";
    std::string data_folder = "/data/";

    std::string data_output = data_folder + dataset_name+".bin";

    bool smooth = (argv[3][0] == '1');
    
    double insert_threshold = std::stod(insert_prop_str);
    bool insert = (insert_threshold > 0.0);

    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_output);
    int values_size = legitimate_data.size();

    KEY_TYPE orignal_P = static_cast<KEY_TYPE>(std::stod(poi_size_str)*values_size);

    unsigned seed = 123;

    KEY_TYPE P = orignal_P;
    KEY_TYPE P_left = P;

    double constant_cost = 1;
    int pow_val = 1;

    double original_poi_thres = 0;
    std::string method = "nt6";

    if(smooth){
        original_poi_thres = 0.5;
        method = "ntp6";
    }
    
    double poi_thres = original_poi_thres;
    long pois_count = poi_thres*legitimate_data.size();
    bool use_new_cost = true;
    int num_total_smooth = 0;
    int max_check = 10;
    int insert_batches = 5;
    
    
    srand(12345);

    //Get the dataname, size and the method to save them
    std::string data_name = argv[1];
    std::string data_name2 = "osm_200M_uint64";
    
    std::string data_size = argv[2];

    //filenames of the benchmark
    std::string model_output = "results/model_output_real";

    std::string output_folder = "";
    std::string performance_output = output_folder+"results/sali/rev_sali_original_performance.csv";
    std::string structure_output = output_folder+"results/sali/rev_sali_original_structure.csv";

    if(smooth){
        performance_output = output_folder+"results/sali/rev_sali_smooth_performance.csv";
        structure_output = output_folder+"results/sali/rev_sali_smooth_structure.csv";
    }

    std::string changed_data_output = "results/changed_data.bin";
    std::string output_model_csv = output_folder+"results/sali/rev_sali_model_output_good.csv";
    

    //To save the changed data
    std::vector<KEY_TYPE> changed_data;

    //Create the values (append new payloads)
    std::mt19937_64 gen_payload(std::random_device{}());
    values_size = legitimate_data.size();

    std::cout << "dataset without duplicates " << values_size << std::endl;
    
    

    //for inserts
    std::vector<int> all_indexes(values_size);
    std::vector<int> bulk_load_indexes;
    std::vector<int> insert_indexes;

    std::iota(all_indexes.begin(), all_indexes.end(), 0);

    auto start_build = std::chrono::high_resolution_clock::now();

    //Depending on the insertion 
    int numberOfIndexes = 0;
    if(insert){
        std::string bulk_index_output = data_folder + "splits/"+  dataset_name +"_bulk.bin";
        bulk_load_indexes = readIndexesFromFile(bulk_index_output);

        std::vector<KEY_TYPE> subvector(bulk_load_indexes.size());

        subvector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(subvector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });
        std::cout << "bulk_load_indexes "  << bulk_load_indexes.size() << std::endl;
        std::cout << "insert_indexes "  << insert_indexes.size() << std::endl;

        std::cout << "subvector "  << subvector.size() << std::endl;
        numberOfIndexes = subvector.size();
        auto values = create_values(subvector,&numberOfIndexes);
        std::cout << "Created values" << std::endl;
        
        index2.bulk_load(values, numberOfIndexes);
        

    }
    else{
        std::cout << "Creating values" << std::endl;
        auto values = create_values(legitimate_data,&values_size);

        std::cout << "Created values" << std::endl;
       
        index2.bulk_load(values, values_size);
        delete[] values; // Delete the dynamically allocated array
        values = nullptr; // Reset the pointer to nullptr
    }
    auto stop_build = std::chrono::high_resolution_clock::now();

    std::vector<int>().swap(all_indexes);
    // std::vector<int>().swap(bulk_load_indexes);

    std::cout << "============================" << std::endl;

    std::cout << "Created the index" << std::endl;
    std::cout << "============================" << std::endl;

    long index_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_build - start_build).count();
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;

    
    auto start_tra1 = std::chrono::high_resolution_clock::now();

    int max_model_height =0;

    std::vector<Node*> nodes;
   
    index2.scan_nodes(&nodes);

    std::sort(nodes.begin(), nodes.end(), compare_level_des);
    max_model_height = nodes[0]->level;
    
    std::cout << "DONE: obtaining models useful and models by level" << std::endl;
    std::cout << "============================" << std::endl;


    int current_level = max_model_height;

    //GET USEFUL NODES BY LEVEL
    //============================
    std::vector<Node*> nodes_by_level;
    current_level = 2;
    get_nodes_by_level_with_child( &nodes, &nodes_by_level, current_level);
    
    std::cout << "Current_level " << current_level << std::endl;
    std::cout << "nodes_by_level " << nodes_by_level.size() << std::endl;

    //Do for all the levels except the root (keep the root as it is)
    //Inner index to iterate through the models in the current level (reset in the outter while loop)
    int inner_idx = 0;
    int current_count = 0;

    std::set<KEY_TYPE>  altered_data_set;
    std::set<KEY_TYPE>  demoted_data_set;

    int size_to_check;

    if(nodes_by_level.size() <= 0){
        current_level = 0;
    }
    
    auto stop_tra1 = std::chrono::high_resolution_clock::now();
    long tra_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra1 - start_tra1).count();

    
    while(current_level == 2 && smooth){
        auto start_tra2 = std::chrono::high_resolution_clock::now();

        Node* best_node;
        Node* new_node;
        bool model_converted = false;
        inner_idx = 0;
    
        size_to_check = nodes_by_level.size();
        int original_size_to_check = size_to_check;

        KEY_TYPE P_i = orignal_P/size_to_check;

        if(P_i <= 0){
            P_i = 1;
        }

        int numThreads = 64; // Set number of threads
        size_t chunkSize = (nodes_by_level.size() + numThreads - 1) / numThreads;
        std::vector<std::thread> threads;
        int total_poi = 0;

        //Finding the virtual points using threads
        for(int num_thread = 0;num_thread< numThreads;num_thread++){
            size_t start = num_thread * chunkSize;
            size_t end = std::min(start + chunkSize, nodes_by_level.size());

            if (start < end) {
                threads.emplace_back(csv_node, std::ref(nodes_by_level),P_i, start, end);
            }
        }

        for (auto& thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        auto stop_tra2 = std::chrono::high_resolution_clock::now();
        tra_time += std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra2 - start_tra2).count();
        //Inserting the virtual points into the nodes and resturcturing
        for(int inner_idx2 = 0; inner_idx2< nodes_by_level.size(); inner_idx2++ ){
            auto start_tra3 = std::chrono::high_resolution_clock::now();

            Node* best_node = nodes_by_level[inner_idx2];
            auto&[level_data1, payload1,keys_poi,values_poi,poisoned_model_data_size,keys,size] = values_for_inner[inner_idx2];

            auto stop_tra3 = std::chrono::high_resolution_clock::now();
            tra_time += std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra3 - start_tra3).count();
            
            if(poisoned_model_data_size > size){
                num_total_smooth = num_total_smooth + poisoned_model_data_size - size;
                //Getting original data
                std::set<KEY_TYPE> node_data_original;
                std::set<KEY_TYPE> node_data_original_level;

                int level_now = best_node->level;

                get_node_data_below_level(best_node,node_data_original, level_now);
                
                get_node_data_by_level(best_node,node_data_original_level, level_now);

                auto start_tra4 = std::chrono::high_resolution_clock::now();

                Node* new_node ;

                new_node = index2.poison_bulk_load_real(index2.root, best_node,level_data1, payload1, keys_poi, values_poi,poisoned_model_data_size,keys, size);

                auto stop_tra4 = std::chrono::high_resolution_clock::now();
                tra_time += std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra4 - start_tra4).count();

                std::set<KEY_TYPE> node_data_changed;
                std::set<KEY_TYPE> node_data_changed_below;
                
                get_node_data_by_level(new_node, node_data_changed, level_now);

                get_node_data_below_level(new_node,node_data_changed_below, level_now);

                std::set<KEY_TYPE> tempIntersection;
                std::set_intersection(node_data_original.begin(), node_data_original.end(),
                                        node_data_changed.begin(), node_data_changed.end(),
                                        std::inserter(tempIntersection, tempIntersection.begin()));

                altered_data_set.insert(tempIntersection.begin(), tempIntersection.end());


                std::set<KEY_TYPE> tempIntersection2;
                std::set_difference(node_data_original_level.begin(), node_data_original_level.end(),
                                        node_data_changed.begin(), node_data_changed.end(),
                                        std::inserter(tempIntersection2, tempIntersection2.begin()));

              
                demoted_data_set.insert(tempIntersection2.begin(), tempIntersection2.end());

                

            }
        }
        
        current_level--;

        //Not needed for LIPP and SALI
        if(current_level  > 1){
        }
        
        std::cout << "=============="  <<std::endl;
        std::cout << "Done: Iteration " << current_level <<std::endl;
            
    }
    

    int max_model_height_fin = 0;
    
    //GET NEW NODE STATS 
    //=====================

    //auto stop_tra = std::chrono::high_resolution_clock::now();

    //The total size of the index after virtual points
    std::cout << " Number " << index2.root->size <<std::endl;
    std::cout << "Virtual points number " << num_total_smooth <<std::endl;
    std::cout << "Number of possible nodes " << size_to_check << std::endl;

    //long tra_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra - start_tra).count();

    //Loading the altered data or saving it
    std::vector<KEY_TYPE> lookup_changed;
    std::vector<KEY_TYPE>  altered_data(altered_data_set.begin(), altered_data_set.end());
    std::vector<KEY_TYPE>  demoted_data(demoted_data_set.begin(), demoted_data_set.end());;
    int altered_data_size = altered_data.size();
    int demoted_data_size = demoted_data.size();
    std::cout << "Size promoted " <<  altered_data_size << std::endl;
    std::cout << "Size demoted " <<  demoted_data_size << std::endl;

    if(smooth){
        lookup_changed = get_search_keys(altered_data, altered_data_size, 1000000, seed);
        save_data_bin(lookup_changed,changed_data_output);
    }
    if(!smooth&benchmark){
        lookup_changed.clear();
        lookup_changed = read_data_bin(changed_data_output);
        
    }

//BENCHMARKING FOR INSERTS & READ-ONLY
//========================

//To hold the lookup keys

    std::vector<KEY_TYPE> lookup_keys;
    std::vector<KEY_TYPE> lookup_keys_zipf;

    //If inserts then select only from the original bulkloaded keys
    if(insert){

        std::vector<KEY_TYPE> bulk_vector(bulk_load_indexes.size());
        bulk_vector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(bulk_vector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });
        int bulk_size = bulk_load_indexes.size();

        lookup_keys = get_search_keys(bulk_vector, bulk_size, 1000000, seed);
        lookup_keys_zipf = get_search_keys_zipf(bulk_vector, bulk_size, 1000000, seed);

        //Performing the query to measure the original times
         benchmark_sali_real(index2,lookup_changed, lookup_changed, "SALI", data_name+ "_changed", poi_thres, performance_output,
         smooth,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index2,lookup_keys, lookup_keys, "SALI", data_name+"_lookup", poi_thres, performance_output,
         smooth,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index2,lookup_keys_zipf, lookup_keys_zipf, "SALI", data_name+ "_lookup_zipf", poi_thres, performance_output,
         smooth,insert,orignal_P,insert_threshold);

    }
    //If read-only then just collect the keys.
    else{
        lookup_keys = get_search_keys(legitimate_data, values_size, 1000000, seed);
        lookup_keys_zipf = get_search_keys_zipf(legitimate_data, values_size, 1000000, seed);
    }
    

    //BENCHMARKING LEVEL DETAILS

    std::vector<Node*> nodes2;
    index2.scan_nodes(&nodes2);
    std::sort(nodes2.begin(), nodes2.end(), compare_level_des);

    max_model_height = nodes2[0]->level;

    size_t total_node_count = 0;
    size_t total_data_count = 0;

    std::string structure_results = "SALI;" + data_name + ";" + std::to_string(smooth) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_smooth)+ ";" + std::to_string(altered_data_size)+ ";" + std::to_string(demoted_data_size)+ ";Data";

    std::string node_info = "Nodes";
    
    //For each of the levels of the index get the relavant data
    for(int i = 1; i <= max_model_height ; i++ ){

        std::cout << "Performance For Level  " << i << std::endl;

        std::vector<KEY_TYPE> level_data;
        std::vector<PAYLOAD_TYPE> level_payload;
        int count = 0;
        for(int j = 0; j < nodes2.size();j++){
            if(nodes2[j]->level==i){
                count++;
                scan_node(nodes2[j], &level_data, &level_payload);
            }
        }

        if(i > 2){
            total_node_count += count;
            total_data_count += level_data.size();
        }
       

        std::cout << "Nodes " << count <<std::endl;
        std::cout << "Data " << level_data.size() <<std::endl;
        structure_results = structure_results + ";" + std::to_string(level_data.size());

         node_info = node_info + ";" + std::to_string(count);

    }

    structure_results = structure_results + ";" + node_info + ";" + std::to_string(altered_data_size);

    saveToCSV(structure_results,structure_output);

    size_t index_size_after =  index2.total_size();

    std::string model_results = "SALI;" + data_name + ";" + std::to_string(smooth) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_smooth) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" +
        std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        saveToCSV(model_results,output_model_csv);


    if(insert){
        //For the 5 insert batches
        for(int j = 0; j < insert_batches; j++){
            insert_threshold = (j+1)*0.1;
            std::string structure_results2 = "SALI;" + dataset_name + "_insert_"+std::to_string(j) +";" + std::to_string(smooth) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_smooth);

            std::vector<KEY_TYPE>  inserted_data;
        
            std::string insert_index_output = data_folder+ "splits/"+ data_name +"_insert"+"_"+std::to_string(j)+".bin";
            std::vector<int> insert_indexes2 = readIndexesFromFile(insert_index_output);


            auto start = std::chrono::high_resolution_clock::now();
            for(int i : insert_indexes2){
                index2.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
                //inserted_data.push_back(legitimate_data[i]);
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count()/(1.0*insert_indexes2.size());

            std::cout << "Inserted " << std::endl;
        //  benchmark_lipp_real(index,inserted_data, inserted_data, "LIPP", dataset_name+ "_changed", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_lipp_real(index,lookup_keys, lookup_keys, "LIPP", dataset_name+"_lookup", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_lipp_real(index,lookup_keys_zipf, lookup_keys_zipf, "LIPP", dataset_name+ "_lookup_zipf", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
            std::vector<Node*> nodes3;
            index2.scan_nodes(&nodes3);
            std::sort(nodes3.begin(), nodes3.end(), compare_level_des);
            max_model_height = nodes3[0]->level;

            //index_original.scan_and_destory_tree(nodes[0], keys, values,false);
            size_t total_node_count2 = 0;
            size_t total_data_count2 = 0;
            std::string node_info2 = "Nodes";

            //For each insert batch check the performance of each level
            for(int i = 1; i <= max_model_height ; i++ ){
                
                std::cout << "Performance For Level  " << i << std::endl;

                std::vector<KEY_TYPE> level_data2;
                std::vector<PAYLOAD_TYPE> level_payload2;
                int count2 = 0;
                for(int j = 0; j < nodes3.size();j++){
                    if(nodes3[j]->level==i){
                        count2++;
                        scan_node(nodes3[j], &level_data2, &level_payload2);
                    }
                }

                if(i > 2){
                    total_node_count2 += count2;
                    total_data_count2 += level_data2.size();
                }
            

                std::cout << "Nodes " << count2 <<std::endl;
                std::cout << "Data " << level_data2.size() <<std::endl;
                node_info2 = node_info2 + ";" + std::to_string(count2);
                structure_results2 = structure_results2 + ";" + std::to_string(level_data2.size());
                //benchmark_alex_real_seperate(index_original, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
                //benchmark_lipp_real(index, level_data, level_data, "LIPP", data_name+" "+std::to_string(i), poi_thres, original_changed_output);
                //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_smooth, &changed_data);

            //     if (check_duplicates(level_data)) {
            //     std::cout << "The original data has duplicate values." << std::endl;
            // } else {
            //     std::cout << "The original data does not have duplicate values." << std::endl;
            // }
                // std::vector<KEY_TYPE> lookup_keys_level2 = get_search_keys(level_data2, level_data2.size(), 1000000, seed);
                // benchmark_sali_real(index,lookup_keys_level2, lookup_keys_level2, "SALI", dataset_name+"_lookup_level_"+std::to_string(i)+ "_insert_"+std::to_string(j), 
                //     poi_thres, performance_output,poison,insert,orignal_P,insert_threshold);
                
            }
            // benchmark_sali_real(index2,altered_data, altered_data, "SALI", dataset_name+"_changed"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
            //  smooth,insert,orignal_P,insert_threshold);
            benchmark_sali_real(index2,lookup_keys, lookup_keys, "SALI", dataset_name+"_lookup"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
                smooth,insert,orignal_P,insert_threshold);
            benchmark_sali_real(index2,lookup_keys_zipf, lookup_keys_zipf, "SALI", dataset_name+ "_lookup_zipf"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
            smooth,insert,orignal_P,insert_threshold);

            std::vector<KEY_TYPE> insert_vector(insert_indexes2.size());

            //subvector.resize(bulk_load_indexes.size());
            insert_vector.clear();

            std::transform(insert_indexes2.begin(), insert_indexes2.end(), std::back_inserter(insert_vector),
                    [&legitimate_data](int index) { return legitimate_data[index]; });

            benchmark_sali_real(index2,insert_vector, insert_vector, "SALI", dataset_name+ "_insert"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
            smooth,insert,orignal_P,insert_threshold);

            
            structure_results2 = structure_results2 + ";" + node_info2 +";" +std::to_string(duration);
            saveToCSV(structure_results2,structure_output);

            size_t index_size_after2 =  index2.total_size();

            std::string model_results = "SALI;" + dataset_name + ";" + std::to_string(smooth) + ";"+ std::to_string(insert) + ";" +
            std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
            + std::to_string(num_total_smooth) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after2) + ";" +
            std::to_string(total_node_count2) + ";" + std::to_string(total_data_count2) + ";" + std::to_string(tra_time)  ;
            
            std::cout << "Index size " << index_size_after2 <<std::endl;


            
            saveToCSV(model_results,output_model_csv);
        }
        
    }

    std::cout << "Virtual number " << num_total_smooth <<std::endl;
    

    if(benchmark){

        if(!insert){
        
            //=================
             benchmark_sali_real(index2,lookup_changed, lookup_changed, "SALI", data_name+ "_changed", poi_thres, performance_output,
             smooth,insert,orignal_P,insert_threshold);
            benchmark_sali_real(index2,lookup_keys, lookup_keys, "SALI", data_name+"_lookup", poi_thres, performance_output,
            smooth,insert,orignal_P,insert_threshold);
            benchmark_sali_real(index2,lookup_keys_zipf, lookup_keys_zipf, "SALI", data_name+ "_lookup_zipf", poi_thres, performance_output,
            smooth,insert,orignal_P,insert_threshold);
            
        }

    }

    std::cout << "Max model node height  " << max_model_height << std::endl;
    std::cout << "Time  " << tra_time/1000000000.0 << "s" << std::endl;
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;

    std::cout << "END " <<std::endl;
    
}
//FUNCTIONS NEEDED FOR LIPP
std::pair<uint64_t, PAYLOAD_TYPE> * create_values(std::vector<KEY_TYPE> data,int * size){

    std::mt19937_64 gen_payload(std::random_device{}());
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[*size];

    int i = 0;
    
    for (KEY_TYPE key : data) {
                values[i].first = key;
                values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                i++;
        }

    return values;

}

//To create the values from poisoned data
std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_poisoned_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, int size){

    std::mt19937_64 gen_payload(std::random_device{}());

    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[size];
    
        int idx = 0;
        for (int i = 0; i < size; i++) {
            //Check if the current legitimate data value and poisoned data value are the same (meaning this key is an original key)
            if(data_leg[idx] == data_poi[i]){
                //Then use the same payload
                  values[i].second = payload[idx];

                  //Go to the next position in the legitimate keys vector.
                  idx++;
            }
            //If they are different meaning this is a poisoned key. Get a new random payload
            else{
                  values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            }   

            //Either way just use the key as the key from the poisoned vector
             values[i].first = data_poi[i];
        }

  

    return values;

}

bool compare_level_des(const Node* a, const Node* b) {
        return a->level > b->level;
}
void scan_node(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload){
    for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2) {
                        level_data->push_back(node->items[i].comp.data.key);
                        level_payload->push_back(node->items[i].comp.data.value);
                    } 
            }

}
void scan_node_all(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload){
    for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2 || node->items[i].entry_type == 3) {
                        level_data->push_back(node->items[i].comp.data.key);
                        level_payload->push_back(node->items[i].comp.data.value);
                    } 
            }

}
void count_data_node(Node* node, int * count){
    int count_temp = 0;
    for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2) {
                        count_temp++;
                    } 
            }
    *count = count_temp;

}


void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1, std::vector<KEY_TYPE> vec2) {
    // Convert vectors to sets
    std::copy_if(vec1.begin(), vec1.end(), std::back_inserter(*existingVector), [&vec2](KEY_TYPE element) {
        return std::find(vec2.begin(), vec2.end(), element) == vec2.end();
    });
}

void get_nodes_by_level_with_child( std::vector<Node*> *nodes, std::vector<Node*> *nodes_by_level, int level)
{

    (*nodes_by_level).clear();
    for(Node* node : *nodes){
            bool has_child = false;
            //get child nodes and check if they are data nodes

            //Also save models by their level (so I can find the parent much easily)
            if(level != node->level){
                continue;
            }
            for (int i = 0; i < node->num_items; i ++) {
                if (node->items[i].entry_type == 1) {
                    has_child = true;
                    break;
                }
            }
            if(has_child){
                nodes_by_level->push_back(node);
            }
       
        }
}

void get_nodes_by_level( std::vector<Node*> *nodes, std::vector<Node*> *nodes_by_level, int level)
{

    (*nodes_by_level).clear();
    for(Node* node : *nodes){
            //get child nodes and check if they are data nodes

            //Also save models by their level (so I can find the parent much easily)
            if(level == node->level ){
                nodes_by_level->push_back(node);
            }
       
        }
}

void get_poisoned_keys(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, KEY_TYPE* keys, PAYLOAD_TYPE* payloads, int size){

    std::mt19937_64 gen_payload(std::random_device{}());

    int idx = 0;
    for (int i = 0; i < size; i++) {
        //Check if the current legitimate data value and poisoned data value are the same (meaning this key is an original key)
        if(data_leg[idx] == data_poi[i]){
            //Then use the same payload
                payloads[i] = payload[idx];

                //Go to the next position in the legitimate keys vector.
                idx++;
        }
        //If they are different meaning this is a poisoned key. Get a new random payload
        else{
                payloads[i] = static_cast<PAYLOAD_TYPE>(gen_payload());
        }   

        //Either way just use the key as the key from the poisoned vector
            keys[i] = data_poi[i];;
    }


}

void shift_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key) {
    if ((*vec).empty()) {
        return; // If the vector is empty, no need to shift
    }

    // Find the minimum value in the vector
    // Shift all elements by subtracting the minimum value
    KEY_TYPE offset = -key;

    // Add the offset to each element in the vector
    std::transform((*vec).begin(), (*vec).end(), (*vec).begin(), [offset](KEY_TYPE& element) {
        return element + offset;
    });
}

void shift_back_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key) {
    if ((*vec).empty()) {
        return; // If the vector is empty, no need to shift
    }

    // Find the minimum value in the vector
    // Shift all elements by subtracting the minimum value
    KEY_TYPE offset = key;

    // Add the offset to each element in the vector
    
    std::transform((*vec).begin(), (*vec).end(), (*vec).begin(), [offset](KEY_TYPE& element) {
        return element + offset;
    });
}

void get_nodes_children(Node* node, std::vector<Node*> *nodes_by_level)
{

   for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 1) {
                        nodes_by_level->push_back(node->items[i].comp.child);
                       
                    } 
            }
}

//To get max height of a data node
void get_data_node_max_height(Node* node, int * height){
    int count_temp = 0;
    
    typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, node));
        while (!s.empty()) {
            


            int level = s.top().first;
            Node* node = s.top().second;
            
            s.pop();

            if(count_temp < level){
                count_temp = level;
            }
            
            for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2) {
                       
                    } else if (node->items[i].entry_type == 1) {
                        s.push(Segment(level+1, node->items[i].comp.child));
                        
                    }
            }
        }
    
    *height = count_temp;

}

//To get data by level
void get_data_node_data_by_level(Node* node, std::vector<std::vector<KEY_TYPE>*>& node_data){
    //int begin = 0;
    
    typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, node));
        while (!s.empty()) {
            int level = s.top().first;
            Node* node = s.top().second;
            
            s.pop();
            
            for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2) {
                        
                        node_data[level]->push_back(node->items[i].comp.data.key);
                        
                    } else if (node->items[i].entry_type == 1) {
                        s.push(Segment(level+1, node->items[i].comp.child));
                        
                    }
            }
        }

}

//To get data by level
void get_node_data_by_level(Node* node, std::set<KEY_TYPE>& node_data, int level){
    
    typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, node));
        while (!s.empty()) {
            int level = s.top().first;
            Node* node = s.top().second;
            
            s.pop();
            
            for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2) {
                       
                        node_data.insert(node->items[i].comp.data.key);
                        
                     
                    } 
                    
            }
        }

}
//To get data by level
void get_node_data_below_level(Node* node, std::set<KEY_TYPE>& node_data, int level){
    
    typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, node));
        while (!s.empty()) {
            int level2 = s.top().first;
            Node* node = s.top().second;
            
            s.pop();
            
            for (int i = 0; i < node->num_items; i ++) {
                    if (node->items[i].entry_type == 2) {
                       
                        if(level2 > level){
                            node_data.insert(node->items[i].comp.data.key);
                         }
                        
                    } else if (node->items[i].entry_type == 1) {
                        s.push(Segment(level+1, node->items[i].comp.child));
                       
                    }
            }
        }

}


int csv_node(std::vector<Node*> nodes_by_level, KEY_TYPE P_i,size_t start, size_t end){
    int current_poi = 0;
    
    for (size_t inner_idx = start; inner_idx < end; ++inner_idx) {
        
        Node* best_node;
        Node* new_node ;
        
        int node_height_original_=0;
        std::vector<std::vector<KEY_TYPE>*> node_data_original;

        best_node = nodes_by_level[inner_idx];
        double best_children_cost = 0;


        //GET ALL DATA FROM SUBTREE
        //===========================
        int size = best_node->size;
        KEY_TYPE* keys = new KEY_TYPE[size];
        PAYLOAD_TYPE* values = new PAYLOAD_TYPE[size];

        index2.scan_subtree2(best_node, keys, values,false);

        //Shift the data to avoid large values
        std::vector<KEY_TYPE>model_node_data(keys, keys + size);
        KEY_TYPE first_val = model_node_data[0];
        shift_vector(&model_node_data, first_val);

        std::vector<PAYLOAD_TYPE>model_node_payload(values, values + size);

        std::vector<KEY_TYPE>  level_data_node;
        std::vector<PAYLOAD_TYPE>  level_payload_node;

        scan_node_all(best_node, &level_data_node, &level_payload_node);

        if(level_data_node.size()== 0 || model_node_data.size()>1000){
           
            continue;
        }
       

        //smooth using the model data or till the cost is less than the children nodes sum of costs.
        //Vector to hold the poisoned model data or if no poisoning use this as well. 
        std::vector<KEY_TYPE> poisoned_model_data;
        

        double poi_thres = P_i/(1.0*model_node_data.size());
        int poisoned_model_data_size = model_node_data.size();
        
        //If lamdha is more than 0 then smooth
        if(model_node_data.size()>2){
        
            //Perform smoothing
            
            poisoned_model_data = perform_poisoning(model_node_data, poi_thres);

            shift_back_vector(&poisoned_model_data, first_val);
            shift_back_vector(&model_node_data, first_val);

            poisoned_model_data_size = poisoned_model_data.size();

            KEY_TYPE* keys_poi = new KEY_TYPE[poisoned_model_data_size];
            PAYLOAD_TYPE* values_poi = new PAYLOAD_TYPE[poisoned_model_data_size];

            //get new values for the poisoning data while using the old ones for the existing
            get_poisoned_keys(model_node_data, poisoned_model_data, model_node_payload, keys_poi, values_poi,poisoned_model_data_size);
            
            // Save result in a thread-safe manner
            {
                std::lock_guard<std::mutex> lock(resultsMutex);
                values_for_inner[inner_idx] = std::make_tuple(level_data_node[0], level_payload_node[0], keys_poi, values_poi,poisoned_model_data_size,keys, size);

            }
            
            current_poi = poisoned_model_data_size - model_node_data.size();
            
        }
        else{
                {
                KEY_TYPE emptyVec;
                PAYLOAD_TYPE emptyp;
                //std::cout <<  "tuple" << std::endl;
                std::lock_guard<std::mutex> lock(resultsMutex);
                values_for_inner[inner_idx] = std::make_tuple(emptyVec,emptyp,nullptr,nullptr,0,nullptr, 0);

            }
        }
    }

    return current_poi;
            
}

