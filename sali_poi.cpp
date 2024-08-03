// RUN using
// g++ sali_poi.cpp -std=c++17 -o sali_poi -march=native -mpopcnt


//./sali_poi ../../../ext/Data/fb_200M_uint64 50000

// g++ -fopenmp -I /opt/intel/oneapi/tbb/2021.12/include sali_poi.cpp -std=c++17 -o sali_poi -L /opt/intel/oneapi/tbb/2021.12/lib -ltbb -march=native -mpopcnt
// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi ../../../ext/Data/fb_200M_uint64 10000000
// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi ../../../ext/Data/osm_cellids_200M_uint64 1000000
// nohup LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi test 1000000 1 100000 0 &
//  LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH nohup ./sali_poi test 1000000 1 100000 0 &

//  LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi fb 200000000 1 20000000 0
// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi test 1000000 0 100000 0


// I changed regression_benchmark and fast_brute force io helper
// I changed exponentional search to key_type2 from double
// I removed regression_linear and put it in pgm
//I added number_of_searches to alex_nodes.h and then ++ it into the exponential searches (2) also =0 at their top

//Set the sizes appopretly
//Check the number of keys and nodes below X for ALEX and LIPP
//Save all of these to a file
//Save index size, keys, query time, search time, build time, numer of poisoning, save change too.
//Create datasets with write heavy, ready heavy, normal. Do the poisoning after bulk_loading.
//Test lookup, inserts performance.

//Write to a csv file.
//Save the changed data then immediately run the non-poisoned one.
//Reduce to just 2 workloads for now.
//run for 4 datasets. 
//That will be 16 times. or 8 * 700 s + 8*200

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


//UNCOMMENT FOR REAL DATA
//======================
// g++ test_alex.cpp -std=c++17 -o test_alex
// ./test_alex ../../../ext/Data/fb_200M_uint64 1000
#include "src/helpers/io_handler_real.h"
//#include "src/poisoning_new_min_real_new_long.h"
// #include "src/poisoning_new_min_real_new_derv.h"
// #include "src/poisoning_new_min_real_new_long_no_outs.h"
//#include "src/poisoning_new_min_real_max_gap.h"
#include "src/poisoning_distribution.h"
#include "src/poisoning_simple.h"

//#include "src/poisoning_new_min_real_new.h"

//I can delete the unneccessary ones if needed in pgm benchmark

//#include "src/helpers/alex_benchmark_real_new.h"
#include "src/helpers/alex_benchmark_real_new_ind.h"
#include "src/helpers/lipp_benchmark_real_new_ind.h"
#include "src/helpers/sali_benchmark_real_new_ind.h"

//I CHANGED NFL BY ADDING MEAN AND VAR FOR THE NUMERICAL_FLOW THERE. BECUASE IT WAS USING IT WITHOUT CALCULATING IT
// AND GIVING ME A SEGMENTATION ERROR WHEN I USE IT BECAUSE THIS RESULTS IN A 0 FOR ALL.

//NFL check conflicts, numercial_flow and nfl
//Because something is wrong with it, it is not finding the correct key (probably the transformation thing)
//#include "src/helpers/nfl_benchmark_real_new2.h"

#include "src/fast_brute_force_real.h"

#include "src/theil_sen.h"

#include <unordered_set>
#include <unistd.h>
// #include <omp.h>
// #include "src/helpers/competitors/SALI/src/core/sali.h"

//#include "src/helpers/nfl_benchmark.h"

// USE THIS TO RUN THE CODE
// source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include test_alex_no_outs.cpp -std=c++17 -o test_alex_no_outs -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

// source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include benchmark.cpp -std=c++17 -o benchmark -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
// typedef LIPP<KEY_TYPE, PAYLOAD_TYPE, true>::Node Node;
typedef sali::SALI<KEY_TYPE, PAYLOAD_TYPE, true>::Node Node;
//typedef sali::SALI<KEY_TYPE, PAYLOAD_TYPE, true>::Node Node_S;

bool compare_level_des(const Node* a, const Node* b);
void scan_node(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload);
void scan_node_all(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload);

void get_nodes_by_level( std::vector<Node*> *nodes, std::vector<Node*> *nodes_by_level, int level);
void get_poisoned_keys(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi, std::vector<PAYLOAD_TYPE> payload, KEY_TYPE* keys, PAYLOAD_TYPE* payloads, int size);
void get_nodes_by_level_with_child( std::vector<Node*> *nodes, std::vector<Node*> *nodes_by_level, int level);
void shift_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
void shift_back_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
void get_nodes_children(Node* node, std::vector<Node*> *nodes_by_level);
bool compare_descending_index(const int& a, const int& b, double *vec);
void count_data_node(Node* node, int * count);
void normalize_array(double* arr, int size);
void get_data_node_max_height(Node* node, int * height);
void get_data_node_data_by_level(Node* node, std::vector<std::vector<KEY_TYPE>*> & node_data);
void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1,  std::vector<KEY_TYPE> vec2);

std::pair<uint64_t, PAYLOAD_TYPE> * create_values(std::vector<KEY_TYPE> data,int * size);

//To create the values from poisoned data
std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_poisoned_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, int size);


std::vector<KEY_TYPE> sampleRandom(std::vector<KEY_TYPE> data, double percent);

std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_new_rank_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi, std::vector<int>& ranks,
std::vector<PAYLOAD_TYPE> payload, int size);

bool check_duplicates(std::vector<KEY_TYPE> data);

void clearMemoryCache() {
      if (system("sync && sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'") != 0) {
        std::cerr << "Failed to clear memory cache." << std::endl;
    }
}
int main(int argc, char *argv[]){
    clearMemoryCache();

    std::string dataset_name = argv[1];
    std::string poi_size_str = argv[4];
    std::string insert_prop_str = argv[5];
    
    bool benchmark = 1;
    std::string data_folder = "../../../mnt2/Data/";
    //std::string data_output = "../../../mnt2/Data/";
    std::string data_output = data_folder + dataset_name+".bin";

    bool poison = (argv[3][0] == '1');

    

    KEY_TYPE orignal_P = static_cast<KEY_TYPE>(std::stoi(poi_size_str));

    
    double insert_threshold = std::stod(insert_prop_str);
    bool insert = (insert_threshold > 0.0);

    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_output);

    unsigned seed = 123;

    bool poison_distribute = false;
    
    KEY_TYPE P = orignal_P;
    KEY_TYPE P_left = P;

    double constant_cost = 1;
    int pow_val = 1;

    double original_poi_thres = 0;
    std::string method = "nt6";

    if(poison){
        original_poi_thres = 0.5;
        method = "ntp6";
    }
    
    double poi_thres = original_poi_thres;
    long pois_count = poi_thres*legitimate_data.size();
    bool use_new_cost = true;
    int num_total_poisoning = 0;
    int max_check = 10;
    
    
    srand(12345);

    //Get the dataname, size and the method to save them
    std::string data_name = argv[1];
    std::string data_name2 = "osm_200M_uint64";
    
    std::string data_size = argv[2];

    //filenames of the benchmark
    std::string model_output = "results/model_output_real";
    std::string original_output = "results/original_index_full.csv";
    std::string original_changed_output = "results/original_index_changed";
    std::string poisoned_output = "results/poisoned_index_full";
    std::string poisoned_changed_output = "results/poisoned_index_changed.csv";

    std::string output_folder = "../../../mnt2/";
    std::string performance_output = output_folder+"results/sali_100/sali2_original_performance.csv";
    std::string structure_output = output_folder+"results/sali_100/sali2_original_structure.csv";

    if(poison){
        performance_output = output_folder+"results/sali_100/sali2_poisoned_performance.csv";
        structure_output = output_folder+"results/sali_100/sali2_poisoned_structure.csv";
    }

    std::string changed_data_output = "results/changed_data.bin";
    std::string bulk_data_output = "../../../ext/Data/fb_bulk_load_data.bin";
    std::string insert_data_output = "../../../ext/Data/fb_insert_data.bin";

    

    std::string lookups_output = "results/lipp_lookups.bin";
    std::string lookups_zip_output = "results/lipp_lookups_zip.bin";

    std::string output_model_csv = output_folder+"results/sali_100/sali2_model_output_good.csv";
    

    //To save the changed data
    std::vector<KEY_TYPE> changed_data;
    std::vector<KEY_TYPE> old_data;

    sali::SALI<KEY_TYPE, PAYLOAD_TYPE> index;

    //Create the values (append new payloads)
    std::mt19937_64 gen_payload(std::random_device{}());
    int values_size = legitimate_data.size();

    values_size = legitimate_data.size();
    std::cout << "removed duplicates." << values_size << std::endl;
    
    

    //for inserts
    std::vector<int> all_indexes(values_size);
    std::vector<int> bulk_load_indexes;
    std::vector<int> insert_indexes;

    std::iota(all_indexes.begin(), all_indexes.end(), 0);

    auto start_build = std::chrono::high_resolution_clock::now();

    //Depending on the insertion 
    int numberOfIndexes = 0;
    if(insert){
        std::string bulk_index_output = data_folder + "Splits/"+  dataset_name +"_bulk.bin";
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
        
        index.bulk_load(values, numberOfIndexes);
        

    }
    else{
        std::cout << "Creating values" << std::endl;
        auto values = create_values(legitimate_data,&values_size);

        std::cout << "Created values" << std::endl;
       
        index.bulk_load(values, values_size);
        delete[] values; // Delete the dynamically allocated array
        values = nullptr; // Reset the pointer to nullptr

        std::cout << "Created index" << std::endl;
    }
    auto stop_build = std::chrono::high_resolution_clock::now();

    std::vector<int>().swap(all_indexes);
    std::vector<int>().swap(bulk_load_indexes);

   
    std::cout << "============================" << std::endl;

    std::cout << "Created the index" << std::endl;
    std::cout << "============================" << std::endl;

    long index_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_build - start_build).count();
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;

    
    auto start_tra = std::chrono::high_resolution_clock::now();

    int max_model_height =0;

    std::vector<Node*> nodes;
   
    index.scan_nodes(&nodes);

    std::sort(nodes.begin(), nodes.end(), compare_level_des);
    max_model_height = nodes[0]->level;
    

    int replaced_models_number = 0;

    
    std::cout << "DONE: models useful and models by level" << std::endl;
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
   
    //To keep count of the number of models we have tried without use
    int tried = 0;
    int sucess =  0 ;

    int current_count = 0;
    int over_count = 0;


    std::vector<KEY_TYPE>  altered_data;
    int size_to_check;

    if(nodes_by_level.size() <= 0){
        current_level = 0;
    }
    
    while(current_level == 2 && poison){
        Node* best_node;
        Node* new_node ;
        bool model_converted = false;
        inner_idx = 0;
    
        size_to_check = nodes_by_level.size();
        int original_size_to_check = size_to_check;

        tried = 0;
        sucess =  0; 
        
         while(inner_idx < size_to_check && tried < max_check  ){

            if(inner_idx%10000 ==0){
                std::cout << "starting " << inner_idx<<std::endl;
            }

            int node_height_original_=0;
            std::vector<std::vector<KEY_TYPE>*> node_data_original;
    
            best_node = nodes_by_level[inner_idx];
            double best_children_cost = 0;

            //GET ALL DATA FROM SUBTREE
            //===========================
            int size = best_node->size;
            KEY_TYPE* keys = new KEY_TYPE[size];
            PAYLOAD_TYPE* values = new PAYLOAD_TYPE[size];

            index.scan_subtree2(best_node, keys, values,false);

            std::vector<KEY_TYPE>model_node_data(keys, keys + size);
            KEY_TYPE first_val = model_node_data[0];
            shift_vector(&model_node_data, first_val);

            std::vector<PAYLOAD_TYPE>model_node_payload(values, values + size);

            std::vector<KEY_TYPE>  level_data_node;
            std::vector<PAYLOAD_TYPE>  level_payload_node;

            scan_node(best_node, &level_data_node, &level_payload_node);
            
            if(level_data_node.size()== 0 || model_node_data.size()>1000){
                 get_nodes_children(best_node, &nodes_by_level);
                 over_count++;
                 size_to_check = nodes_by_level.size();
                 inner_idx++;
                 tried++;
                continue;
            }

            
            //Vector to hold the poisoned model data or if no poisoning use this as well.
            std::vector<KEY_TYPE> poisoned_model_data;
            
            KEY_TYPE P_i;
            if(inner_idx < original_size_to_check){
                P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(original_size_to_check  - current_count))));
            }
            else{
                P_i = std::max(static_cast<KEY_TYPE>(model_node_data.size()*0.1),
                static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current_count))))) ;
            }
            current_count++;
         
            
            

            poi_thres = P_i/(1.0*model_node_data.size());
            int cutoff = -100 ;
            int poisoned_model_data_size = model_node_data.size();
            int current_poi =0;

            sucess++;
            if(P_i > 0){
                
                poisoned_model_data = perform_poisoning(model_node_data, poi_thres);
                shift_back_vector(&poisoned_model_data, first_val);
                shift_back_vector(&model_node_data, first_val);

                poisoned_model_data_size = poisoned_model_data.size();

                KEY_TYPE* keys_poi = new KEY_TYPE[poisoned_model_data_size];
                PAYLOAD_TYPE* values_poi = new PAYLOAD_TYPE[poisoned_model_data_size];

                //get new values for the poisoning data while using the old ones for the existing
                get_poisoned_keys(model_node_data, poisoned_model_data, model_node_payload, keys_poi, values_poi,poisoned_model_data_size);
                
                P_left -= (poisoned_model_data_size -  model_node_data.size());
                
                //Cost condition
                if(poisoned_model_data_size > model_node_data.size()){
                    replaced_models_number++;

                    //Getting the original data in each level of it
                    get_data_node_max_height(best_node, &node_height_original_);
                   
                    node_data_original.resize(node_height_original_+1);

                    for(int i = 0; i <= node_height_original_; i++){
                        node_data_original[i] = new std::vector<KEY_TYPE>();
                    }
                  
                    
                    get_data_node_data_by_level(best_node, node_data_original);
                    
                    new_node = index.poison_bulk_load_real(index.root, best_node,level_data_node[0], level_payload_node[0], keys_poi, values_poi,poisoned_model_data_size,keys, size);

                }


                current_poi = poisoned_model_data_size - model_node_data.size();
            }
            else{
              
            }
            
            int best_level = best_node->level;
            int best_node_idx = 0;
            bool parent_found = false;
            int parent_model_idx = 0;

                 // Getting the changed data details
                if(poisoned_model_data_size > model_node_data.size()){

                    //Getting new data after changing in each level.
                    int node_height_changed = 0;
                    get_data_node_max_height(new_node, &node_height_changed);
                    std::vector<std::vector<KEY_TYPE>*> node_data_changed;
                    node_data_changed.resize(node_height_changed+1);
                    
                    

                    for(int i = 0; i <= node_height_changed; i++){
                        node_data_changed[i] = new std::vector<KEY_TYPE>();
                    }
                    
                    get_data_node_data_by_level(new_node, node_data_changed);

                    int start = 1;
                    if(inner_idx >= original_size_to_check){
                        start = 0;
                    }
                    for(int i = start; i <= node_height_original_; i++ ){
                        //If the node is higher than the current
                        if(i > node_height_changed){
                            // std::cout << "Inside " << i << std::endl;
                            (altered_data).insert((altered_data).end(), (*node_data_original[i]).begin(), (*node_data_original[i]).end());
                        }
                        else{
                            find_difference(&altered_data,*node_data_original[i], *node_data_changed[i]);
                        }
                      
                    }



                    for (const auto& element : model_node_data) {
                // Check if the element is not in the same level node (meaning it is from the bottom)
                        if (std::find(level_data_node.begin(), level_data_node.end(), element) == level_data_node.end()) {
                            // If not found, add it to the original vector
                            old_data.push_back(element);
                        }
                
                    }
                }
              
            tried = 0;

            model_converted = true;

            
            num_total_poisoning += current_poi;
            //If model is converted then set tried back to 0 and increase the inner_idx
            if(model_converted){
                inner_idx++;
                tried = 0;
            }
        }
        
        current_level--;
        if(current_level  > 1){
        }
        
        std::cout << "=============="  <<std::endl;
        std::cout << "Done: Iteration " << current_level <<std::endl;
            
    }
    

    int max_model_height_fin = 0;
    
    //GET NEW NODE STATS 
                //=====================

    auto stop_tra = std::chrono::high_resolution_clock::now();

    
    std::cout << " Number " << index.root->size <<std::endl;
     std::cout << "Replaced number " << replaced_models_number <<std::endl;

    
    std::cout << "Poisoned number " << num_total_poisoning <<std::endl;
    std::cout << "Size to check " << size_to_check << std::endl;
    long tra_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra - start_tra).count();

    
    if(poison){
        save_data_bin(altered_data,changed_data_output);
    }
    if(!poison&benchmark){
        altered_data.clear();
        altered_data = read_data_bin(changed_data_output);
    }

    if(insert){
   
        for(int j = 0; j < 5; j++){

        
            std::string insert_index_output = data_folder+ "Splits/"+ data_name +"_insert"+"_"+std::to_string(j)+".bin";
            std::vector<int> insert_indexes2 = readIndexesFromFile(insert_index_output);
            for(int i : insert_indexes2){
                index.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
            }

            benchmark_sali_real(index,altered_data, altered_data, "SALI", dataset_name+ "_changed", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
        }
    }

    //TO BENCHMARK LEVELS
    //index.insert(5,10);
    std::vector<Node*> nodes2;
    index.scan_nodes(&nodes2);
    std::sort(nodes2.begin(), nodes2.end(), compare_level_des);
    max_model_height = nodes2[0]->level;

    size_t total_node_count = 0;
    size_t total_data_count = 0;

    std::string structure_results = "SALI;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning)+ ";Data";;

        // std::cout << test_output <<std::endl; 
    std::string node_info = "Nodes";
    
    std::vector<KEY_TYPE> lookup_keys_level_prop_lower;
    std::vector<KEY_TYPE> lookup_keys_level_prop;
    int level_one_size = 0;
    int level_two_size = 0;
    // for(int i = max_model_height; i >0 ; i-- ){
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

    structure_results = structure_results + ";" + node_info;

    saveToCSV(structure_results,structure_output);

    std::cout << "Poisoned number " << num_total_poisoning <<std::endl;

    //GETTING THE LOOK UPS

    std::vector<KEY_TYPE> lookup_keys = get_search_keys(legitimate_data, values_size, 1000000, seed);

    std::vector<KEY_TYPE> lookup_keys_zipf = get_search_keys_zipf(legitimate_data, values_size, 1000000, seed);
    

    if(benchmark){
         
        size_t index_size_after =  index.total_size();
        std::ofstream file;
        file.open(model_output+".txt", std::ios_base::app);
    
        file << "SALI" << ";" << data_name << ";"  << insert << ";" << insert_threshold << ";"<< legitimate_data.size() << ";" << orignal_P << ";"
        << num_total_poisoning << ";"<< max_model_height << ";" << index_size_after << ";" <<
        total_node_count << ";" << total_data_count << ";" << std::endl;

        file.close();

         std::string model_results = "SALI;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" +
        std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        saveToCSV(model_results,output_model_csv);

    
    if(!insert){

    
        //=================
         benchmark_sali_real(index,altered_data, altered_data, "SALI", data_name+ "_changed", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index,lookup_keys, lookup_keys, "SALI", data_name+"_lookup", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index,lookup_keys_zipf, lookup_keys_zipf, "SALI", data_name+ "_lookup_zipf", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);

         
        //  benchmark_sali_real(index,lookup_keys_level_prop, lookup_keys_level_prop, "SALI", data_name+ "_level_prop", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_sali_real(index,lookup_keys_level_prop_lower, lookup_keys_level_prop_lower, "SALI", data_name+ "_level_prop_lower", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
         
    }

    }

    std::cout << "Max model node height  " << max_model_height << std::endl;
    std::cout << "Time  " << tra_time/1000000000.0 << "s" << std::endl;
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;

    std::cout << "END " <<std::endl;
    
}
//FUNCTIONS NEEDED FOR LIPP
std::pair<uint64_t, PAYLOAD_TYPE> * create_values(std::vector<KEY_TYPE> data,int * size){
    //std::sort(data.begin(), data.end());

    std::mt19937_64 gen_payload(std::random_device{}());
    // std::set<KEY_TYPE> uniqueKeys(data.begin(), data.end());
    // *size = uniqueKeys.size();

    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[*size];

    int i = 0;
    // for (KEY_TYPE key : uniqueKeys) {
    //             values[i].first = key;
    //             values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    //             i++;
    //     }

    for (KEY_TYPE key : data) {
                values[i].first = key;
                values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                i++;
        }

    // KEY_TYPE prev = data[0];
    // values[0].first = data[0];
    // values[0].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    // //std::mt19937_64 gen_payload(std::random_device{}());
    // int dup_count = 0;

    //     for (int i = 1; i < *size; i++) {
    //         if(data[i] != prev){
    //             values[i].first = data[i];
    //             values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    //             prev = data[i];
    //         }
    //         else{
    //             dup_count++;
    //         }
            
    //     }
    //     auto nonZeroValues = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[*size - dup_count];
    //     int index = 0;
    //     for (int i = 0; i < *size; ++i) {
    //         if (values[i].first != 0) {
    //             nonZeroValues[index++] = values[i];
    //         }
    //     }
    //     //std::cout << "Size " << *size << std::endl;
    //     *size = *size -  dup_count; 
        //std::cout << "Size " << dup_count << std::endl;
        //std::cout << "Size " << *size << std::endl;
    

    return values;

}

//To create the values from poisoned data
std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_poisoned_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, int size){

    std::mt19937_64 gen_payload(std::random_device{}());

    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[size];
    //int index = 0;
    //std::mt19937_64 gen_payload(std::random_device{}());

        // for (int i = 0; i < size; i++) {

        //     auto it = std::find(data_leg.begin(), data_leg.end(), data_poi[i]);

        //      values[i].first = data_poi[i];

        //     if (it != data_leg.end()) {
        //         values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        //     } else {
        //         index = std::distance(data_leg.begin(), it);
        //         values[i].second = payload[index];
        //     }
           
            
        // }
        int idx = 0;
        for (int i = 0; i < size; i++) {

            //auto it = std::find(data_leg.begin(), data_leg.end(), data_poi[i]);
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

            // if (it != data_leg.end()) {
            //     values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            // } else {
            //     index = std::distance(data_leg.begin(), it);
            //     values[i].second = payload[index];
            // }
           
            
        }

  

    return values;

}

bool compare_level_des(const Node* a, const Node* b) {
        return a->level > b->level;
}
void scan_node(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload){
    //int begin = 0;
    for (int i = 0; i < node->num_items; i ++) {
                //if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (node->items[i].entry_type == 2) {
                        level_data->push_back(node->items[i].comp.data.key);
                        level_payload->push_back(node->items[i].comp.data.value);
                       // begin ++;
                    } 
                //}
            }

}
void scan_node_all(Node* node, std::vector<KEY_TYPE> * level_data, std::vector<PAYLOAD_TYPE> * level_payload){
    //int begin = 0;
    for (int i = 0; i < node->num_items; i ++) {
                //if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (node->items[i].entry_type == 2 || node->items[i].entry_type == 3) {
                        level_data->push_back(node->items[i].comp.data.key);
                        level_payload->push_back(node->items[i].comp.data.value);
                       // begin ++;
                    } 
                //}
            }

}
void count_data_node(Node* node, int * count){
    //int begin = 0;
    int count_temp = 0;
    for (int i = 0; i < node->num_items; i ++) {
                //if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (node->items[i].entry_type == 2) {
                        count_temp++;
                       // begin ++;
                    } 
                //}
            }
    *count = count_temp;

}

//To get max height of a data node
void get_data_node_max_height(Node* node, int * height){
    //int begin = 0;
    int count_temp = 0;
    // for (int i = 0; i < node->num_items; i ++) {
    //             if (BITMAP_GET(node->none_bitmap, i) == 0) {
    //                 if (BITMAP_GET(node->child_bitmap, i) == 1) {
    //                     count_temp++;
    //                    // begin ++;
    //                 } 
    //             }
    //         }

    typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, node));
        while (!s.empty()) {
            


            int level = s.top().first;
            Node* node = s.top().second;
            //int level = s.top().
            //const int SHOULD_END_POS = begin + node->size;
            s.pop();

            if(count_temp < level){
                count_temp = level;
            }
            //node_data[level] = new std::vector<KEY_TYPE>();
            //std::cout << "Here scan1 " << begin <<" " << SHOULD_END_POS<< std::endl;
            for (int i = 0; i < node->num_items; i ++) {
                //if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (node->items[i].entry_type == 2) {
                        // keys[begin] = node->items[i].comp.data.key;
                        // values[begin] = node->items[i].comp.data.value;
                        //std::cout << "Here scan1 " << level << std::endl;
                        //node_data[level]->push_back(node->items[i].comp.data.key);
                        
                        //begin ++;
                    } else if (node->items[i].entry_type == 1) {
                        s.push(Segment(level+1, node->items[i].comp.child));
                        //begin += node->items[i].comp.child->size;
                        //std::cout << "Here scan2 " << level+1 << std::endl;
                    }
                //}
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
            //int level = s.top().
            //const int SHOULD_END_POS = begin + node->size;
            s.pop();
            
            //std::cout << "Here scan1 " << begin <<" " << SHOULD_END_POS<< std::endl;
            for (int i = 0; i < node->num_items; i ++) {
                //if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (node->items[i].entry_type == 2) {
                        // keys[begin] = node->items[i].comp.data.key;
                        // values[begin] = node->items[i].comp.data.value;
                        //std::cout << "Here scan1 " << level << " "<< node->items[i].comp.data.key << std::endl;
                        node_data[level]->push_back(node->items[i].comp.data.key);
                        
                        //begin ++;
                    } else if (node->items[i].entry_type == 1) {
                        s.push(Segment(level+1, node->items[i].comp.child));
                        //begin += node->items[i].comp.child->size;
                        //std::cout << "Here scan3 " << level+1 << std::endl;
                    }
                //}
            }
        }

}

void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1, std::vector<KEY_TYPE> vec2) {
    // Convert vectors to sets
    //std:: cout << "Inside diff" <<std::endl;

    // for (KEY_TYPE data : vec1) {
    //                             std::cout << data << " ";
    //                         }
    // for (KEY_TYPE data : vec2) {
    //                             std::cout << data << " ";
    //                         }
    // std::set<KEY_TYPE> set1(vec1.begin(), vec1.end());
    // std::set<KEY_TYPE> set2(vec2.begin(), vec2.end());

    // std::vector<KEY_TYPE> difference;
    
    // // Find elements in set1 but not in set2
    // std::set_difference(set1.begin(), set1.end(),
    //                     set2.begin(), set2.end(),
    //                     std::inserter((*existingVector), (*existingVector).end()));

    // return difference;
    // std:: cout << "Inside diff" <<std::endl;
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
            // if(level == node->level && (*node->child_bitmap & 0xFF) ){
            //     nodes_by_level->push_back(node);
            // }
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

    //auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[size];
    //int index = 0;
    //std::mt19937_64 gen_payload(std::random_device{}());

        // for (int i = 0; i < size; i++) {

        //     auto it = std::find(data_leg.begin(), data_leg.end(), data_poi[i]);

        //      values[i].first = data_poi[i];

        //     if (it != data_leg.end()) {
        //         values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        //     } else {
        //         index = std::distance(data_leg.begin(), it);
        //         values[i].second = payload[index];
        //     }
           
            
        // }
        int idx = 0;
        for (int i = 0; i < size; i++) {

            //auto it = std::find(data_leg.begin(), data_leg.end(), data_poi[i]);
            //Check if the current legitimate data value and poisoned data value are the same (meaning this key is an original key)
            if(data_leg[idx] == data_poi[i]){
                //Then use the same payload
                  //values[i].second = payload[idx];
                  payloads[i] = payload[idx];

                  //Go to the next position in the legitimate keys vector.
                  idx++;
            }
            //If they are different meaning this is a poisoned key. Get a new random payload
            else{
                  //values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                  payloads[i] = static_cast<PAYLOAD_TYPE>(gen_payload());
            }   

            //Either way just use the key as the key from the poisoned vector
            // values[i].first = data_poi[i];
             keys[i] = data_poi[i];;

            // if (it != data_leg.end()) {
            //     values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            // } else {
            //     index = std::distance(data_leg.begin(), it);
            //     values[i].second = payload[index];
            // }
           
            
        }


}
std::vector<KEY_TYPE> sampleRandom(std::vector<KEY_TYPE> data, double percent) {
    std::vector<KEY_TYPE> shuffledData = data; // Copy the original data
    std::random_device rd;
    std::mt19937 gen(rd());
    if(data.size()>1000){
        std::cout<< "SAMPLING NOW" <<std::endl;
        std::shuffle(shuffledData.begin(), shuffledData.end(), gen); // Shuffle the data

        // Calculate the number of elements to sample (10% of the size)
        std::size_t sampleSize = int(shuffledData.size() * percent);

        return {shuffledData.begin(), shuffledData.begin() + sampleSize};
    }
    else{
        return data;

    }
    

    
}

std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_new_rank_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi, std::vector<int>& ranks,
std::vector<PAYLOAD_TYPE> payload, int size){

    std::mt19937_64 gen_payload(std::random_device{}());

    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[size];
    //int index = 0;
    //std::mt19937_64 gen_payload(std::random_device{}());

        // for (int i = 0; i < size; i++) {

        //     auto it = std::find(data_leg.begin(), data_leg.end(), data_poi[i]);

        //      values[i].first = data_poi[i];

        //     if (it != data_leg.end()) {
        //         values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        //     } else {
        //         index = std::distance(data_leg.begin(), it);
        //         values[i].second = payload[index];
        //     }
           
            
        // }
        int idx = 0;
        for (int i = 0; i < data_poi.size(); i++) {

            //auto it = std::find(data_leg.begin(), data_leg.end(), data_poi[i]);
            //Check if the current legitimate data value and poisoned data value are the same (meaning this key is an original key)
            if(data_leg[idx] == data_poi[i]){
                //Then use the same payload
                  values[idx].first = data_poi[i];
                  values[idx].second = payload[idx];
                  ranks.push_back(i);
                idx++;
                  //Go to the next position in the legitimate keys vector.
                  
            }
            //If they are different meaning this is a poisoned key. Get a new random payload
            // else{
            //       values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            // }   

            //Either way just use the key as the key from the poisoned vector
             

            // if (it != data_leg.end()) {
            //     values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            // } else {
            //     index = std::distance(data_leg.begin(), it);
            //     values[i].second = payload[index];
            // }
           
            
        }

    return values;

}
//======

void shift_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key) {
    if ((*vec).empty()) {
        return; // If the vector is empty, no need to shift
    }

    // Find the minimum value in the vector
    //KEY_TYPE minVal = vec[0];

    // Shift all elements by subtracting the minimum value
    KEY_TYPE offset = -key;

    // Add the offset to each element in the vector
    // for (KEY_TYPE& num : (*vec)) {
    //     num += offset;
    // }
    std::transform((*vec).begin(), (*vec).end(), (*vec).begin(), [offset](KEY_TYPE& element) {
        return element + offset;
    });
}

void shift_back_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key) {
    if ((*vec).empty()) {
        return; // If the vector is empty, no need to shift
    }

    // Find the minimum value in the vector
    //KEY_TYPE minVal = vec[0];

    // Shift all elements by subtracting the minimum value
    KEY_TYPE offset = key;

    // Add the offset to each element in the vector
    // for (KEY_TYPE& num : (*vec)) {
    //     num += offset;
    // }
    std::transform((*vec).begin(), (*vec).end(), (*vec).begin(), [offset](KEY_TYPE& element) {
        return element + offset;
    });
}

void get_nodes_children(Node* node, std::vector<Node*> *nodes_by_level)
{

   for (int i = 0; i < node->num_items; i ++) {
                //if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (node->items[i].entry_type == 1) {
                        nodes_by_level->push_back(node->items[i].comp.child);
                        //level_payload->push_back(node->items[i].comp.data.value);
                       // begin ++;
                    } 
                //}
            }
}

bool compare_descending_index(const int& a, const int& b, double *vec) {
    return vec[a] < vec[b]; // Sort in descending order based on value
}

void normalize_array(double* arr, int size) {
    // Find the minimum and maximum values in the array
    double min_val = *std::min_element(arr, arr + size);
    double max_val = *std::max_element(arr, arr + size);
    
    // Calculate the range of values
    double range = max_val - min_val;
    
    // Normalize each element of the array
    for (int i = 0; i < size; ++i) {
        arr[i] = (arr[i] - min_val) / range;
    }
}

