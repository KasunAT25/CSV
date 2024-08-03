// RUN using
// g++ sali_poi.cpp -std=c++17 -o sali_poi -march=native -mpopcnt


//./sali_poi ../../../ext/Data/fb_200M_uint64 50000

// g++ -fopenmp -I /opt/intel/oneapi/tbb/2021.12/include sali_poi2.cpp -std=c++17 -o sali_poi2 -L /opt/intel/oneapi/tbb/2021.12/lib -ltbb -march=native -mpopcnt
// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi ../../../ext/Data/fb_200M_uint64 10000000
// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi ../../../ext/Data/osm_cellids_200M_uint64 1000000
// nohup LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi2 test 1000000 1 100000 0 &
//  LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH nohup ./sali_poi2 test 1000000 1 100000 0 &

//  LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi2 fb 200000000 0 10000000 1

// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./sali_poi2 test 1000000 0 100000 1


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
#include "src/helpers/competitors/PGM/include/pgm/pgm_index_dynamic.hpp"
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

typedef alex::AlexNode<KEY_TYPE, PAYLOAD_TYPE> node_type;

typedef alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE> model_node_type;
  //This is the data node
typedef alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> data_node_type;

bool check_duplicates(std::vector<KEY_TYPE> data);

std::vector<std::pair<model_node_type*,double>> calculate_best_node_cost_prev(std::vector<model_node_type*>& model_nodes_useful);

void calculate_children_new_cost(std::vector<model_node_type*>& model_nodes_useful,alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index);
void calculate_children_new_cost(std::vector<model_node_type*>& model_nodes_useful, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);

void calculate_data_node_new_cost(std::vector<data_node_type*>& data_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);
void calculate_data_node_new_cost_new(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);

std::vector<double>  calculate_model_node_new_cost(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);

std::vector<double> calculate_model_node_poisoning(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double pois_thres,double search_const, double traversal_const);

std::vector<std::pair<model_node_type*,double>> calculate_best_node_cost_reduction(std::vector<model_node_type*>& model_nodes_useful, 
std::vector<double> new_costs,
std::vector<double> pois_costs );

void get_constants(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<KEY_TYPE> data, int max_data, double *const_search, double *const_traversal);

void get_all_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index, std::vector<node_type*> *nodes,std::vector<model_node_type*> *model_nodes,
std::vector<data_node_type*> *data_nodes,int *max_model_height);

bool check_all_data_nodes(model_node_type* node);

void get_useful_model_nodes( std::vector<model_node_type*> *model_nodes, std::vector<model_node_type*> *model_nodes_useful, 
std::vector<model_node_type*> *model_nodes_by_level);

void get_children_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,model_node_type* model, std::vector<KEY_TYPE> *model_node_data, 
std::vector<PAYLOAD_TYPE> *model_node_payload, double *children_cost,
double const_search, double const_traversal, bool calculate_childern_costs);

std::pair<uint64_t, PAYLOAD_TYPE> * create_values(std::vector<KEY_TYPE> data,int * size);

//To create the values from poisoned data
std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_poisoned_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, int size);

data_node_type * create_new_data_node(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index , model_node_type* model, std::pair<uint64_t, double> * values, std::vector<KEY_TYPE> leg,int size, int cur_size,
double const_search, double const_traversal);

std::vector<model_node_type*> find_parent(model_node_type* model,std::vector<model_node_type*> *model_nodes_by_level, int *model_idx
,int *parent_model_idx);

void update_data_structure(model_node_type* best_node,std::vector<model_node_type*> *model_nodes,
std::vector<model_node_type*> *model_nodes_by_level, std::vector<model_node_type*>* model_nodes_useful,
std::vector<model_node_type*>* model_nodes_useful_by_level, std::vector<model_node_type*> parents,
alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const,std::vector<model_node_type*>* new_useful_parents);

std::vector<model_node_type*> replace_parents(model_node_type* model,data_node_type *new_data_node, std::vector<model_node_type*> *model_nodes_by_level, int *model_idx
,int *parent_model_idx);

bool compare_models(const model_node_type* a, const model_node_type* b);
bool compare_models_old(const model_node_type* a, const model_node_type* b);
bool compare_models_new(const model_node_type* a, const model_node_type* b);

void calculate_cost_difference(std::vector<model_node_type*>& model_nodes_useful);

//void get_children_cost(model_node_type* model, double *children_cost);
void get_children_cost(model_node_type* model, double *children_cost,alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);

void get_useful_model_nodes_by_level( std::vector<model_node_type*> *model_nodes, std::vector<model_node_type*> *model_nodes_useful, int level);

void poison_data_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal, int *poisoned_count,std::vector<KEY_TYPE>* changed_data, double percentile);

void poison_all_data_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal,int *poisoned_count,std::vector<KEY_TYPE>* changed_data);

bool check_parents_removed_correct(model_node_type* model,data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level);
bool check_parents_added_correct(model_node_type* model,data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level);
void replace_parents_data_node(data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level, int *model_idx,
int *parent_model_idx);

bool compare_searches(const data_node_type* a, const data_node_type* b);
void get_data_nodes_by_level( std::vector<data_node_type*> *data_nodes, std::vector<data_node_type*> *data_nodes_by_level, 
int level);

std::vector<model_node_type*> can_be_done(std::vector<model_node_type*>& model_nodes_useful);
//void shiftVector(std::vector<KEY_TYPE>& vec);
void get_data_by_level(std::vector<data_node_type*> *data_nodes_by_level, std::vector<KEY_TYPE> * data);

void calculate_expected_searches_all_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal);

std::vector<KEY_TYPE> sampleRandom(std::vector<KEY_TYPE> data, double percent);

std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_new_rank_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi, std::vector<int>& ranks,
std::vector<PAYLOAD_TYPE> payload, int size);

bool compare_searches_des(const data_node_type* a, const data_node_type* b);

void clearMemoryCache() {
      if (system("sync && sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'") != 0) {
        std::cerr << "Failed to clear memory cache." << std::endl;
    }
}

int main(int argc, char *argv[]){

    clearMemoryCache();
//     sali::SALI<int, int> sali;

//   int key_num = 1000;
//   std::pair<int, int> *keys = new std::pair<int, int>[key_num];
//   for (int i = 0; i < 1000; i++) {
//     keys[i]={i,i};
//   }
//   sali.bulk_load(keys, 1000);

//   omp_set_num_threads(12);

// #pragma omp parallel for schedule(static, 12)
//   for (int i = 1000; i < 2000; i++) {
//     sali.insert(i,i);
//   }
// #pragma omp parallel for schedule(static, 12)
//   for (int i = 0; i < 2000; i++) {
//     std::cout<<"value at "<<i<<": "<<sali.at(i,i)<<std::endl;
//   }

    // bool benchmark = true;
    // std::string data_output = "../../../ext/Data/fb.bin";

    // bool poison = true;
    // KEY_TYPE orignal_P = 200000;

    // bool insert = true;
    // double insert_threshold = 0.5;
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

   // std::vector<KEY_TYPE> legitimate_data = parse_arguments(argc, argv);
    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_output);

    //shiftVector(legitimate_data);

    // for(int i =0; i< legitimate_data.size();i++){
    //     std::cout<< legitimate_data[i] << " , ";
    // }

    //save_data(legitimate_data, "legitimate_data.txt");
    
    


    //PARAMETERS
    //===============
    //REMEBER TO CHANGE THE INDEX FILE NAMES
    // TO DISABLE ADJUSTMENT
    //ENABLE insert
    //FOR OSM CHANGE CREATE_VALUES
   
    unsigned seed = 123;
    
   

    //bool other = false;
    

    bool poison_distribute = false;
    
    KEY_TYPE P = orignal_P;
    // KEY_TYPE P = 120856;
    //KEY_TYPE P = 20000000;
    KEY_TYPE P_left = P;

    double constant_cost = 1;
    int pow_val = 1;

    double original_poi_thres = 0;
    std::string method = "nt6";

    if(poison){
        original_poi_thres = 0.5;
        method = "ntp6";
    }
    
    //double original_poi_thres = 0;
    //std::string method = "nt";
    double poi_thres = original_poi_thres;
    long pois_count = poi_thres*legitimate_data.size();
    bool use_new_cost = true;
    int num_total_poisoning = 0;
    int max_check = 10;
    
    
     //shift_vector(&legitimate_data, legitimate_data[0]);
     //save_data(legitimate_data, "legitimate_data.txt");
     //double derv = 0;
     //perform_poisoning2(legitimate_data, 20, &derv);
     //std::cout<< derv <<std::endl;
    
    srand(12345);
   
   //Check if the dataset has any duplicates
    // if (check_duplicates(legitimate_data)) {
    //     std::cout << "The original data has duplicate values." << std::endl;
    // } else {
    //     std::cout << "The original data does not have duplicate values." << std::endl;
    // }
    // orignal_P = 5;
    // poi_thres = orignal_P/(1.0*legitimate_data.size());
    //poi_thres = 0.2;

    //perform_poisoning(legitimate_data, poi_thres);

    //Get the dataname, size and the method to save them
    std::string data_name = argv[1];
    std::string data_name2 = "osm_200M_uint64";
    
    // std::string method = "their_cost";
    // std::string method = "thier_cost_size";
    std::string data_size = argv[2];

    //filenames of the benchmark
    std::string model_output = "results/model_output_real";
    std::string original_output = "results/original_index_full.csv";
    std::string original_changed_output = "results/original_index_changed";
    std::string poisoned_output = "results/poisoned_index_full";
    std::string poisoned_changed_output = "results/poisoned_index_changed.csv";

    // std::string performance_output = "results/sali_original_performance.csv";
    // std::string structure_output = "results/sali_original_structure.csv";

    // if(poison){
    //     performance_output = "results/sali_poisoned_performance.csv";
    //     structure_output = "results/sali_poisoned_structure.csv";
    // }
    std::string output_folder = "../../../mnt2/";
    std::string performance_output = output_folder+"results/sali3/insert/sali2_insert_original_performance.csv";
    std::string structure_output = output_folder+"results/sali3/insert/sali2_insert_original_structure.csv";

    if(poison){
        performance_output = output_folder+"results/sali3/insert/sali2_insert_poisoned_performance.csv";
        structure_output = output_folder+"results/sali3/insert/sali2_insert_poisoned_structure.csv";
    }
 


    //original_output = original_output+"_"+method+"_"+data_name2+"_"+data_size;
    //original_changed_output = original_changed_output+ "_"+method+"_"+data_name2+"_"+data_size;

    //poisoned_output = poisoned_output+"_"+method+"_"+data_name2+"_"+data_size;
    //poisoned_changed_output = poisoned_changed_output+ "_"+method+"_"+data_name2+"_"+data_size;

    // std::string bulk_index_output = "data/bulk_load_indexes_05.bin";
    // std::string insert_index_output = "data/insert_indexes_05.bin";
    //std::string bulk_index_output = "data/bulk_load_indexes_08.bin";
    //std::string insert_index_output = "data/insert_indexes_08.bin";

    std::string changed_data_output = "results/changed_data.bin";
    std::string bulk_data_output = "../../../ext/Data/fb_bulk_load_data.bin";
    std::string insert_data_output = "../../../ext/Data/fb_insert_data.bin";

    

    std::string lookups_output = "results/lipp_lookups.bin";
    std::string lookups_zip_output = "results/lipp_lookups_zip.bin";

    std::string output_model_csv = output_folder+"results/sali3/insert/sali2_insert_model_output_good.csv";
    

    //To save the changed data
    std::vector<KEY_TYPE> changed_data;
    std::vector<KEY_TYPE> old_data;

    //To save two Alex indexes for testing (original and edited) I only actually need the index
    sali::SALI<KEY_TYPE, PAYLOAD_TYPE> index;
    //sali::SALI<KEY_TYPE, PAYLOAD_TYPE> index_original;

    //sali::SALI<KEY_TYPE, PAYLOAD_TYPE> sali;

    // pgm::DynamicPGMIndex<KEY_TYPE, PAYLOAD_TYPE,index_pgm> index_pgm;
    // alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index_alex;

    //BULK LOAD IT (std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()], data.size())
    //Create the values (append new payloads)
    std::mt19937_64 gen_payload(std::random_device{}());
    int values_size = legitimate_data.size();
    // auto values = create_values(legitimate_data,&values_size);

    //std::set<KEY_TYPE> uniqueKeys(legitimate_data.begin(), legitimate_data.end());
    //values_size = uniqueKeys.size();

    //legitimate_data.clear();
    //legitimate_data.assign(uniqueKeys.begin(), uniqueKeys.end());

    //CAN REMOVE FOR FB
    //std::sort(legitimate_data.begin(), legitimate_data.end());

    // legitimate_data.erase(std::unique(legitimate_data.begin(), legitimate_data.end()), legitimate_data.end());
    // save_data_bin(legitimate_data,data_output);

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
        std::cout << "Insert." << values_size << std::endl;
        //Instead of this, maybe have a file written.
    //     std::random_device rd;
    //     std::mt19937 gen(seed);

    //     std::shuffle(all_indexes.begin(), all_indexes.end(), gen);
    //     std::cout << "Shuffuled." << values_size << std::endl;

    //     numberOfIndexes = values_size*(1-insert_threshold);
    //     bulk_load_indexes.resize(numberOfIndexes);
    //     insert_indexes.resize(values_size - numberOfIndexes);
       
    //    //Saving part
    //     bulk_load_indexes.assign(all_indexes.begin(), all_indexes.begin() + numberOfIndexes);
    //     std::sort(bulk_load_indexes.begin(), bulk_load_indexes.end());
    //      std::cout << "Shuffuled2." << values_size << std::endl;
    //     insert_indexes.assign(all_indexes.begin() + numberOfIndexes, all_indexes.end());
    //     std::sort(insert_indexes.begin(), insert_indexes.end());

        // //SAVING THE INDEXES TO FILES
        // saveIndexesToFile(bulk_load_indexes, bulk_index_output);
        // saveIndexesToFile(insert_indexes, insert_index_output);

        // bulk_load_indexes = readIndexesFromFile(bulk_index_output);
        // insert_indexes = readIndexesFromFile(insert_index_output);
        std::string bulk_index_output = data_folder + "Splits/"+  dataset_name +"_bulk.bin";
        bulk_load_indexes = readIndexesFromFile(bulk_index_output);

        std::cout << "Shuffuled3." << values_size << std::endl;

        //Removing indexes that have duplicates
        // auto it = std::upper_bound(bulk_load_indexes.begin(), bulk_load_indexes.end(), values_size-1);
        // bulk_load_indexes.erase(it, bulk_load_indexes.end());

        // auto it2 = std::upper_bound(insert_indexes.begin(), insert_indexes.end(), values_size-1);
        // insert_indexes.erase(it2, insert_indexes.end());

        // bulk_load_indexes.erase(std::remove_if(bulk_load_indexes.begin(), bulk_load_indexes.end(), [values_size](int num) { return num >= values_size; }), bulk_load_indexes.end());
        // insert_indexes.erase(std::remove_if(insert_indexes.begin(), insert_indexes.end(), [values_size](int num) { return num >= values_size; }), insert_indexes.end());


        // for(int i : bulk_load_indexes ){
        //      std::cout << i  << " "; 
        // }

        std::vector<KEY_TYPE> subvector(bulk_load_indexes.size());

        //subvector.resize(bulk_load_indexes.size());
        subvector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(subvector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });

        //  std::sort(subvector.begin(), subvector.end());
        
        // for(int i : subvector ){
        //      std::cout << i  << " "; 
        // }
       

        //SAVING THE DATA TO FILES
        // save_data_bin(subvector, bulk_data_output);
        
        // subvector.clear();
        
        // subvector = read_data_bin(bulk_data_output);

        // std::sort(subvector.begin(), subvector.end());

        // if (check_duplicates(subvector)) {
        // std::cout << "The original data has duplicate values." << std::endl;
        // } else {
        //     std::cout << "The original data does not have duplicate values." << std::endl;
        // }


        // for(int i : subvector ){
        //      std::cout << i  << " "; 
        // }
        std::cout << "bulk_load_indexes "  << bulk_load_indexes.size() << std::endl;
        std::cout << "insert_indexes "  << insert_indexes.size() << std::endl;

        std::cout << "subvector "  << subvector.size() << std::endl;

        // std::cout << "Right"  << std::endl;
        numberOfIndexes = subvector.size();
        auto values = create_values(subvector,&numberOfIndexes);
        std::cout << "Created values" << std::endl;
        //std::pair<KEY_TYPE, PAYLOAD_TYPE> *values2 = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[numberOfIndexes];

        // std::copy_if(values, values + values_size, values2,
        //          [&](const std::pair<uint64_t, double>& value) {
        //              return std::find(bulk_load_indexes.begin(), bulk_load_indexes.end(),
        //                               static_cast<KEY_TYPE>(value.first)) != bulk_load_indexes.end();
        //          });

        //index_original.bulk_load(values, numberOfIndexes);
        
        index.bulk_load(values, numberOfIndexes);
        // sali.bulk_load(values, numberOfIndexes);
        // index_original.bulk_load(values, values_size);
        // index.bulk_load(values, values_size);

    }
    else{
        std::cout << "Creating values" << std::endl;
        auto values = create_values(legitimate_data,&values_size);

        std::cout << "Created values" << std::endl;
        //index_original.bulk_load(values, values_size);
         //std::cout << "Created index1" << std::endl;
        index.bulk_load(values, values_size);
        //sali.bulk_load(values, values_size);
        delete[] values; // Delete the dynamically allocated array
        values = nullptr; // Reset the pointer to nullptr

        std::cout << "Created index" << std::endl;
    }
    auto stop_build = std::chrono::high_resolution_clock::now();

    std::vector<int>().swap(all_indexes);
    // std::vector<int>().swap(bulk_load_indexes);

   
    std::cout << "============================" << std::endl;
    
    // for(int i =0 ; i < values_size; i ++ ){
    //     std::cout << " " << values[i].first << std::endl;
    // }
    //bulk_load indexes and measure their times per index build
    // auto start_build = std::chrono::high_resolution_clock::now();

    // index_original.bulk_load(values, values_size);
    // index.bulk_load(values, values_size);
    

    // auto stop_build = std::chrono::high_resolution_clock::now();

    std::cout << "Created the index" << std::endl;
    std::cout << "============================" << std::endl;

    long index_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_build - start_build).count();
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;

    

    //Calculating the new constants by using the created original index
    auto start_tra = std::chrono::high_resolution_clock::now();

    int max_model_height =0;

    std::vector<Node*> nodes;
   
    //Need to change this scan_nodes
    index.scan_nodes(&nodes);

    //Testing area
    //std::vector<Node_S*> nodes_sali;
    //sali.scan_nodes(&nodes_sali);

    //int size = 50000;
            

    // KEY_TYPE* keys = new KEY_TYPE[size];
    // PAYLOAD_TYPE* values = new PAYLOAD_TYPE[size];

    // sali.scan_subtree2(sali.root, keys, values,false);

    std::sort(nodes.begin(), nodes.end(), compare_level_des);
    max_model_height = nodes[0]->level;
    

    //num_pois = 0;
    //num_pois = max_model_height;

    int replaced_models_number = 0;

    //std::vector<model_node_type*> *model_nodes_by_level = new std::vector<model_node_type*>[max_model_height + 1];
    
    std::cout << "DONE: models useful and models by level" << std::endl;
    //std::cout << "models useful " << model_nodes_useful.size() <<std::endl;
    std::cout << "============================" << std::endl;


    int current_level = max_model_height;

    //GET USEFUL NODES BY LEVEL
    //============================
    std::vector<Node*> nodes_by_level;
    //get_useful_model_nodes_by_level(&model_nodes_useful, &model_nodes_useful_by_level, current_level);
    current_level = 2;
    //Change this
    //Changed
    get_nodes_by_level_with_child( &nodes, &nodes_by_level, current_level);
    //get_nodes_by_level(&nodes, &nodes_by_level, current_level);

    
    std::cout << "Current_level " << current_level << std::endl;
    std::cout << "nodes_by_level " << nodes_by_level.size() << std::endl;

    //POISONING POINTS DISTRIBUTION
    // Get the nodes by level with children (otherwise what to combine)
    // poison all nodes until max gaps reached OR threshold is reached OR gradual 4 is reached
    // Save the derivative of the poisoning part and the derivative of the gradual part. To find the gradual and non gradual use thresholds from the points
    // Calculate the cost for each
    // Divide them by the size.
    // Assign points by sorting them by cost. Readjust the available points too 
    double *dervs = new double[nodes_by_level.size()];

    double *grad_dervs = new double[nodes_by_level.size()];
    
    double *maxs = new double[nodes_by_level.size()];

    double *elbows = new double[nodes_by_level.size()];

    double *costs = new double[nodes_by_level.size()];

    double *grad_costs = new double[nodes_by_level.size()];

    double *propotions = new double[nodes_by_level.size()];

    int *counts = new int[nodes_by_level.size()];

    double sum_costs = 0 ;
    double sum_costs_completed = 0 ;

    int cost_0 = 0 ;
    int grad_cost_0 = 0 ;
    int propotions_0 = 0 ;
    int total_size = 0;

    std::vector<int> indices(nodes_by_level.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    // if(poison_distribute){
        
    //     double derv = 0;
    //     double grad_derv = 0;

    //     int max = 0;
    //     int elbow = 0;

    //     int count_node = 0;

    //     for(int i =0 ; i < nodes_by_level.size();i++){
    //         count_node = 0;
    //         count_data_node(nodes_by_level[i],&count_node);
    //         // counts[i] = count_node;
    //         counts[i] = nodes_by_level[i]->size - count_node;

    //     }

    //     // for(int i =0 ; i < nodes_by_level.size();i++){
    //     //     total_size += nodes_by_level[i]->size;
    //     // }

    //     for(int i =0 ; i < nodes_by_level.size();i++){
    //         int size = nodes_by_level[i]->size;
    //         total_size += size;

    //         KEY_TYPE* keys = new KEY_TYPE[size];
    //         PAYLOAD_TYPE* values = new PAYLOAD_TYPE[size];

    //         //MAYBE CHANGE POISONED TO WORK WITH POINTERS TOO
    //         // index.scan_subtree(best_node, &model_node_data, &model_node_payload,false);
    //         index.scan_subtree2(nodes_by_level[i], keys, values,false);
    //         //index.scan_and_destory_tree(best_node, keys, values,false);

    //         std::vector<KEY_TYPE>model_node_data(keys, keys + size);
    //         KEY_TYPE first_val = model_node_data[0];
    //         shift_vector(&model_node_data, first_val);

    //         //std::cout<< "Size " << size << std::endl;
            
    //         perform_poisoning2(model_node_data, 0.5, &derv, &grad_derv, &elbow, &max);
    //         //perform_poisoning(model_node_data, 20);
    //         dervs[i] = derv;
    //         maxs[i] = max;

    //         grad_dervs[i] = grad_derv;
    //         elbows[i] = elbow;
    //     }
        

    //     normalize_array(dervs, nodes_by_level.size());
    //     normalize_array(grad_dervs, nodes_by_level.size());
    //     normalize_array(costs, nodes_by_level.size());
        
    //     for(int i =0 ; i < nodes_by_level.size();i++){
    //         //std::cout<< "Derv "<< dervs[i] <<std::endl;
    //         //std::cout<< "Max "<< maxs[i] <<std::endl;
    //         // costs[i] = elbows[i]*dervs[i];

    //         // if(grad_dervs[i] <= 0){
    //         //     costs[i] = (nodes_by_level[i]->size - counts[i])*dervs[i];
    //         // }
            
    //         // grad_costs[i] = (nodes_by_level[i]->size -counts[i]- elbows[i])*grad_dervs[i];

    //         // costs[i] = pow(nodes_by_level[i]->size - counts[i],pow_val)+constant_cost*dervs[i];
    //         // grad_costs[i] = pow(nodes_by_level[i]->size - counts[i],pow_val)+constant_cost*grad_dervs[i];
    //         costs[i] = pow(counts[i],pow_val)*constant_cost*dervs[i];
    //         grad_costs[i] = pow(counts[i],pow_val)*constant_cost*grad_dervs[i];

    //         // costs[i] = dervs[i];
    //         //  grad_costs[i] = grad_dervs[i];

    //         sum_costs += (costs[i] + grad_costs[i]);

    //         // std::cout<< "Derv "<< dervs[i] <<std::endl;
    //         // std::cout<< "Max "<< maxs[i] <<std::endl;

    //         // std::cout<< "Grad Derv "<< grad_dervs[i] <<std::endl;
    //         // std::cout<< "elbows "<< elbows[i] <<std::endl;
    //     }
    //     int P_i = 0;

    //     for(int i =0 ; i < nodes_by_level.size();i++){
    //         costs[i] = costs[i]/sum_costs; 

    //          grad_costs[i] = grad_costs[i]/sum_costs;
    //         //std::cout<< "Costs "<< costs[i] <<std::endl;
    //         if(costs[i] <= 0 ){
    //             cost_0++;
    //         }
    //         if(grad_costs[i] <= 0 ){
    //             grad_cost_0++;
    //         }
            
    //         P_i = std::ceil(costs[i]*P);

    //         if(P_i > elbows[i] && grad_costs[i] > 0){
    //             P_i = elbows[i];
    //         }
    //         P_i += std::ceil(grad_costs[i]*P);

    //          propotions[i] = P_i/static_cast<double>(nodes_by_level[i]->size);
    //         //propotions[i] = P_i;

    //         if(propotions[i] <= 0 ){
    //             propotions_0++;
    //         }
  
    //     }

    //     // std::sort(indices.begin(), indices.end(), [&](int a, int b) {
    //     // return compare_descending_index(a, b, propotions);
    //     //  });

    //       std::sort(indices.begin(), indices.end(), [&costs, &grad_costs](size_t i, size_t j) {
    //     // Compare the combined values (values1[i], values2[i]) and (values1[j], values2[j]) in descending order
    //     return std::max(costs[i], grad_costs[i]) < std::max(costs[j], grad_costs[j]);
    //     });

            

    //     // std::sort(indices.begin(), indices.end(), [&](int a, int b) {
    //     // return compare_descending_index(a, b, costs);
    //     //  });

    //     //  for(int i =0 ; i < indices.size();i++){ 
    //     //     std::cout<< "indices "<< indices[i] << " Cost " << dervs[indices[i]] << " Cost " << (counts[indices[i]])<< " propotion " << propotions[indices[i]] << std::endl;
    //     // }

    // }
    //NEXT sort according to the costs. Highest first. Get the new poisoning number and poison. If there is extra.
    // Recalculate the new P accordingly.
    //and go for the next one.
    
    //Do for all the levels except the root (keep the root as it is)
    //Inner index to iterate through the models in the current level (reset in the outter while loop)
    int inner_idx = 0;
    //int inner_idx = propotions_0 - 0;
    //To keep count of the number of models we have tried without use
    int tried = 0;
    int sucess =  0 ;

    std::cout << "Counts 0 " << cost_0 << std::endl;
    std::cout << "Counts 0 " << grad_cost_0 << std::endl;
    std::cout << "Counts 0 " << propotions_0 << std::endl;

    int current_count = 0;
    int over_count = 0;

    //TESTING 
    //========
    // Node* test_node = nodes_by_level[indices[0]];
    // int node_height=0;
    // //Change
    // get_data_node_max_height(test_node, &node_height);
    // std::vector<std::vector<KEY_TYPE>*> node_data;
    // //std::vector<KEY_TYPE> * node_data[node_height];
    // // std::cout << "Done1 " <<  std::endl;
    // // std::cout << "Level " << node_height << std::endl;
    // node_data.resize(node_height+1);
    
    // for(int i = 0; i <= node_height; i++){
    //     node_data[i] = new std::vector<KEY_TYPE>();
    // }
    // std::cout << "Level "  << std::endl;
    // node_data[0]->push_back(2);
    // std::cout << "Done2 " <<  std::endl;
    
    // //Change
    // get_data_node_data_by_level(test_node, node_data);
    // // std::cout << "Done2 " <<  std::endl;
    // //  std::cout << "Size " << test_node->size << std::endl;

    // for(int i = 0; i <= node_height; i++ ){
    //     std::cout << "Level " << i << std::endl;
    //     // for(int j = 0 ; j < node_data[i]->size(); j++){
    //     //     // if(node_data[i]->size()==0){
    //     //     //     break;
    //     //     // }
    //     //     std::cout << "Size " << node_data[i]->size() << std::endl;
    //     //     std::cout << (*node_data[i])[j] << " , ";
    //     // }
    //     for (KEY_TYPE data : *(node_data[i])) {
    //         std::cout << data << " ";
    //     }
    // }

    
    //  std::cout << "Done3 " <<  std::endl;

    std::vector<KEY_TYPE>  altered_data;
    int size_to_check;

    if(nodes_by_level.size() <= 0){
        current_level = 0;
    }
    
    while(current_level == 2 && poison){
    //while(current_level > 0){
        std::cout << "while outter " << current_level << std::endl;
        

        Node* best_node;
        Node* new_node ;
        //std::sort(nodes_by_level.begin(), nodes_by_level.end(), compare_models_new);
        bool model_converted = false;
        inner_idx = 0;
    
        size_to_check = nodes_by_level.size();
        int original_size_to_check = size_to_check;

        tried = 0;
        sucess =  0; 

       // inner_idx = 6;
        
         while(inner_idx < size_to_check && tried < max_check  ){

            if(inner_idx%10000 ==0){
                std::cout << "starting " << inner_idx<<std::endl;
            }
            //std::cout << "starting " << inner_idx<<std::endl;

            int node_height_original_=0;
            std::vector<std::vector<KEY_TYPE>*> node_data_original;
        //while(inner_idx < 10 && tried < max_check  ){

            // if(sucess > 2000){
            //     break;
            // }

            // best_node = nodes_by_level[inner_idx];
            //best_node = nodes_by_level[indices[inner_idx]];
            best_node = nodes_by_level[inner_idx];
            double best_children_cost = 0;

            //std::cout << "Got node " << inner_idx<<std::endl;

            // std::vector<KEY_TYPE> model_node_data;
            // std::vector<PAYLOAD_TYPE> model_node_payload;

            //GET ALL DATA FROM SUBTREE
            //===========================
            int size = best_node->size;
            KEY_TYPE* keys = new KEY_TYPE[size];
            PAYLOAD_TYPE* values = new PAYLOAD_TYPE[size];

            //std::cout << "Setup arrays " << inner_idx<<std::endl;
            //MAYBE CHANGE POISONED TO WORK WITH POINTERS TOO
            // index.scan_subtree(best_node, &model_node_data, &model_node_payload,false);
            index.scan_subtree2(best_node, keys, values,false);
            //index.scan_and_destory_tree(best_node, keys, values,false);
            //std::cout << "Scanned " << inner_idx<<std::endl;

            std::vector<KEY_TYPE>model_node_data(keys, keys + size);
            KEY_TYPE first_val = model_node_data[0];
            shift_vector(&model_node_data, first_val);

            std::vector<PAYLOAD_TYPE>model_node_payload(values, values + size);

            //std::sort(model_node_data.begin(),model_node_data.end());
            std::vector<KEY_TYPE>  level_data_node;
            std::vector<PAYLOAD_TYPE>  level_payload_node;

            //Change
            scan_node(best_node, &level_data_node, &level_payload_node);
            //std::cout << "Scanned node" << inner_idx<<std::endl;
            //if(level_data_node.size()== 0){
            if(level_data_node.size()== 0 || model_node_data.size()>1000){
                // if(level_data_node.size()== 0){
                //     get_nodes_children(best_node, &nodes_by_level);
                // }
                //Changed
                 get_nodes_children(best_node, &nodes_by_level);
                 over_count++;
                 size_to_check = nodes_by_level.size();
                 inner_idx++;
                 tried++;
                continue;
            }

            //get_children_data(index, best_node, &model_node_data,&model_node_payload, &best_children_cost, const_search, const_traversal, false);

            // std::cout << "DONE: model & children node data" << std::endl;
            // std::cout << "============================" << std::endl;


            //Poison using the model data or till the cost is less than the children nodes sum of costs.
            //Vector to hold the poisoned model data or if no poisoning use this as well. I put these outside the if becasue I need to use it in else as well
            std::vector<KEY_TYPE> poisoned_model_data;
            //poi_thres = (pois_count_model/(1.0*model_node_data.size()));
            // poi_thres = original_poi_thres;

            //FOR DISTRIBUTION
            //===============
            // KEY_TYPE P_i = static_cast<KEY_TYPE>(std::round(costs[indices[inner_idx]]*P));

            // if(P_i > elbows[indices[inner_idx]] && grad_costs[indices[inner_idx]] > 0){
            //     P_i = static_cast<KEY_TYPE>(elbows[indices[inner_idx]]);
            // }
            // P_i += static_cast<KEY_TYPE>(std::round(grad_costs[indices[inner_idx]]*P));

            // int current = inner_idx;

            // // if(P_i < 1 && costs[indices[inner_idx]] > 0 || P_i < 1 && grad_costs[indices[inner_idx]] > 0 ){
            // //     // P_i = P/(1.0*(size_to_check - current));
            // //     P_i = 1;
            // //     //P_i = std::ceil(static_cast<long double>(P_left)/(1.0*(size_to_check - current)));
            // // }
            
            // if(P_i < static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current))))){
            //     // if(current > propotions_0-1){
            //     //     P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current))));
            //     // }
            //     // if(current < size_to_check - propotions_0-1){
            //     //     P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check - propotions_0 - current))));
            //     // }
            //      if(propotions[indices[inner_idx]] > 0){
            //          P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check   - current))));
            //         // P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - propotions_0 - current_count))));
            //         // current_count++;
            //     }
            // }
            // if(P_i < std::ceil(static_cast<long double>(P_left)/(1.0*(576 - current)))){
            //     P_i = std::ceil(static_cast<long double>(P_left)/(1.0*(576 - current)));
            // }
            // if(P_i < model_node_data.size()*original_poi_thres){
            //     P_i =  model_node_data.size()*original_poi_thres;
            // }
            // if(current>propotions_0 - 1){
            //    // P_i = std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current)));
            //     P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current))));
            // }
            // if(propotions[indices[inner_idx]] > 0){
                    
            //         P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - propotions_0 - current_count))));
            //         current_count++;
            //     }
            //==============

            //FOR NO POISONING DISTRIBUTION AKA EQUALLY
            KEY_TYPE P_i;
            if(inner_idx < original_size_to_check){
                P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(original_size_to_check  - current_count))));
            }
            else{
                //P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check - original_size_to_check - current_count))));
                P_i = std::max(static_cast<KEY_TYPE>(model_node_data.size()*0.1),
                static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current_count))))) ;
            }
            //KEY_TYPE P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check  - current_count))));
            current_count++;
            
            // if(current < size_to_check - propotions_0-1){
            //         P_i = static_cast<KEY_TYPE>(std::round(static_cast<long double>(P_left)/(1.0*(size_to_check - propotions_0 - current))));
            //     }
            
            

             poi_thres = P_i/(1.0*model_node_data.size());
            //poi_thres = 0.1;
            //poi_thres = original_poi_thres*(max_model_height-current_level+1);
            int cutoff = -100 ;
            int poisoned_model_data_size = model_node_data.size();
            data_node_type* new_data_node;
            int current_poi =0;

            sucess++;
            //poi_thres = original_poi_thres;
            //If the cost differebce is greater than 0 meaning model cost > children cost try poisoning
            if(P_i > 0){
            //if(true){
                
                
                // std::cout << "while inner POISONING " << inner_idx << std::endl;
                // std::cout << "Size to check " << size_to_check << std::endl;
                // std::cout << "Keys : "<< model_node_data.size()  <<std::endl;
                // std::cout << "P : "<< P  <<std::endl;
                // std::cout << "P left : "<< P_left  <<std::endl;
                // std::cout << "Pi : "<< P_i  <<std::endl;
                //std::cout << "poi_thres : "<< poi_thres  <<std::endl;

                // std::cout << "cost : "<< costs[indices[inner_idx]]  <<std::endl;
                // std::cout << "grad : "<< grad_costs[indices[inner_idx]]  <<std::endl;
                // std::cout << "elbow : "<< elbows[indices[inner_idx]]  <<std::endl;

                // for(int i = 0 ; i < model_node_data.size(); i++){
                //     std::cout << model_node_data[i] << " , ";
                // }
                
                poisoned_model_data = perform_poisoning(model_node_data, poi_thres);
                shift_back_vector(&poisoned_model_data, first_val);
                shift_back_vector(&model_node_data, first_val);

                poisoned_model_data_size = poisoned_model_data.size();

                KEY_TYPE* keys_poi = new KEY_TYPE[poisoned_model_data_size];
                PAYLOAD_TYPE* values_poi = new PAYLOAD_TYPE[poisoned_model_data_size];

                //get new values for the poisoning data while using the old ones for the existing
                //auto poisoned_values = get_poisoned_values(model_node_data, poisoned_model_data, model_node_payload, poisoned_model_data_size);
                get_poisoned_keys(model_node_data, poisoned_model_data, model_node_payload, keys_poi, values_poi,poisoned_model_data_size);
                
                //std::cout << "Keys : "<< poisoned_model_data_size  <<std::endl;
               //ALso check the max data node number here itself
               //Add a constant that if there were no new points also save
                //  if(index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size ){
                //     inner_idx++;
                //     tried++;
                //     //replace_parents_data_node(new_data_node,model_nodes_by_level, &best_node_idx, &parent_model_idx);
                //     //std::cout << "SKIP" << std::endl;
                //     continue;
                // }
                // std::cout << "Creating new data" << std::endl;
                

                //CREATE THE NEW FMCD
                //=====================
                // for(int i = 0 ; i < model_node_data.size(); i++){
                //     std::cout << model_node_data[i] << " , ";
                // }
                //  std::cout  << std::endl;

                // for(int i = 0 ; i < poisoned_model_data_size; i++){
                //     std::cout << keys_poi[i] << " , ";
                // }

                //I need to find the node of the best node and set it to the new one. Not the root.
                //std::cout << level_data_node[0] << std::endl;
                // std::cout << "Model a : "<< best_node->model.a  <<std::endl;
                // std::cout << "Model b : "<< best_node->model.b  <<std::endl;
                // double c =  elbows[indices[inner_idx]]*dervs[indices[inner_idx]];

                // if(grad_dervs[indices[inner_idx]] <= 0){
                //     c = (best_node->size - counts[indices[inner_idx]])*dervs[indices[inner_idx]];
                // }
                // double g = (best_node->size -  counts[indices[inner_idx]]- elbows[indices[inner_idx]])*grad_dervs[indices[inner_idx]];

                // double c =  pow(best_node->size -  counts[indices[inner_idx]],pow_val)+constant_cost*dervs[indices[inner_idx]];
                // double g = pow(best_node->size -  counts[indices[inner_idx]],pow_val)+constant_cost*grad_dervs[indices[inner_idx]];

                // double c =  pow(counts[indices[inner_idx]],pow_val)*constant_cost*dervs[indices[inner_idx]];
                // double g = pow(counts[indices[inner_idx]],pow_val)*constant_cost*grad_dervs[indices[inner_idx]];

                // double c =  dervs[indices[inner_idx]];
                // double g = grad_dervs[indices[inner_idx]];
                
                //sum_costs_completed += (c + g) ;
                
                P_left -= (poisoned_model_data_size -  model_node_data.size());
                //P = static_cast<KEY_TYPE>((P_left)*sum_costs/(sum_costs - sum_costs_completed));

                //grad_costs[i] = (nodes_by_level[i]->size - elbows[i])*grad_dervs[i]/(static_cast<double>(total_size));

                //sum_costs_completed += (dervs[indices[inner_idx]]+ grad_dervs[indices[inner_idx]]) ;
                // sum_costs_completed += (best_node->size*dervs[indices[inner_idx]]/(static_cast<double>(total_size)) + (best_node->size)*grad_dervs[indices[inner_idx]]/(static_cast<double>(total_size)));
                // sum_costs_completed += (best_node->size*dervs[indices[inner_idx]] + best_node->size*grad_dervs[indices[inner_idx]]) ;
                
                //  P = static_cast<KEY_TYPE>((P_left)*sum_costs/(sum_costs - sum_costs_completed));
                //std::cout << "P_left : "<< P_left  <<std::endl;

                //For the original node, keep the level data.
                //Return the new node instead of the root. 
                //get the data one level by one and then check the difference between them and add them to the changed data 
                //std::vector<KEY_TYPE> * node_data_original;
                

                if(poisoned_model_data_size > model_node_data.size()){
                    replaced_models_number++;

                    //Getting the original data in each level of it
                    
                    get_data_node_max_height(best_node, &node_height_original_);
                    //std::vector<KEY_TYPE> * node_data_original[node_height_original_];
                    // std::cout << "Done1 " <<  std::endl;
                   // std::cout << "Height " << node_height_original_ << std::endl;

                    node_data_original.resize(node_height_original_+1);

                    for(int i = 0; i <= node_height_original_; i++){
                        node_data_original[i] = new std::vector<KEY_TYPE>();
                    }
                    //node_data[0]->push_back(2);
                    //std::cout << "Done2 " <<  std::endl;
                    
                    get_data_node_data_by_level(best_node, node_data_original);
                    // std::cout << "Done2 " <<  std::endl;
                    //  std::cout << "Size " << test_node->size << std::endl;

                    // for(int i = 0; i <= node_height_original_; i++ ){
                    //     std::cout << "Level " << i << std::endl;
                    //     // for(int j = 0 ; j < node_data[i]->size(); j++){
                    //     //     // if(node_data[i]->size()==0){
                    //     //     //     break;
                    //     //     // }
                    //     //     std::cout << "Size " << node_data[i]->size() << std::endl;
                    //     //     std::cout << (*node_data[i])[j] << " , ";
                    //     // }
                    //     for (KEY_TYPE data : *(node_data_original[i])) {
                    //         std::cout << data << " ";
                    //     }
                    // }
                   
                    //std::cout << "Creating new data" << std::endl;
                    //index.root = index.poison_bulk_load(index.root, best_node,level_data_node[0], level_payload_node[0], keys_poi, values_poi,poisoned_model_data_size);
                    new_node = index.poison_bulk_load_real(index.root, best_node,level_data_node[0], level_payload_node[0], keys_poi, values_poi,poisoned_model_data_size,keys, size);

                    //std::cout << "Created new data" << std::endl;
                }

                

                // std::cout << "Model a : "<< best_node->model.a  <<std::endl;
                // std::cout << "Model b : "<< best_node->model.b  <<std::endl;
                //index.insert(level_data_node[0]+100,10);

                current_poi = poisoned_model_data_size - model_node_data.size();
                //std::cout << "==========================" << inner_idx << std::endl;
            }
            else{
                //break;
            }
            //Can directly create a new data node instead of a model
            // else{
                
            //     std::cout << "while inner " << inner_idx << std::endl;

            //         auto values_real = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[model_node_data.size()];
            //         //std::mt19937_64 gen_payload(std::random_device{}());

            //             for (int i = 0; i < model_node_data.size(); i++) {
            //                 values_real[i].first = model_node_data[i];
            //                 values_real[i].second = model_node_payload[i];
            //                 }
                
            //     // if(index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= model_node_data.size()){
            //     //     inner_idx++;
            //     //     tried++;
            //     //     //replace_parents_data_node(new_data_node,model_nodes_by_level, &best_node_idx, &parent_model_idx);
            //     //     //std::cout << "SKIP" << std::endl;
            //     //     continue;
            //     // }
                
            //    //CREATE THE NEW FMCD
            //     //=====================
                

            //     poisoned_model_data_size = model_node_data.size();
            // }

            //ADD some condition to change the nodes
            int best_level = best_node->level;
            int best_node_idx = 0;
            bool parent_found = false;
            int parent_model_idx = 0;

            //If the cost of the new data is higher than the children or there is no more space then try the next best node
            //Change the condition here
            // if((new_data_node->cost_ - best_node->children_cost) > (cutoff) || index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size){

            //ADD NEW CONDITION BASED ON CONFLICTS
            //======================
             if(false){
                inner_idx++;
                tried++;
                //replace_parents_data_node(new_data_node,model_nodes_by_level, &best_node_idx, &parent_model_idx);
                //std::cout << "SKIP" << std::endl;
                continue;
            }

            //Get the new node. 

            //if(current_level == max_model_height){
                //  std::cout << "Adding"  <<std::endl;
                if(poisoned_model_data_size > model_node_data.size()){

                    //Getting new data after changing in each level.
                    int node_height_changed = 0;
                    get_data_node_max_height(new_node, &node_height_changed);
                    //std::vector<KEY_TYPE> * node_data_changed[node_height_changed];
                    std::vector<std::vector<KEY_TYPE>*> node_data_changed;
                    // std::cout << "Done1 " <<  std::endl;
                    
                    //std::cout << "Height " << node_height_changed << std::endl;

                    node_data_changed.resize(node_height_changed+1);
                    
                    

                    for(int i = 0; i <= node_height_changed; i++){
                        node_data_changed[i] = new std::vector<KEY_TYPE>();
                    }
                    
                    //node_data[0]->push_back(2);
                    //std::cout << "Done2 " <<  std::endl;
                    get_data_node_data_by_level(new_node, node_data_changed);
                    // std::cout << "Done2 " <<  std::endl;
                    //  std::cout << "Size " << test_node->size << std::endl;

                    

                    // for(int i = 0; i <= node_height_changed; i++ ){
                    //     std::cout << "Level " << i << std::endl;
                    //     // for(int j = 0 ; j < node_data[i]->size(); j++){
                    //     //     // if(node_data[i]->size()==0){
                    //     //     //     break;
                    //     //     // }
                    //     //     std::cout << "Size " << node_data[i]->size() << std::endl;
                    //     //     std::cout << (*node_data[i])[j] << " , ";
                    //     // }
                    //     for (KEY_TYPE data : *(node_data_changed[i])) {
                    //         std::cout << data << " ";
                    //     }
                    // }
                    //Getting the differences
                        // std::cout << "Getting difference " << std::endl;
                        // std::cout << "inner  " << inner_idx << std::endl;
                        // std::cout << node_height_original_<< " "<< node_height_changed<< std::endl;

                    int start = 1;
                    if(inner_idx >= original_size_to_check){
                        start = 0;
                    }
                    for(int i = start; i <= node_height_original_; i++ ){
                        // std::cout << "Level " << i << std::endl;
                        // for(int j = 0 ; j < node_data[i]->size(); j++){
                        //     // if(node_data[i]->size()==0){
                        //     //     break;
                        //     // }
                        //     std::cout << "Size " << node_data[i]->size() << std::endl;
                        //     std::cout << (*node_data[i])[j] << " , ";
                        // }
                        
                        //If the node is higher than the current
                        if(i > node_height_changed){
                            // std::cout << "Inside " << i << std::endl;
                            (altered_data).insert((altered_data).end(), (*node_data_original[i]).begin(), (*node_data_original[i]).end());
                        }
                        else{
                            find_difference(&altered_data,*node_data_original[i], *node_data_changed[i]);
                        }
                        // for (KEY_TYPE data : *(node_data_changed[i])) {
                        //     std::cout << data << " ";
                        // }
                        // std::cout << "Done " << i << std::endl;
                        // for (KEY_TYPE data : *(node_data_original[i])) {
                        //     std::cout << data << " ";
                        // }
                    }



                    for (const auto& element : model_node_data) {
                // Check if the element is not in the same level node (meaning it is from the bottom)
                        if (std::find(level_data_node.begin(), level_data_node.end(), element) == level_data_node.end()) {
                            // If not found, add it to the original vector
                            old_data.push_back(element);
                        }
                
                    }
                }
                // else{
                //     //This wont work if I do the reconstruction because now everything is new. because now the best node does not exist.
                //    // std::cout << "Here" <<std::endl;
                //     get_nodes_children(best_node, &nodes_by_level);
                //     size_to_check = nodes_by_level.size();
                //     //std::cout << "Here Done" <<std::endl;
                // }
                
            // }
                
            
            
            tried = 0;

            std::vector<model_node_type*> parents;
           
            new_data_node = static_cast<data_node_type*>(new_data_node);
            
            //REPLACE THE DATA STRUCTURES
            //=====================


            // parents = replace_parents(best_node, new_data_node, model_nodes_by_level, &best_node_idx, &parent_model_idx);
            
            // update_data_structure(best_node, &model_nodes, model_nodes_by_level, &model_nodes_useful,&model_nodes_useful_by_level, parents,
            // index,const_search,const_traversal,&new_useful_parents);


            

            model_converted = true;

            
            num_total_poisoning += current_poi;
            //If model is converted then set tried back to 0 and increase the inner_idx
            if(model_converted){
                inner_idx++;
                tried = 0;
            }
            // node_height=0;
            // get_data_node_max_height(new_node, &node_height);
            // std::vector<KEY_TYPE> * node_data2[node_height];
            // // std::cout << "Done1 " <<  std::endl;
            // std::cout << "Height " << node_height << std::endl;

            // for(int i = 0; i <= node_height; i++){
            //     node_data2[i] = new std::vector<KEY_TYPE>();
            // }
            // //node_data[0]->push_back(2);
            // //std::cout << "Done2 " <<  std::endl;
            // get_data_node_data_by_level(new_node, node_data2);
            // // std::cout << "Done2 " <<  std::endl;
            // //  std::cout << "Size " << test_node->size << std::endl;

            // for(int i = 0; i <= node_height; i++ ){
            //     std::cout << "Level " << i << std::endl;
            //     // for(int j = 0 ; j < node_data[i]->size(); j++){
            //     //     // if(node_data[i]->size()==0){
            //     //     //     break;
            //     //     // }
            //     //     std::cout << "Size " << node_data[i]->size() << std::endl;
            //     //     std::cout << (*node_data[i])[j] << " , ";
            //     // }
            //     for (KEY_TYPE data : *(node_data2[i])) {
            //         std::cout << data << " ";
            //     }
            // }
            
            // std::cout << over_count <<std::endl;
            //  std::cout << nodes_by_level.size()  <<std::endl;
        }
        
        //After a level go to the upper level
        current_level--;
        if(current_level  > 1){
            //Get the new models by level from model nodes useful (which we update in the update_data_structures
            //get_useful_model_nodes_by_level(&model_nodes_useful, &model_nodes_useful_by_level, current_level);

            //GET USEFULL NODES BY LEVEL
            //get_nodes_by_level_with_child( &nodes, &nodes_by_level, current_level);

        }
        
        std::cout << "=============="  <<std::endl;
        std::cout << "Done: Iteration " << current_level <<std::endl;
        std::cout << "Prop 0 " << propotions_0 <<std::endl;
            
    }
    
    // std::sort(changed_data.begin(), changed_data.end());
    // std::vector<node_type*> nodes_fin;
    // std::vector<model_node_type*> model_nodes_fin;
    // std::vector<data_node_type*> data_nodes_fin;

    //Getting all the nodes from the index and the maximum height (default is 0 aka root node only)
    //Printing information at the end
    int max_model_height_fin = 0;
    
    //GET NEW NODE STATS 
                //=====================
    //get_all_nodes(index, &nodes_fin, &model_nodes_fin,&data_nodes_fin, &max_model_height_fin);


    auto stop_tra = std::chrono::high_resolution_clock::now();

    //Original information
    // std::cout << "====================" << std::endl;

    // std::cout << "original nodes number "<< original_nodes_number << std::endl;
    // std::cout << "original model nodes number "<< original_model_nodes_number << std::endl;
    // std::cout << "original data nodes number "<< original_data_nodes_number << std::endl;

    // std::cout << "====================" << std::endl;

    // //afterwards
    // std::cout << "nodes number "<< nodes_fin.size() << std::endl;
    // std::cout << "model nodes number "<< model_nodes_fin.size() << std::endl;
    // std::cout << "data nodes number "<< data_nodes_fin.size() << std::endl;

    // std::cout << "====================" << std::endl;

    // std::cout << "model_nodes : "<< model_nodes.size() << std::endl;
    // std::cout << "model_nodes_useful " << model_nodes_useful.size() << std::endl;
    std::cout << " Number " << index.root->size <<std::endl;
   // std::cout << " Number " << index_original.root->size <<std::endl;
     std::cout << "Replaced number " << replaced_models_number <<std::endl;

    
    std::cout << "Poisoned number " << num_total_poisoning <<std::endl;
    std::cout << "Size to check " << size_to_check << std::endl;
    long tra_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra - start_tra).count();

    //TO POISON THE DATA NODES
    //==================
    // for(int i = max_model_height+1; i > max_model_height  ; i-- ){
    //     std::vector<data_node_type*> data_nodes_level;
    //     get_data_nodes_by_level(&data_nodes_fin,&data_nodes_level,i);
    //     poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);
    // }
    //poison_data_nodes(index,&data_nodes_fin, poi_thres, const_search, const_traversal,&num_total_poisoning,&changed_data,0.90);
    

    // std::vector<data_node_type*> data_nodes_level;
    // get_data_nodes_by_level(&data_nodes_fin,&data_nodes_level,max_model_height-1);
    
    // std::vector<KEY_TYPE> level_data;
    // get_data_by_level(&data_nodes_level, &level_data);

    // benchmark_alex_real_seperate(index_original, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
    //Out putting the performance by level

    //Now insert the rest
    if(poison){
        save_data_bin(altered_data,changed_data_output);
    }
    // std::vector<KEY_TYPE> new_altered = read_data_bin(changed_data_output);

    if(!poison&benchmark){
        altered_data.clear();
        altered_data = read_data_bin(changed_data_output);
    }

    std::vector<KEY_TYPE> bulk_vector(bulk_load_indexes.size());

        //subvector.resize(bulk_load_indexes.size());
        bulk_vector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(bulk_vector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });
    int bulk_size = bulk_load_indexes.size();
    std::vector<KEY_TYPE> lookup_keys = get_search_keys(bulk_vector, bulk_size, 1000000, seed);
    // save_data_bin(lookup_keys,lookups_output);
    // lookup_keys.clear();
    // lookup_keys = read_data_bin(lookups_output);

    std::vector<KEY_TYPE> lookup_keys_zipf = get_search_keys_zipf(bulk_vector, bulk_size, 1000000, seed);

    benchmark_sali_real(index,altered_data, altered_data, "SALI", dataset_name+"_changed", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);

    benchmark_sali_real(index,lookup_keys, lookup_keys, "SALI", dataset_name+"_lookup", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index,lookup_keys_zipf, lookup_keys_zipf, "SALI", dataset_name+ "_lookup_zipf", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);

    //Benchmark normal
    std::vector<Node*> nodes2;
    index.scan_nodes(&nodes2);
    std::sort(nodes2.begin(), nodes2.end(), compare_level_des);
    max_model_height = nodes2[0]->level;

    //index_original.scan_and_destory_tree(nodes[0], keys, values,false);
    size_t total_node_count = 0;
    size_t total_data_count = 0;

    std::string structure_results = "SALI;" + dataset_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";Data";

        // std::cout << test_output <<std::endl; 
    
    std::string node_info = "Nodes";
    std::vector<KEY_TYPE> lookup_keys_level_prop_lower;
    std::vector<KEY_TYPE> lookup_keys_level_prop;

    int level_one_size = 0;
    int level_two_size = 0;

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

        //Saving the node info into the 
        node_info = node_info + ";" + std::to_string(count);
        //benchmark_alex_real_seperate(index_original, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
        //benchmark_lipp_real(index, level_data, level_data, "LIPP", data_name+" "+std::to_string(i), poi_thres, original_changed_output);
        //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);

    //     if (check_duplicates(level_data)) {
    //     std::cout << "The original data has duplicate values." << std::endl;
    // } else {
    //     std::cout << "The original data does not have duplicate values." << std::endl;
    // }
        //Selecting and benchmarking the level information.
        std::vector<KEY_TYPE> lookup_keys_level = get_search_keys(level_data, level_data.size(), 1000000, seed);

        if(i > 0){
            //int e = int(1000000/(1.0*(max_model_height-1)));
            int p = int((1000000/(1.0*values_size-level_one_size))*level_data.size());

            //std::copy(lookup_keys_level.begin(), lookup_keys_level.begin() + e, std::back_inserter(lookup_keys_level_equal));
            std::copy(lookup_keys_level.begin(), lookup_keys_level.begin() + p, std::back_inserter(lookup_keys_level_prop));
            // lookup_keys_level_equal.
            // lookup_keys_level_prop
        }
        else{
            level_one_size = level_data.size();
        }

        if(i > 1){
            int p = int((1000000/(1.0*values_size-level_one_size -level_two_size))*level_data.size());

            //std::copy(lookup_keys_level.begin(), lookup_keys_level.begin() + e, std::back_inserter(lookup_keys_level_equal));
            std::copy(lookup_keys_level.begin(), lookup_keys_level.begin() + p, std::back_inserter(lookup_keys_level_prop_lower));
        }
        if(i == 1){
            level_two_size = level_data.size();
        }

        // benchmark_sali_real(index,lookup_keys_level, lookup_keys_level, "SALI", dataset_name+"_lookup_level_"+std::to_string(i), 
        // poi_thres, performance_output,poison,insert,orignal_P,insert_threshold);
    }
    structure_results = structure_results + ";" + node_info +";" +"0";

    saveToCSV(structure_results,structure_output);

    size_t index_size_after =  index.total_size();

         std::string model_results = "SALI;" + dataset_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" +
        std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        
        saveToCSV(model_results,output_model_csv);

     if(insert){
    //if(0){      
        for(int j = 0; j < 5; j++){
            insert_threshold = (j+1)*0.1;
            std::string structure_results2 = "SALI;" + dataset_name + "_insert_"+std::to_string(j) +";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning);

            std::vector<KEY_TYPE>  inserted_data;
        
            std::string insert_index_output = data_folder+ "Splits/"+ data_name +"_insert"+"_"+std::to_string(j)+".bin";
            std::vector<int> insert_indexes2 = readIndexesFromFile(insert_index_output);


            auto start = std::chrono::high_resolution_clock::now();
            for(int i : insert_indexes2){
                //index_original.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
                index.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
                //inserted_data.push_back(legitimate_data[i]);
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count()/(1.0*insert_indexes2.size());

            std::cout << "Inserted " << std::endl;
        //     benchmark_lipp_real(index,inserted_data, inserted_data, "LIPP", dataset_name+ "_changed", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_lipp_real(index,lookup_keys, lookup_keys, "LIPP", dataset_name+"_lookup", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_lipp_real(index,lookup_keys_zipf, lookup_keys_zipf, "LIPP", dataset_name+ "_lookup_zipf", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        std::vector<Node*> nodes3;
    index.scan_nodes(&nodes3);
    std::sort(nodes3.begin(), nodes3.end(), compare_level_des);
    max_model_height = nodes3[0]->level;

    //index_original.scan_and_destory_tree(nodes[0], keys, values,false);
    size_t total_node_count2 = 0;
    size_t total_data_count2 = 0;
    std::string node_info2 = "Nodes";

    // std::string structure_results = "LIPP;" + dataset_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
    //       std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
    //     + std::to_string(num_total_poisoning);

        // std::cout << test_output <<std::endl; 
    
    
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
        //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);

    //     if (check_duplicates(level_data)) {
    //     std::cout << "The original data has duplicate values." << std::endl;
    // } else {
    //     std::cout << "The original data does not have duplicate values." << std::endl;
    // }
        // std::vector<KEY_TYPE> lookup_keys_level2 = get_search_keys(level_data2, level_data2.size(), 1000000, seed);
        // benchmark_sali_real(index,lookup_keys_level2, lookup_keys_level2, "SALI", dataset_name+"_lookup_level_"+std::to_string(i)+ "_insert_"+std::to_string(j), 
        //     poi_thres, performance_output,poison,insert,orignal_P,insert_threshold);
        
        }
        benchmark_sali_real(index,altered_data, altered_data, "SALI", dataset_name+"_changed"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
        benchmark_sali_real(index,lookup_keys, lookup_keys, "SALI", dataset_name+"_lookup"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
            poison,insert,orignal_P,insert_threshold);
        benchmark_sali_real(index,lookup_keys_zipf, lookup_keys_zipf, "SALI", dataset_name+ "_lookup_zipf"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
        poison,insert,orignal_P,insert_threshold);

        std::vector<KEY_TYPE> insert_vector(insert_indexes2.size());

        //subvector.resize(bulk_load_indexes.size());
        insert_vector.clear();

        std::transform(insert_indexes2.begin(), insert_indexes2.end(), std::back_inserter(insert_vector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });

        benchmark_sali_real(index,insert_vector, insert_vector, "SALI", dataset_name+ "_insert"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
        poison,insert,orignal_P,insert_threshold);

        
        structure_results2 = structure_results2 + ";" + node_info2 +";" +std::to_string(duration);
        saveToCSV(structure_results2,structure_output);

        size_t index_size_after2 =  index.total_size();

         std::string model_results = "SALI;" + dataset_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after2) + ";" +
        std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        
        saveToCSV(model_results,output_model_csv);
        }
        
    }

    //TO BENCHMARK LEVELS
    //index.insert(5,10);
    // std::vector<Node*> nodes2;
    // index.scan_nodes(&nodes2);
    // std::sort(nodes2.begin(), nodes2.end(), compare_level_des);
    // max_model_height = nodes2[0]->level;

    // //index_original.scan_and_destory_tree(nodes[0], keys, values,false);
    // size_t total_node_count = 0;
    // size_t total_data_count = 0;

    // std::string structure_results = "SALI;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
    //       std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
    //     + std::to_string(num_total_poisoning);

    //     // std::cout << test_output <<std::endl; 
    
    
    // for(int i = max_model_height; i >0 ; i-- ){

    //     std::cout << "Performance For Level  " << i << std::endl;

    //     std::vector<KEY_TYPE> level_data;
    //     std::vector<PAYLOAD_TYPE> level_payload;
    //     int count = 0;
    //     for(int j = 0; j < nodes2.size();j++){
    //         if(nodes2[j]->level==i){
    //             count++;
    //             scan_node(nodes2[j], &level_data, &level_payload);
    //         }
    //     }

    //     if(i > 2){
    //         total_node_count += count;
    //         total_data_count += level_data.size();
    //     }
       

    //     std::cout << "Nodes " << count <<std::endl;
    //     std::cout << "Data " << level_data.size() <<std::endl;
    //     structure_results = structure_results + ";" + std::to_string(level_data.size());

    // //     if (check_duplicates(level_data)) {
    // //     std::cout << "The original data has duplicate values." << std::endl;
    // // } else {
    // //     std::cout << "The original data does not have duplicate values." << std::endl;
    // // }
    // // for(int j = 0; j < level_data.size();j++){
    // //         std::cout << level_data[j] <<" " ;
    // //     }
    //     //benchmark_alex_real_seperate(index_original, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
    //     //benchmark_lipp_real(index, level_data, level_data, "LIPP", data_name+" "+std::to_string(i), poi_thres, original_changed_output);
    //     //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);
    // }

    //Save Loss and MSE.
    double loss_value = 0.0;
    double mse = 0.0;

    // for(Node* node : nodes2){
    //     std::vector<KEY_TYPE> x;
        
    //     std::vector<PAYLOAD_TYPE> node_payload;
    //     scan_node(node, &x, &node_payload);

    //     int n = x.size();

    //     if(n <= 1){
    //         continue;
    //     }
    //     shift_vector(&x, x[0]);

    //     std::vector<int> y(n);
    //     std::iota(y.begin(), y.end(), 0);
    //     long double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXSquare = 0.0;

    //     for(int i =0; i < n;i++){
    //         sumX += static_cast<long double>(x[i]);
    //         sumY += y[i];
    //         sumXY += static_cast<long double>(x[i]) * y[i];
    //         sumXSquare += static_cast<long double>(x[i]) * static_cast<long double>(x[i]);
    //     }

    //     double numerator = n * sumXY - sumX * sumY;
    //     double denominator = n * sumXSquare - sumX * sumX;

    //     double a = numerator / denominator;
    //     double b = (sumY - a * sumX) / n;

        
    //     for (int i = 0; i < n; ++i) {
    //         double predictedY = b + a * x[i];
    //         double error = predictedY - y[i];
    //         loss_value += error * error;
    //     }
    //     mse = loss_value/(1.0*n);

    // }
    // structure_results = structure_results + ";" + std::to_string(loss_value)+ ";" + std::to_string(mse);

    // saveToCSV(structure_results,structure_output);
    // std::vector<KEY_TYPE> new_level_data;
    // std::vector<PAYLOAD_TYPE> new_level_payload;
    // for(int j = 0; j < nodes2.size();j++){
    //         // if(nodes2[j]->level==3 || nodes2[j]->level==4 || nodes2[j]->level==5){
    //         if(nodes2[j]->level > 2){
    //             scan_node(nodes2[j], &new_level_data, &new_level_payload);
    //         }
    //     }
    // std::vector<Node*> nodes3;
    // //UNCOMMENT THIS
    // // index_original.scan_nodes(&nodes3);
    // index.scan_nodes(&nodes3);
    // std::vector<KEY_TYPE> old_level_data;
    // std::vector<PAYLOAD_TYPE> old_level_payload;
    // for(int j = 0; j < nodes3.size();j++){
    //         // if(nodes3[j]->level==3 || nodes3[j]->level==4 || nodes3[j]->level==5){
    //         if(nodes3[j]->level > 2){
    //             scan_node(nodes3[j], &old_level_data, &old_level_payload);
    //         }
    //     }
    // std::sort(old_level_data.begin(), old_level_data.end());
    // std::sort(new_level_data.begin(), new_level_data.end());

    // // for(int j = 0; j < old_level_data.size();j++){
    // //         std::cout << old_level_data[j] << " , ";
    // //     }
    // // std::cout <<  std:: endl;
    // // for(int j = 0; j < new_level_data.size();j++){
    // //         std::cout << new_level_data[j] << " , ";
    // //     }
    // std::sort(old_data.begin(), old_data.end());

    // std::cout << "Data Old " << old_data.size() <<std::endl;
    // std::cout << "Data New " << new_level_data.size() <<std::endl;
    // std::set_difference(old_data.begin(), old_data.end(), new_level_data.begin(), new_level_data.end(), std::back_inserter(changed_data));

    std::cout << "Poisoned number " << num_total_poisoning <<std::endl;


    // for(int i = max_model_height+1; i >=0 ; i--  ){
    //     std::cout << "Performance For Level  " << i << std::endl;
    //     std::vector<data_node_type*> data_nodes_level;
    //     get_data_nodes_by_level(&data_nodes_fin,&data_nodes_level,i);
    //     std::vector<KEY_TYPE> level_data;
    //     get_data_by_level(&data_nodes_level, &level_data);
    //      std::cout << "Before  " << i << std::endl;
    //     benchmark_alex_real_seperate(index_original, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
    //      std::cout << "After  " << i << std::endl;
    //     benchmark_alex_real_seperate(index, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
    //     //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);
    // }
    
    // calculate_expected_searches_all_data(index,&data_nodes_fin, poi_thres, const_search, const_traversal);
    // poison_data_nodes(index,&data_nodes_fin, poi_thres, const_search, const_traversal,&num_total_poisoning,&changed_data,0.99);

    //std::vector<data_node_type*> data_nodes_level;
    //get_data_nodes_by_level(&data_nodes_fin,&data_nodes_level,max_model_height+1);
    //std::vector<data_node_type*> first_10(data_nodes_level.begin(), data_nodes_level.begin() + 100);
    //poison_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal,&num_total_poisoning,&changed_data,0.9);
    // std::vector<Node*> nodes_by_level2;

    // get_nodes_by_level( &nodes, &nodes_by_level2, 2);

    // for(int i =0; i < nodes_by_level2.size(); i ++){
    //     std::cout<< nodes_by_level2[i]->size << " , ";
    // }

    //GETTING THE LOOP UPS

    // std::vector<KEY_TYPE> lookup_keys = get_search_keys(legitimate_data, values_size, 1000000, seed);
    // // save_data_bin(lookup_keys,lookups_output);
    // // lookup_keys.clear();
    // // lookup_keys = read_data_bin(lookups_output);

    // std::vector<KEY_TYPE> lookup_keys_zipf = get_search_keys_zipf(legitimate_data, values_size, 1000000, seed);
    // save_data_bin(lookup_keys_zipf,lookups_zip_output);
    // lookup_keys_zipf.clear();
    // lookup_keys_zipf = read_data_bin(lookups_zip_output);

    // for(int i :lookup_keys ){
    //     std::cout<< i << " , ";
    // }
    // for(int i :lookup_keys_zipf ){
    //     std::cout<< i << " | ";
    // }
    // for(KEY_TYPE i : altered_data){
    //     std::cout << i << ", ";
    // }
    // if(poison){
    //     save_data_bin(altered_data,changed_data_output);
    // }
    // // std::vector<KEY_TYPE> new_altered = read_data_bin(changed_data_output);

    // if(!poison&benchmark){
    //     altered_data.clear();
    //     altered_data = read_data_bin(changed_data_output);
    // }
    // for(KEY_TYPE i : new_altered){
    //     std::cout << i << ", ";
    // }

    if(false){
         
        //benchmark_alex_real(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);
        //size_t index_size_before =  index_original.index_size(true,false);

        //Change this size later
        size_t index_size_after =  index.total_size();
        //size_t index_size_after = 50000;
        //size_t index_size_change =  index_size_before - index_size_after;
        // std::ofstream file;
        // file.open(model_output+".txt", std::ios_base::app);
    
        // file << "SALI" << ";" << data_name << ";"  << insert << ";" << insert_threshold << ";"<< legitimate_data.size() << ";" << orignal_P << ";"
        // << num_total_poisoning << ";"<< max_model_height << ";" << index_size_after << ";" <<
        // total_node_count << ";" << total_data_count << ";" << std::endl;

        // file.close();

         std::string model_results = "SALI;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" + std::to_string(tra_time)
         + ";" + std::to_string(loss_value)+ ";" + std::to_string(mse) ;

        // std::cout << test_output <<std::endl; 
        saveToCSV(model_results,output_model_csv);

        


    //     std::cout << "Performance Before " <<std::endl;
    //      benchmark_lipp_real(index_original, changed_data, changed_data, "LIPP", data_name+ "_changed", poi_thres, original_output, 
    //      poison,insert,orignal_P,insert_threshold);
    //      benchmark_lipp_real(index_original, lookup_keys, lookup_keys, "LIPP", data_name+"_lookup", poi_thres, original_output,
    //      poison,insert,orignal_P,insert_threshold);
    //      benchmark_lipp_real(index_original, lookup_keys_zipf, lookup_keys_zipf, "LIPP", data_name+"_lookup_zipf", poi_thres, original_output,
    //      poison,insert,orignal_P,insert_threshold);
    //     // benchmark_lipp_real(index_original, lookup_keys_zipf, lookup_keys_zipf, "LIPP", data_name, poi_thres, original_output);
        
    //     //benchmark_lipp_real(index_original,changed_data, changed_data, "LIPP", data_name, poi_thres, original_changed_output);
    //     //index_original.bulk_load(values, legitimate_data.size());
       
    // //    benchmark_alex_real_seperate(index_original, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, original_output);
    // //     benchmark_alex_real_seperate(index_original, changed_data, changed_data, "ALEX", data_name, poi_thres, original_changed_output);
       
    //     // benchmark_alex_real(index_original, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, original_output);
    //     // benchmark_alex_real(index_original, changed_data, changed_data, "ALEX", data_name, poi_thres, original_changed_output);
        
    //     // //Something is happening when I enter these data to the nfl the data changes
    //     // //ISSUE is that after the transformation all the data becomes 0 
    //     // //Some issue with the data type
    //      std::cout << "Index size " << index_original.index_size(true,false) <<std::endl;

        std::cout << "Performance " << std::endl;
        //benchmark_lipp_real(index, legitimate_data, legitimate_data, "LIPP", data_name, poi_thres, poisoned_output);
        // benchmark_lipp_real(index,lookup_keys_zipf, lookup_keys_zipf, "LIPP", data_name, poi_thres, poisoned_changed_output);
    if(!insert){

    
        //=================
         benchmark_sali_real(index,altered_data, altered_data, "SALI", data_name+ "_changed", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index,lookup_keys, lookup_keys, "SALI", data_name+"_lookup", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
         benchmark_sali_real(index,lookup_keys_zipf, lookup_keys_zipf, "SALI", data_name+ "_lookup_zipf", poi_thres, performance_output,
         poison,insert,orignal_P,insert_threshold);
         
        // std::cout << "Index size " << index.index_size(true,false) <<std::endl;
    }
        //==================

        // benchmark_alex_real_seperate(index, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, original_output);
        // benchmark_alex_real_seperate(index, changed_data, changed_data, "ALEX", data_name, poi_thres, original_changed_output);

    //    benchmark_alex_real_seperate(index, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, poisoned_output);
    //     benchmark_alex_real_seperate(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);

         
        // benchmark_alex_real(index, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, poisoned_output);
        // benchmark_alex_real(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);

        // std::ofstream file;
        // file.open(model_output+".txt", std::ios_base::app);
    
        // file << "LIPP" << ";" << data_name << ";"  << insert << ";" << insert_threshold << ";"<< legitimate_data.size() << ";" << orignal_P << ";"
        // << num_total_poisoning << ";"<< max_model_height << ";" << index_size_after << ";" <<
        // total_node_count << ";" << total_data_count << ";" << std::endl;

        // file.close();
       
        //std::cout << "Performance Other " <<std::endl;

       

        // std::ofstream file;
        // file.open(model_output+".txt", std::ios_base::app);
        // //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
        // file << "ALEX " << ";" << data_name2 << ";" << data_size << ";"<< method <<";" << original_poi_thres << ";" << num_total_poisoning <<";"<< original_nodes_number << ";" << original_model_nodes_number << ";" << original_data_nodes_number << ";" 
        // << nodes_fin.size() << ";" << model_nodes_fin.size() << ";" << data_nodes_fin.size() <<  ";"  << replaced_models_number << ";" << tra_time << std::endl;

        // file.close();
        
    //     std::vector<KEY_TYPE> temp_child;
    //     std::vector<PAYLOAD_TYPE> temp_child_payload;
    //      for(int i =0; i < static_cast<int>((data_nodes_fin).size());i++){
    //    // for(int i =0; i < 500;i++){
    // //for(data_node_type* node : (*data_nodes)){
    //          data_node_type* node = (data_nodes_fin)[i];
    //     //std::cout << "Poisoning node " << std::endl;

    //         data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
    //             data_it(node,0);

            

    //         while(!data_it.is_end()){
    //             //nodes[cur_in] = node_it.current();

    //             temp_child.push_back(data_it.key());
    //             temp_child_payload.push_back(data_it.payload());

    //             //check if 0 if ok
    //             data_it.operator++(0);
    //             // if(node_it.is_end()){
    //             //     break;
    //             // }
    //         }
    //      }
    //      std::sort(temp_child.begin(), temp_child.end());
    //      std::sort(legitimate_data.begin(), legitimate_data.end());

    //     std::cout << "Data in nodes " << temp_child.size() << std::endl;
    //     std::cout << "Original size " << legitimate_data.size() << std::endl;
    //      for(int i =0; i < temp_child.size();i++){
    //         if(legitimate_data[i] != temp_child[i]){
    //             std::cout << "Not same " << i << " L " << legitimate_data[i] << " N " << temp_child[i] << std::endl;
    //         }
    //      }

      //   std::cout << "Data in nodes " << temp_child.size() << std::endl;

    // std::vector<size_t> indices;

    // for (size_t i = 0; i < legitimate_data.size(); ++i) {
    //     if (legitimate_data[i] == 44040735) {
    //         indices.push_back(i);
    //     }
    // }
    // std::cout << "Data in nodes " << indices.size() << std::endl;

        
    //DynamicPGMIndex
        // I need to create the index here and then evaluate
        // //For PGM
        // //===========
    //pgm::PGMIndex<KEY_TYPE, 4,4> index_pgm;
    //pgm::DynamicPGMIndex<KEY_TYPE, PAYLOAD_TYPE,index_pgm> index_pgmd(legitimate_data);

    // std::vector<std::pair<uint32_t, uint32_t>> data(1000000);
    // std::generate(data.begin(), data.end(), [] { return std::make_pair(std::rand(), std::rand()); });
    // std::sort(data.begin(), data.end());
    // std::vector<std::pair<KEY_TYPE, PAYLOAD_TYPE>> data(1000000);
    // std::generate(data.begin(), data.end(), [] { return std::make_pair(std::rand(), std::rand()); });
    // std::sort(data.begin(), data.end());
    // pgm::DynamicPGMIndex<KEY_TYPE, PAYLOAD_TYPE> dynamic_pgm(data.begin(),data.end());
    
    // dynamic_pgm.insert_or_assign(2, 4);
    
        //benchmark_pgm<4>(index_pgm,legitimate_data, legitimate_data, "PGM_4", data_name, poi_thres, original_output,false);
    // benchmark_pgm<4>(index_pgm,changed_data, changed_data, "PGM_4", data_name, poi_thres, original_changed_output,false);

        // // //FOR LIPP
        // // //===========
        //LIPP<KEY_TYPE, PAYLOAD_TYPE> lipp;
        //lipp.bulk_load(values, legitimate_data.size()); 

        //benchmark_lipp_real(lipp, legitimate_data, legitimate_data, "LIPP", data_name, poi_thres, original_output);
        // benchmark_lipp_real(lipp, changed_data, changed_data, "LIPP", data_name, poi_thres, original_changed_output);


        // benchmark_nfl(legitimate_data,legitimate_data,"NFL",data_name, poi_thres,original_output);
        // benchmark_nfl(changed_data,changed_data,"NFL",data_name, poi_thres,original_changed_output);

    }

    std::cout << "Max model node height  " << max_model_height << std::endl;
    std::cout << "Time  " << tra_time/1000000000.0 << "s" << std::endl;
    std::cout << "Index Build Time  " << index_time/2000000000.0 << "s" << std::endl;
    // uint64_t myUint64 = 12345678901234567890;

    // // Use the uint64_t directly
    // double myDouble = static_cast<double>(myUint64);

    // std::cout << "Uint64_t: " << myUint64 << std::endl;
    // std::cout << "Double: " << myDouble << std::endl;
    // for (size_t i = 0; i < legitimate_data.size(); ++i) {
    //     std::cout  << values[i].first << " , ";
    // }

    std::cout << "END " <<std::endl;
    
}
//FUNCTIONS NEEDED FOR LIPP
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

//The functions required
//==============================================
//==============================================
//==============================================
bool compare_searches_des(const data_node_type* a, const data_node_type* b) {
        return a->expected_exp_search_iterations_ > b->expected_exp_search_iterations_;
}

bool compare_searches(const data_node_type* a, const data_node_type* b) {
        return a->expected_exp_search_iterations_ < b->expected_exp_search_iterations_;
}

bool compare_models_new(const model_node_type* a, const model_node_type* b) {
        return a->cost_diff_ < b->cost_diff_;
}

bool compare_models(const model_node_type* a, const model_node_type* b) {
    // Sort by level in descending order
    if (a->level_ != b->level_) {
        return a->level_ > b->level_;
    }
    // If levels are equal, sort by cost in ascending order
    else {
        return a->cost_diff_ < b->cost_diff_;
    }
}

bool compare_models_old(const model_node_type* a, const model_node_type* b) {
    // Sort by level in descending order
    if (a->level_ != b->level_) {
        return a->level_ > b->level_;
    }
    // If levels are equal, sort by cost in ascending order
    else {
        return a->cost_ < b->cost_;
    }
}


std::vector<std::pair<model_node_type*,double>> calculate_best_node_cost_prev(std::vector<model_node_type*>& model_nodes_useful)
{
        // model_node_type* best_node;
        // std::vector<model_node_type*> model_nodes_useful;
        // std::vector<double> children_sum;
       // int cur_in3=0;
        
        double cost_diff;
        double temp_cost_diff;

        model_node_type* best_node;
        std::vector<std::pair<model_node_type*,double>> costs_diff_models;
        
        // //Get max height
        // int max_height = 0;
        // for(int i =0 ; i< model_nodes_useful.size(); i++){

        //     if(max_height<model_nodes_useful[i]->level_){
        //         max_height = model_nodes_useful[i]->level_;
        //     }
            
        // }
        //Find the best node to use (the one with the least cost diff)
        //Create a method get best model node
        for(int i =0 ; i< model_nodes_useful.size(); i++){
        //for(model_node_type* node : model_nodes_useful){
            //get child nodes and check if they are data nodes
            // alex::AlexNode<KEY_TYPE, PAYLOAD_TYPE>** children = node->children_;
        
            // int num_child = node->num_children_;
            if(i == 0){
                cost_diff = model_nodes_useful[i]->cost_ - model_nodes_useful[i]->children_cost;
                best_node = model_nodes_useful[i];
            }

            double temp_cost_diff = model_nodes_useful[i]->cost_ - model_nodes_useful[i]->children_cost;
            
            costs_diff_models.push_back(std::make_pair(model_nodes_useful[i], temp_cost_diff));

            if(cost_diff>temp_cost_diff){
                cost_diff = temp_cost_diff;
                best_node = model_nodes_useful[i];
                //I should probably get the index or something (How can I find the parent?)
                //Need to change the NodeIterator to get parents as well
            }


            // for(int i =0; i <num_child;i++){
            //     std::cout << children[i]->is_leaf_<< ", " ;
            // }
            //std::cout << "Model cost: " <<node->cost_<< " Childs cost: " <<children_sum[cur_in3]<< " diff " << temp_cost_diff<< std::endl;
            //cur_in3++;
        }
    return  costs_diff_models;
}
// model_node_type* calculate_best_node_cost_new(std::vector<model_node_type*>& model_nodes_useful, std::vector<double>& children_sum, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index,auto values,int size)
// {
//         // model_node_type* best_node;
//         // std::vector<model_node_type*> model_nodes_useful;
//         // std::vector<double> children_sum;
//         int cur_in3=0;
        
//         double cost_diff;
//         double temp_cost_diff;

//         model_node_type* best_node;


//         //Find the best node to use (the one with the least cost diff)
//         //Create a method get best model node
//         for(int i =0 ; i< model_nodes_useful.size(); i++){
//         //for(model_node_type* node : model_nodes_useful){
//             //get child nodes and check if they are data nodes
//             // alex::AlexNode<KEY_TYPE, PAYLOAD_TYPE>** children = node->children_;
        
//             // int num_child = node->num_children_;
//             auto new_data_node = new (data_node_type::alloc_type(index.get_allocator()).allocate(1))
//             data_node_type(model_nodes_useful[i]->level_, index.derived_params_.max_data_node_slots,
//                             index.key_comp(), index.get_allocator());

//             alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

//             //Need to change values
//             //bulk_load this node
//             new_data_node->bulk_load(values, size, data_node_model,
//                                 index.params_.approximate_model_computation);

//         std::cout << "DONE: New data node created" << std::endl;
//         std::cout << "============================" << std::endl;

        
//         // alex::LinearModel<KEY_TYPE> data_node_model2;

//         // data_node_type::build_model(poisoned_values, poisoned_model_data_size, &data_node_model2,
//         //                             index.params_.approximate_model_computation);
//         //I didnt use &stats
//         //Check the cost etc.

//             new_data_node->cost_ = data_node_type::compute_new_cost(
//             values, size, best_node->level_, data_node_type::kInitDensity_,
//             index.params_.expected_insert_frac, data_node_model,
//             index.params_.approximate_cost_computation);
            
//             if(cur_in3 == 0){
//                 cost_diff = model_nodes_useful[i]->cost_ - children_sum[cur_in3];
//                 best_node = model_nodes_useful[i];
//             }

//             double temp_cost_diff = model_nodes_useful[i]->cost_ - children_sum[cur_in3];
//             if(cost_diff>temp_cost_diff){
//                 cost_diff = temp_cost_diff;
//                 best_node = model_nodes_useful[i];
//                 //I should probably get the index or something (How can I find the parent?)
//                 //Need to change the NodeIterator to get parents as well
//             }


//             // for(int i =0; i <num_child;i++){
//             //     std::cout << children[i]->is_leaf_<< ", " ;
//             // }
//             //std::cout << "Model cost: " <<node->cost_<< " Childs cost: " <<children_sum[cur_in3]<< " diff " << temp_cost_diff<< std::endl;
//             cur_in3++;
//         }
//     return best_node;
// }


void calculate_children_new_cost(std::vector<model_node_type*>& model_nodes_useful, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const)
{

    for(int m =0; m < model_nodes_useful.size(); m++){
        int num_child_best = model_nodes_useful[m]->num_children_;
        node_type** best_children = model_nodes_useful[m]->children_;
        std::mt19937_64 gen_payload(std::random_device{}());
    
        for(int i = 0; i < num_child_best;i++){
            //alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>::Iterator data_it(&best_children[i]);
            // if(best_children[i] == nullptr){
            //     continue;
            // }
                data_node_type* child = static_cast<data_node_type*>(best_children[i]);

                // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                // data_it(static_cast<data_node_type*>(best_children[i]),0);
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;


            while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                //model_node_data.push_back(data_it.key());
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(data_it.key());

                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
            }
            //child_node_data.push_back(temp_child);
            //static_cast<data_node_type*>(best_children[i]);

            auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

        for (int i = 0; i < temp_child.size(); i++) {
            child_values[i].first = temp_child[i];
            child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        }

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            child->cost_ = child->compute_new_cost(
            child_values, temp_child.size(),temp_child.size(),model_nodes_useful[m]->num_keys_model_,child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            index.params_.expected_insert_frac, &child->model_,
            index.params_.approximate_cost_computation);

            // child->cost_ = data_node_type::compute_expected_cost(
            // child_values, temp_child.size(), data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            //best_children_cost += child->cost_;
        }
    }
}

void calculate_data_node_new_cost(std::vector<data_node_type*>& data_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const)
{

    std::mt19937_64 gen_payload(std::random_device{}());
    for(int m =0; m < data_nodes.size(); m++){
        
        
    
        
            //alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>::Iterator data_it(&best_children[i]);
            // if(best_children[i] == nullptr){
            //     continue;
            // }
                data_node_type* child = data_nodes[m];

                // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                // data_it(static_cast<data_node_type*>(best_children[i]),0);
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;


            while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                //model_node_data.push_back(data_it.key());
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(data_it.key());

                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
            }
            //child_node_data.push_back(temp_child);
            //static_cast<data_node_type*>(best_children[i]);

            auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

        for (int i = 0; i < temp_child.size(); i++) {
            child_values[i].first = temp_child[i];
            child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        }

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

             child->cost_ = child->compute_new_cost_new(
            child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            index.params_.expected_insert_frac, &child->model_,
            index.params_.approximate_cost_computation);

            // child->cost_ = data_node_type::compute_expected_cost(
            // child_values, temp_child.size(), data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            //best_children_cost += child->cost_;
        
    }
}
void calculate_data_node_new_cost_new(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const)
{

     std::mt19937_64 gen_payload(std::random_device{}());
     std::vector<double> new_costs;
     data_node_type* temp;
     

    for(int m =0; m < model_nodes.size(); m++){
        
       
        std::vector<node_type*> used_children;
        bool add = false;
        double children_cost =0;
            //alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>::Iterator data_it(&best_children[i]);
            // if(best_children[i] == nullptr){
            //     continue;
            // }
            int num_child = model_nodes[m]->num_children_;
            node_type** children = model_nodes[m]->children_;
           std::vector<KEY_TYPE> child_node_data;

            for(int i =0; i < num_child; i++){
                 

            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                // data_it(static_cast<data_node_type*>(best_children[i]),0);
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                    //std::cout << "NEW CHILD" << std::endl;
                }
                else{
                    //std::cout << "OLD CHILD" << std::endl;
                    add = false;
                }

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;
            if(add){
                while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                //model_node_data.push_back(data_it.key());
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(data_it.key());
                child_node_data.push_back(data_it.key());

                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
                }
                auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                for (int i = 0; i < temp_child.size(); i++) {
                    child_values[i].first = temp_child[i];
                    child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                }

                    // child->cost_ = data_node_type::compute_new_cost(
                    // child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                    // index.params_.expected_insert_frac, &child->model_,
                    // index.params_.approximate_cost_computation);

                    child->cost_ = child->compute_new_cost(
                    child_values, temp_child.size(),temp_child.size(),model_nodes[m]->num_keys_model_ ,child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                    index.params_.expected_insert_frac, &child->model_,
                    index.params_.approximate_cost_computation);

                    children_cost += child->cost_;

            }
            
            //child_node_data.push_back(temp_child);
            //static_cast<data_node_type*>(best_children[i]);

          
        // child_node_data.push_back(temp_child);
        

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),child->level_, data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            // child->cost_ = data_node_type::compute_expected_cost(
            // child_values, temp_child.size(), data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            //best_children_cost += child->cost_;
            }
            model_nodes[m]->children_cost = children_cost;

    }
}

std::vector<double>  calculate_model_node_new_cost(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const)
{

     std::mt19937_64 gen_payload(std::random_device{}());
     std::vector<double> new_costs;
     data_node_type* temp;
     

    for(int m =0; m < model_nodes.size(); m++){
        
       
        std::vector<node_type*> used_children;
        bool add = false;
            //alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>::Iterator data_it(&best_children[i]);
            // if(best_children[i] == nullptr){
            //     continue;
            // }
            int num_child = model_nodes[m]->num_children_;
            node_type** children = model_nodes[m]->children_;
           std::vector<KEY_TYPE> child_node_data;

            for(int i =0; i < num_child; i++){
                 

            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                // data_it(static_cast<data_node_type*>(best_children[i]),0);
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                    //std::cout << "NEW CHILD" << std::endl;
                }
                else{
                    //std::cout << "OLD CHILD" << std::endl;
                    add = false;
                }

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;
            if(add){
                while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                //model_node_data.push_back(data_it.key());
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(data_it.key());
                child_node_data.push_back(data_it.key());

                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
                }
                //auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                // for (int i = 0; i < temp_child.size(); i++) {
                //     child_values[i].first = temp_child[i];
                //     child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                // }

                //     // child->cost_ = data_node_type::compute_new_cost(
                //     // child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                //     // index.params_.expected_insert_frac, &child->model_,
                //     // index.params_.approximate_cost_computation);

                //     child->cost_ = child->compute_new_cost(
                //     child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                //     index.params_.expected_insert_frac, &child->model_,
                //     index.params_.approximate_cost_computation);

            }
            
            //child_node_data.push_back(temp_child);
            //static_cast<data_node_type*>(best_children[i]);

          
        // child_node_data.push_back(temp_child);
        

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),child->level_, data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            // child->cost_ = data_node_type::compute_expected_cost(
            // child_values, temp_child.size(), data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            //best_children_cost += child->cost_;
            }
              auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[child_node_data.size()];

            for (int i = 0; i < child_node_data.size(); i++) {
                child_values[i].first = child_node_data[i];
                child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            }

            auto new_data_node = new (data_node_type::alloc_type(index.get_allocator()).allocate(1))
            data_node_type(model_nodes[m]->level_, index.derived_params_.max_data_node_slots,
                            index.key_comp(), index.get_allocator());

            alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

            new_data_node->bulk_load(child_values, child_node_data.size(), data_node_model,
                                index.params_.approximate_model_computation);


            // model_nodes[m]->cost_ = data_node_type::compute_new_cost(
            //     child_values, child_node_data.size(),child_node_data.size(), model_nodes[m]->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            //     index.params_.expected_insert_frac, data_node_model,
            //     index.params_.approximate_cost_computation);
            model_nodes[m]->cost_ = new_data_node->compute_new_cost(
                child_values, child_node_data.size(),child_node_data.size(),child_node_data.size(), model_nodes[m]->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                index.params_.expected_insert_frac, data_node_model,
                index.params_.approximate_cost_computation);

                model_nodes[m]->num_keys_model_ = child_node_data.size();
            

        new_costs.push_back(model_nodes[m]->cost_);
    }
    return new_costs;
}

std::vector<double> calculate_model_node_poisoning(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double pois_thres,double search_const, double traversal_const)
{
        std::mt19937_64 gen_payload(std::random_device{}());

        std::vector<double> pois_costs; 

        data_node_type* temp;

        
        //model_node_type* best_node;
    for(int m =0; m < model_nodes.size(); m++){
        
        
        std::vector<node_type*> used_children;
        bool add = false;
    
        
            //alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>::Iterator data_it(&best_children[i]);
            // if(best_children[i] == nullptr){
            //     continue;
            // }
            int num_child = model_nodes[m]->num_children_;
            node_type** children = model_nodes[m]->children_;
            std::vector<KEY_TYPE> child_node_data;
            double children_cost =0;

           

            for(int i =0; i < num_child; i++){
                 

            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                // data_it(static_cast<data_node_type*>(best_children[i]),0);
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                    //std::cout << "NEW CHILD" << std::endl;
                }
                else{
                    //std::cout << "OLD CHILD" << std::endl;
                    add = false;
                }

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;

            if(add){
                while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                //model_node_data.push_back(data_it.key());
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(data_it.key());
                
                child_node_data.push_back(data_it.key());

                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
                }
                auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                    for (int i = 0; i < temp_child.size(); i++) {
                        child_values[i].first = temp_child[i];
                        child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                    }

                    // child->cost_ = data_node_type::compute_new_cost(
                    // child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_, search_const, traversal_const,
                    // index.params_.expected_insert_frac, &child->model_,
                    // index.params_.approximate_cost_computation);
                    // children_cost += child->cost_;
                    child->cost_ = child->compute_new_cost(
                    child_values, temp_child.size(),temp_child.size(),model_nodes[m]->num_keys_model_,child->level_, data_node_type::kInitDensity_, search_const, traversal_const,
                    index.params_.expected_insert_frac, &child->model_,
                    index.params_.approximate_cost_computation);
                    children_cost += child->cost_;
            }
            
            //child_node_data.push_back(temp_child);
            //static_cast<data_node_type*>(best_children[i]);

          
        // child_node_data.push_back(temp_child);
        

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),child->level_, data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            // child->cost_ = data_node_type::compute_expected_cost(
            // child_values, temp_child.size(), data_node_type::kInitDensity_,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);

            //best_children_cost += child->cost_;
            }
              
            std::vector<KEY_TYPE> poisoned_model_data = perform_poisoning(child_node_data, pois_thres);


            auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[poisoned_model_data.size()];

            for (int i = 0; i < poisoned_model_data.size(); i++) {
                child_values[i].first = poisoned_model_data[i];
                child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            }

            auto new_data_node = new (data_node_type::alloc_type(index.get_allocator()).allocate(1))
            data_node_type(model_nodes[m]->level_, index.derived_params_.max_data_node_slots,
                            index.key_comp(), index.get_allocator());

            alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

            new_data_node->bulk_load(child_values, poisoned_model_data.size(), data_node_model,
                                index.params_.approximate_model_computation);


            // model_nodes[m]->cost_ = data_node_type::compute_new_cost(
            //     child_values, poisoned_model_data.size(),child_node_data.size(), model_nodes[m]->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            //     index.params_.expected_insert_frac, data_node_model,
            //     index.params_.approximate_cost_computation);

            model_nodes[m]->cost_ = new_data_node->compute_new_cost(
                child_values, poisoned_model_data.size(),child_node_data.size(),child_node_data.size(), model_nodes[m]->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                index.params_.expected_insert_frac, data_node_model,
                index.params_.approximate_cost_computation);

            pois_costs.push_back(model_nodes[m]->cost_);
            model_nodes[m]->children_cost = children_cost;

            // if(m ==0){
            //     best_cost = model_nodes[m]->cost_;
            //     //best_node = model_nodes[m];
            // }
            //if(best_cost < )
    }
    return pois_costs;
}

void calculate_cost_difference(std::vector<model_node_type*>& model_nodes_useful)
{

     std::vector<double> cost_diff;

    for(int i =0; i < model_nodes_useful.size(); i++){


            model_nodes_useful[i]->cost_diff_ = model_nodes_useful[i]->cost_ - model_nodes_useful[i]->children_cost;
            
            //costs_diff_models.push_back(std::make_pair(model_nodes_useful[i], temp_cost_diff))   
    
    }
    //return cost_diff;
}

std::vector<model_node_type*> can_be_done(std::vector<model_node_type*>& model_nodes_useful)
{

     std::vector<double> cost_diff;
     std::vector<model_node_type*> models_reduce;

    for(int i =0; i < model_nodes_useful.size(); i++){

        if(model_nodes_useful[i]->cost_diff_<0){
            models_reduce.push_back(model_nodes_useful[i]);
        }
            //model_nodes_useful[i]->cost_diff_ = model_nodes_useful[i]->cost_ - model_nodes_useful[i]->children_cost;
            
            //costs_diff_models.push_back(std::make_pair(model_nodes_useful[i], temp_cost_diff))   
    
    }
    return models_reduce;
}

std::vector<std::pair<model_node_type*,double>> calculate_best_node_cost_reduction(std::vector<model_node_type*>& model_nodes_useful, std::vector<double> new_costs,
std::vector<double> pois_costs )
{
        model_node_type* best_node;
        // std::vector<model_node_type*> model_nodes_useful;
        // std::vector<double> children_sum;
        double best_reduction;
         std::vector<std::pair<model_node_type*,double>> costs_diff_models;

        for(int i =0 ; i< model_nodes_useful.size(); i++){
        //for(model_node_type* node : model_nodes_useful){
            //get child nodes and check if they are data nodes
            // alex::AlexNode<KEY_TYPE, PAYLOAD_TYPE>** children = node->children_;
        
            // int num_child = node->num_children_;
            if(i == 0){
                best_reduction = pois_costs[i] - new_costs[i];
                best_node = model_nodes_useful[i];
            }

            double temp_cost_diff = pois_costs[i] - new_costs[i];
            costs_diff_models.push_back(std::make_pair(model_nodes_useful[i], temp_cost_diff));
            
            if(best_reduction>temp_cost_diff){
                best_reduction = temp_cost_diff;
                best_node = model_nodes_useful[i];
                //I should probably get the index or something (How can I find the parent?)
                //Need to change the NodeIterator to get parents as well
            }


            // for(int i =0; i <num_child;i++){
            //     std::cout << children[i]->is_leaf_<< ", " ;
            // }
            //std::cout << "Model cost: " <<node->cost_<< " Childs cost: " <<children_sum[cur_in3]<< " diff " << temp_cost_diff<< std::endl;
        }
    return costs_diff_models;
}

// //Method to get the constants of search per search and traversal per level
// void get_constants(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<KEY_TYPE> data, int max_data, double *const_search, double *const_traversal){
//     int const_size = std::min(static_cast<int>(round(data.size()*1)),max_data);

//     //search cost per search
//     // double search_const =0;

//     // //travel cost per level
//     // double traversal_const =0;
//     std::vector<double> const_time_trav(const_size,0);
//     std::vector<double> const_time_search(const_size,0);

//     //To get the random index
//     std::uniform_int_distribution<size_t> rand_dis(0, data.size() - 1);

//     std::random_device rd;
//     std::mt19937 generator(rd());
    

//     //std::cout << "Calculating the constants" << std::endl;
//     for(int i =0; i < const_size; i++){

//         int idx = rand_dis(generator);

//         auto start_q = std::chrono::high_resolution_clock::now();

//         //auto payload = index.get_payload(lookups[i]);
//         // auto payload = index.find(lookups[i]);

//         data_node_type* leaf = index.get_leaf(data[idx]);

//         auto end_q = std::chrono::high_resolution_clock::now();

//         auto start_search = std::chrono::high_resolution_clock::now();

//         //  std::cout << "Leaf :"<< leaf << std::endl;
//         //   std::cout << "Lookup :"<< lookups[i] << std::endl;

//         int location_key = leaf->find_key(data[idx]);

        
        

//         auto end_search = std::chrono::high_resolution_clock::now();

//         double level = leaf->level_;
//         double num_searches = leaf->number_of_searches_; 
//         //std::cout << "Adding to it " << level << std::endl;

//         const_time_trav[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end_q - start_q).count()/(1.0*level);

//         //I added +1 for the initial one as well
//         //const_time_search[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end_search - start_search).count()/(num_searches+1);
//         if(num_searches==0){
//             const_time_search[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end_search - start_search).count();

//         }
//         else{
//             const_time_search[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end_search - start_search).count()/(1.0*num_searches);
//         }
//         //std::cout << "Added" << std::endl;

//     }

//     *const_search = std::accumulate(const_time_search.begin(), const_time_search.end(), 0.0)/static_cast<double>(const_size);
//     *const_traversal = std::accumulate(const_time_trav.begin(), const_time_trav.end(), 0.0)/static_cast<double>(const_size);

// }

//Method to get the constants of search per search and traversal per level
void get_constants(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<KEY_TYPE> data, int max_data, double *const_search, double *const_traversal){
    int const_size = std::min(static_cast<int>(round(data.size()*1)),max_data);

    //search cost per search
    // double search_const =0;

    // //travel cost per level
    // double traversal_const =0;
    // std::vector<double> const_time_trav(const_size,0);
    // std::vector<double> const_time_search(const_size,0);

    //To get the random index
    std::uniform_int_distribution<size_t> rand_dis(0, data.size() - 1);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::vector<KEY_TYPE> data_test;

    //get the random indexes to search
    for(int i =0; i < const_size; i++){
        //int idx = rand_dis(generator);
        data_test.push_back(i);
    }
    

    //std::cout << "Calculating the constants" << std::endl;
    //Calculating just the query
    int levels;
    auto start_q = std::chrono::high_resolution_clock::now();
    
    for(int i =0; i < const_size; i++){

        int idx = data_test[i];
      
        data_node_type* leaf = index.get_leaf(data[idx]);

        levels += static_cast<int>(leaf->level_);
       
    
    }
    auto end_q = std::chrono::high_resolution_clock::now();

    auto time_query = std::chrono::duration_cast<std::chrono::nanoseconds>(end_q - start_q).count();

    //Calculating the query time per level
    std::cout << "const_size " << const_size << std::endl;
    std::cout << "Levels " << levels << std::endl;
    std::cout << "time_query " << time_query << std::endl;
    
    if(levels==0){
            *const_traversal = time_query;
        }
        else{
            *const_traversal = time_query/(1.0*levels);
            //*const_traversal = time_query;
        }

    //Calculating the total time
     long num_searches =0;

     auto start_tot = std::chrono::high_resolution_clock::now();
     for(int i =0; i < const_size; i++){

        int idx = data_test[i];
 

        data_node_type* leaf = index.get_leaf(data[idx]);


        int location_key = leaf->find_key(data[idx]);

        //Get the number of searches too
        num_searches += leaf->number_of_searches_; 

    }
     auto end_tot = std::chrono::high_resolution_clock::now();

     auto time_tot = std::chrono::duration_cast<std::chrono::nanoseconds>(end_tot - start_tot).count();

     auto search_time = time_tot - time_query;

     if(num_searches==0){
            *const_search = search_time;
        }
        else{
            *const_search = search_time/(1.0*num_searches);
        }


    //*const_search = *const_search/static_cast<double>(const_size);
    //*const_traversal = *const_traversal/static_cast<double>(const_size);

}

//To get if the dataset has duplicates or not
bool check_duplicates(std::vector<KEY_TYPE> data){
    std::unordered_set<KEY_TYPE> seen;

    bool hasDuplicates = false;

    for (KEY_TYPE num : data) {
        if (seen.find(num) != seen.end()) {
            hasDuplicates = true;
            break;
        }
        seen.insert(num);
    }

    return hasDuplicates;
    // std::vector<KEY_TYPE> sortedVec = data; // Create a copy to avoid modifying the original vector
    // std::sort(sortedVec.begin(), sortedVec.end()); // Sort the vector

    // for (size_t i = 1; i < sortedVec.size(); ++i) {
    //     if (sortedVec[i] == sortedVec[i - 1]) {
    //         return true; // Found a duplicate
    //     }
    // }
    // return false; // No duplicates found
}

//To get all the nodes from the index as well as the max model height
void get_all_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index, std::vector<node_type*> *nodes,std::vector<model_node_type*> *model_nodes,
std::vector<data_node_type*> *data_nodes,int *max_model_height){
    //Iterator to get all the nodes
    alex::Alex<KEY_TYPE, PAYLOAD_TYPE>::NodeIterator node_it(&index);
   
    int cur_in =0;
    
    //get all nodes and save them in the vectors
    while(!node_it.is_end()){
        
        nodes->push_back(node_it.current());

        if((*nodes)[cur_in]->is_leaf_){
            data_nodes->push_back(static_cast<data_node_type*>((*nodes)[cur_in]));
            
        }
        else{
            model_nodes->push_back(static_cast<model_node_type*>((*nodes)[cur_in]));
            
             //Getting the max model node height

            if(*max_model_height < static_cast<model_node_type*>((*nodes)[cur_in])->level_){

                *max_model_height = static_cast<model_node_type*>((*nodes)[cur_in])->level_;
            }
        }
        cur_in ++;
        node_it.next();
    }

}

// To check if all children of a model node are data nodes and calculates the children_sum
bool check_all_data_nodes(model_node_type* node){

    node_type** children = node->children_;
    bool is_all_data_nodes = true;
    int num_child = node->num_children_;
    double child_sum = 0;


    for(int i =0; i < num_child;i++){
        //If the child is not a leaf break
        if(!children[i]->is_leaf_){
            is_all_data_nodes =false;
            break;
        }
        child_sum += children[i]->cost_;
    }

    if(is_all_data_nodes){
        node->children_cost = child_sum;
    }

    return is_all_data_nodes;
    
}

//To get all the usefull aka all children are data nodes models and the models by level
void get_useful_model_nodes( std::vector<model_node_type*> *model_nodes, std::vector<model_node_type*> *model_nodes_useful, 
std::vector<model_node_type*> *model_nodes_by_level){

    for(model_node_type* node : *model_nodes){
            //get child nodes and check if they are data nodes

            //Also save models by their level (so I can find the parent much easily)
            model_nodes_by_level[node->level_].push_back(node);

            //If all are data nodes
            if(check_all_data_nodes(node)){
                model_nodes_useful->push_back(node);
            }
        }
}

//Get the children data and calculate the children data costs expects a model node with all children data nodes
void get_children_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,model_node_type* model, std::vector<KEY_TYPE> *model_node_data, 
std::vector<PAYLOAD_TYPE> *model_node_payload, double *children_cost,
double const_search, double const_traversal, bool calculate_childern_costs = true){

        std::mt19937_64 gen_payload(std::random_device{}());

        node_type** children = model->children_;
        std::vector< std::vector<KEY_TYPE> > child_node_data;
        
        int num_child = model->num_children_;

        std::vector<node_type*> used_children;

        std::vector<node_type*> used_children2;
        bool add = false;

        bool add2 = false;

        int count_model = 0;
       
        //check for each num_children
        // for(int i = 0; i < num_child;i++){
        //     data_node_type* child = static_cast<data_node_type*>(children[i]);

        //     if (std::find((used_children2).begin(), (used_children2).end(), child) == (used_children2).end()){
        //             add2 = true;
        //             used_children2.push_back(child);
        //             //std::cout << "NEW CHILD" << std::endl;
        //         }
        //         else{
        //             //std::cout << "OLD CHILD" << std::endl;
        //             add2 = false;
        //         }
        //     if(add2){
        //         count_model+=child->num_keys_;
        //     }

        // }
    
        for(int i = 0; i < num_child;i++){
            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                // Get the iterator
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;
                std::vector<PAYLOAD_TYPE> temp_child_payload;

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                    //std::cout << "NEW CHILD" << std::endl;
                }
                else{
                    //std::cout << "OLD CHILD" << std::endl;
                    add = false;
                }


            //Check if there are data
            if(add){
                while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                
                KEY_TYPE key = data_it.key();
                PAYLOAD_TYPE payload = data_it.payload();

                // if (std::find((*model_node_data).begin(), (*model_node_data).end(), key) == (*model_node_data).end()) {
                //         // If not found, add it to the original vector
                //         //changed_data.push_back(element);
                //         model_node_data->push_back(key);
                //         model_node_payload->push_back(payload);
                //     }
                model_node_data->push_back(key);
                model_node_payload->push_back(payload);
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(key);
                temp_child_payload.push_back(payload);

                //check if 0 if ok
                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
                }
                //Add to child data

                if(calculate_childern_costs)
                {
                        child_node_data.push_back(temp_child);
                    //static_cast<data_node_type*>(children[i]);


                    auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                    for (int i = 0; i < temp_child.size(); i++) {
                        child_values[i].first = temp_child[i];
                        child_values[i].second = temp_child_payload[i];
                    }

                    // child->cost_ = data_node_type::compute_new_cost(
                    // child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_, const_search, const_traversal,
                    // index.params_.expected_insert_frac, &child->model_,
                    // index.params_.approximate_cost_computation);

                    child->cost_ = child->compute_new_cost(
                    child_values, temp_child.size(),temp_child.size(),count_model,child->level_, data_node_type::kInitDensity_, const_search, const_traversal,
                    index.params_.expected_insert_frac, &child->model_,
                    index.params_.approximate_cost_computation);


                    *children_cost += child->cost_;
                }
                // else{
                //     *children_cost += child->cost_;
                // }
            }
            // while(!data_it.is_end()){
            //     //nodes[cur_in] = node_it.current();
                
            //     KEY_TYPE key = data_it.key();
            //     PAYLOAD_TYPE payload = data_it.payload();

            //     // if (std::find((*model_node_data).begin(), (*model_node_data).end(), key) == (*model_node_data).end()) {
            //     //         // If not found, add it to the original vector
            //     //         //changed_data.push_back(element);
            //     //         model_node_data->push_back(key);
            //     //         model_node_payload->push_back(payload);
            //     //     }
            //     // model_node_data->push_back(key);
            //     // model_node_payload->push_back(payload);
            //     //std::cout << model_node_data[i] <<", ";
            //     temp_child.push_back(key);
            //     temp_child_payload.push_back(payload);

            //     //check if 0 if ok
            //     data_it.operator++(0);
            //     // if(node_it.is_end()){
            //     //     break;
            //     // }
            // }
            // //Add to child data

            // if(calculate_childern_costs)
            // {
            //         child_node_data.push_back(temp_child);
            //     //static_cast<data_node_type*>(children[i]);


            //     auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

            //     for (int i = 0; i < temp_child.size(); i++) {
            //         child_values[i].first = temp_child[i];
            //         child_values[i].second = temp_child_payload[i];
            //     }

            //     child->cost_ = data_node_type::compute_new_cost(
            //     child_values, temp_child.size(),child->level_, data_node_type::kInitDensity_, const_search, const_traversal,
            //     index.params_.expected_insert_frac, &child->model_,
            //     index.params_.approximate_cost_computation);


            //     *children_cost += child->cost_;
            // }
            // else{
            //      *children_cost += child->cost_;
            // }
            // child_node_data.push_back(temp_child);
            // //static_cast<data_node_type*>(children[i]);


            // auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

            // for (int i = 0; i < temp_child.size(); i++) {
            //     child_values[i].first = temp_child[i];
            //     child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            // }

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),child->level_, data_node_type::kInitDensity_, const_search, const_traversal,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);


           
        }
        *children_cost = model->children_cost;

        // if(check_duplicates(*(model_node_data))){
        //     //std::cout << "THERE ARE DUPLICATES" << std::endl;
        // }
        // else{
        //     //std::cout << "THERE ARE NO DUPLICATES" << std::endl;
        // }

}

//To create the values needed for alex key value pairs
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

//create a new data node from the poisoned data
data_node_type * create_new_data_node(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index , model_node_type* model, std::pair<uint64_t, double> * values, std::vector<KEY_TYPE> leg, int size, int cur_size,
double const_search, double const_traversal){

    node_type* model_data = static_cast<node_type*>(model);
    //data_node_type* new_data_node;

    auto new_data_node = new (data_node_type::alloc_type(index.get_allocator()).allocate(1))
            data_node_type(model->level_, index.derived_params_.max_data_node_slots,
                            index.key_comp(), index.get_allocator());

        alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

        //Need to change values
        //bulk_load this node
        //data_node_type* new_data_node;
        //new_data_node->level_ = model->level_;
        //std::cout << "Data capacity : "<< new_data_node->data_capacity_  <<std::endl;

        new_data_node->bulk_load(values, size, data_node_model,
                            index.params_.approximate_model_computation);
        
         

        // new_data_node->bulk_load_new(values, leg,size, data_node_model,
        //                     index.params_.approximate_model_computation);

        //std::cout << "DONE: New data node created" << std::endl;
        //std::cout << "============================" << std::endl;

        
        // alex::LinearModel<KEY_TYPE> data_node_model2;

        // data_node_type::build_model(poisoned_values, poisoned_model_data_size, &data_node_model2,
        //                             index.params_.approximate_model_computation);
        //I didnt use &stats
        //Check the cost etc.

        // new_data_node->cost_ = data_node_type::compute_new_cost(
        //     values, size,cur_size, model->level_, data_node_type::kInitDensity_,const_search, const_traversal,
        //     index.params_.expected_insert_frac, data_node_model,
        //     index.params_.approximate_cost_computation);
        std:: cout << "model size " << model->num_keys_model_ << std::endl;
        std:: cout << "cur size " << cur_size << std::endl;
        new_data_node->cost_ = new_data_node->compute_new_cost(
            values, size,cur_size,cur_size, model->level_, data_node_type::kInitDensity_,const_search, const_traversal,
            index.params_.expected_insert_frac, data_node_model,
            index.params_.approximate_cost_computation);
            
        new_data_node->is_leaf_ =true;


    return new_data_node;
}

std::vector<model_node_type*> find_parent(model_node_type* model,std::vector<model_node_type*> *model_nodes_by_level, int *model_idx
,int *parent_model_idx){

    bool parent_found = false;
    std::vector<model_node_type*> parents;

    int best_level = model->level_;
    int level = best_level -1;
    int count =0;

        for(int l = level ; level >=0 ; level--){
            for(model_node_type* node : model_nodes_by_level[l]){


                node_type** childrens  = node->children_;
                //node_type* childrens_pt = *childrens;

                int num_child = node->num_children_;

                for(int i =0; i< num_child; i++){
                    //std::cout << "Child "<< childrens[i] << " , ";
                    if(model == childrens[i]){
                        //parent = node;
                        parents.push_back(node);
                        parent_found = true;
                        *model_idx = i;
                        count++;
                        //break;
                    }
                }
                if(count >1){
                    //std::cout<< "MULTIPLE PARENTS FOUND" << std::endl;
                }
                if(parent_found){
                    //std::cout << "Parent "<< parent << std::endl;
                    //std::cout << "Parent index "<< *model_idx << std::endl;
                    //std::cout << "Best node "<< model << std::endl;
                    //std::cout << "Parent child " << parent->children_[*model_idx] << std::endl;
                    break;
                }
                *parent_model_idx++;
            
            }
        }
        return parents;

}

void update_data_structure(model_node_type* best_node,std::vector<model_node_type*> *model_nodes,
std::vector<model_node_type*> *model_nodes_by_level, std::vector<model_node_type*>* model_nodes_useful,
std::vector<model_node_type*>* model_nodes_useful_by_level, std::vector<model_node_type*> parents,
alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const, std::vector<model_node_type*>* new_useful_parents){

    //I have to remove this model node from the data structures
    int best_level = best_node->level_;

    //remove from model_nodes
    auto it = std::find((*model_nodes).begin(), (*model_nodes).end(), best_node);

    //int index = std::distance((*model_nodes).begin(), it);

    //Remove from model nodes and best level
    if (it != (*model_nodes).end()) {
        // Pointer found, erase it from the vector
        (*model_nodes).erase(it);
        
       // std::cout << "Pointer removed from the vector." << std::endl;
    } else {
        //std::cout << "Pointer not found in the vector." << std::endl;
    }

    //std::cout << "model nodes number "<< model_nodes_by_level[best_level].size() << std::endl;

    //remove from model nodes by level
    auto it2 = std::find(model_nodes_by_level[best_level].begin(), model_nodes_by_level[best_level].end(), best_node);

    if (it2 != model_nodes_by_level[best_level].end()) {
        // Pointer found, erase it from the vector
        model_nodes_by_level[best_level].erase(it2);
        //std::cout << "Pointer removed from the vector." << std::endl;
    } else {
        //std::cout << "Pointer not found in the vector." << std::endl;
    }

    //remove from model nodes useful
    auto it3 = std::find((*model_nodes_useful).begin(), (*model_nodes_useful).end(), best_node);

    if (it3 != (*model_nodes_useful).end()) {
        // Pointer found, erase it from the vector
        (*model_nodes_useful).erase(it3);
        //std::cout << "Pointer removed from the vector." << std::endl;
    } else {
        //std::cout << "Pointer not found in the vector." << std::endl;
    }

    //remove from model nodes useful by level
    // auto it4 = std::find((*model_nodes_useful_by_level).begin(), (*model_nodes_useful_by_level).end(), best_node);

    // if (it4 != (*model_nodes_useful_by_level).end()) {
    //     // Pointer found, erase it from the vector
    //     (*model_nodes_useful_by_level).erase(it4);
    //     //std::cout << "Pointer removed from the vector." << std::endl;
    // } else {
    //     //std::cout << "Pointer not found in the vector." << std::endl;
    // }

    //Check if the parent is now has all data nodes or not
    // If so, then add to the useful data nodes
    for(int i =0; i< parents.size();i++){
        model_node_type* parent = parents[i];

        if(check_all_data_nodes(parent)){
            //std::cout << "Updated AKA new useful mode created now." << std::endl;
            //Make sure to do the calculations of costs now
            model_nodes_useful->push_back(parent);
            new_useful_parents->push_back(parent);

            //Get the new cost as well as the difference to the parent
            // The new data node cost should also be here now
            double cost;
            get_children_cost(parent, &cost,index, search_const, traversal_const);
            parent->children_cost = cost;
            parent->cost_diff_ = parent->cost_ - cost;
            
        }
        // else{

        // }
    }

    // If not do nothing

}

std::vector<model_node_type*> replace_parents(model_node_type* model,data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level, int *model_idx,
int *parent_model_idx){

    bool parent_found = false;
    model_node_type* parent;
    std::vector<model_node_type*> parents;

    //Start from the level above of the best node
    int level = model->level_ - 1;
    int count = 0;

        //For all levels I have to replace the parent
        for(int l = level ; level >=0 ; level--){
            //for(int j =0; j < model_nodes_by_level[l].size(); j++){
            //For each model node in the model nodes by level (all models not just useful)
            for(model_node_type* node : model_nodes_by_level[l]){

                //model_node_type* node = model_nodes_by_level[l][j];
                node_type** childrens  = node->children_;
                //node_type* childrens_pt = *childrens;

                int num_child = node->num_children_;

                for(int i =0; i< num_child; i++){
                    //std::cout << "Child "<< childrens[i] << " , ";
                    if(model == childrens[i]){
                        parent = node;
                        parents.push_back(node);
                        parent_found = true;
                        *model_idx = i;
                        count++;
                        // std::cout << "INSIDE IF" << std::endl;
                        // std::cout << "Best node "<< new_data_node << std::endl;
                        // std::cout << "Parent "<< parent->children_[i] << std::endl;
                        //Setting the new parent pointer
                        parent->children_[i] = new_data_node;

                        //Maybe I can add nullptr to all the original pointers because now I have new pointers
                        //parent->children_[*model_idx] == nullptr;
                       // std::cout<<"parent : " << parent << std::endl;
                        //std::cout << "Parent childern size "<< num_child << std::endl;
                        // for(int j = 0; j < num_child; j++){
                        //     std::cout<< parent->children_[j] << " ";
                        // }
                        //break;
                    }
                    // if(parent_found){
                    //     break;
                    // }
                }
                // if(count >1){
                //     //std::cout<< "MULTIPLE PARENTS FOUND" << std::endl;
                // }
                // if(parent_found){
                //     std::cout << "Parent "<< parent << std::endl;
                //     std::cout << "Parent index "<< *model_idx << std::endl;
                //     std::cout << "Best node "<< model << std::endl;
                //     std::cout << "Parent child " << parent->children_[*model_idx] << std::endl;
                //     break;
                // }
                *parent_model_idx++;
            
            }
        }
        //std::cout << "Parents "<< count << std::endl;
        
        return parents;

}

void replace_parents_data_node(data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level, int *model_idx,
int *parent_model_idx){

    //bool parent_found = false;
    model_node_type* parent;

    //Start from the level above of the best node
    int level = new_data_node->level_ - 1;
    int count =0;

        //For all levels I have to replace the parent
        for(int l = level ; level >=0 ; level--){
            //for(int j =0; j < model_nodes_by_level[l].size(); j++){
            //For each model node in the model nodes by level (all models not just useful)
            for(model_node_type* node : model_nodes_by_level[l]){

                //model_node_type* node = model_nodes_by_level[l][j];
                node_type** childrens  = node->children_;
                //node_type* childrens_pt = *childrens;

                int num_child = node->num_children_;

                for(int i =0; i< num_child; i++){
                    //std::cout << "Child "<< childrens[i] << " , ";
                    if(new_data_node == childrens[i]){
                        parent = node;
                        //parents.push_back(node);
                        //parent_found = true;
                        *model_idx = i;
                        count++;
                        // std::cout << "INSIDE IF" << std::endl;
                        // std::cout << "Best node "<< new_data_node << std::endl;
                        // std::cout << "Parent "<< parent->children_[i] << std::endl;
                        //Setting the new parent pointer
                        //parent->children_[*model_idx] = new_data_node;

                        //Maybe I can add nullptr to all the original pointers because now I have new pointers
                        parent->children_[*model_idx] == nullptr;
                        //std::cout<<"parent : " << parent << std::endl;
                        //std::cout << "Parent childern size "<< num_child << std::endl;
                        //break;
                    }
                }
                // if(count >1){
                //     //std::cout<< "MULTIPLE PARENTS FOUND" << std::endl;
                // }
                // if(parent_found){
                //     std::cout << "Parent "<< parent << std::endl;
                //     std::cout << "Parent index "<< *model_idx << std::endl;
                //     std::cout << "Best node "<< model << std::endl;
                //     std::cout << "Parent child " << parent->children_[*model_idx] << std::endl;
                //     break;
                // }
                *parent_model_idx++;
            
            }
        }
        //std::cout << "Parents "<< count << std::endl;

}

//Expects the children cost to be already calculated
 void get_children_cost(model_node_type* model, double *children_cost,alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const){

        std::mt19937_64 gen_payload(std::random_device{}());

        node_type** children = model->children_;
        std::vector< std::vector<KEY_TYPE> > child_node_data;
        
        int num_child = model->num_children_;

        std::vector<node_type*> used_children;
        bool add = false;

        std::vector<node_type*> used_children2;

        std::vector<KEY_TYPE> model_data;
        std::vector<PAYLOAD_TYPE> model_payload;

        bool add2 = false;

        int count_model = 0;
       
        //check for each num_children
        for(int i = 0; i < num_child;i++){
            data_node_type* child = static_cast<data_node_type*>(children[i]);

            if (std::find((used_children2).begin(), (used_children2).end(), child) == (used_children2).end()){
                    add2 = true;
                    used_children2.push_back(child);
                    //std::cout << "NEW CHILD" << std::endl;
                }
                else{
                    //std::cout << "OLD CHILD" << std::endl;
                    add2 = false;
                }
            if(add2){
                count_model+=child->num_keys_;
            }

        }
       
        //check for each num_children
    
        for(int i = 0; i < num_child;i++){
            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                // Get the iterator
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;
                std::vector<PAYLOAD_TYPE> temp_child_payload;

                //std::cout << static_cast<data_node_type*>(best_children[i])->num_keys_ << std::endl;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                    //std::cout << "NEW CHILD" << std::endl;
                }
                else{
                    //std::cout << "OLD CHILD" << std::endl;
                    add = false;
                }


            //Check if there are data
            if(add){
                while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();
                
                // KEY_TYPE key = data_it.key();
                // PAYLOAD_TYPE payload = data_it.payload();

                // if (std::find((*model_node_data).begin(), (*model_node_data).end(), key) == (*model_node_data).end()) {
                //         // If not found, add it to the original vector
                //         //changed_data.push_back(element);
                //         model_node_data->push_back(key);
                //         model_node_payload->push_back(payload);
                //     }
                // model_node_data->push_back(key);
                // model_node_payload->push_back(payload);
                //std::cout << model_node_data[i] <<", ";
                temp_child.push_back(data_it.key());
                temp_child_payload.push_back(data_it.payload());

                model_data.push_back(data_it.key());
                model_payload.push_back(data_it.payload());

                //check if 0 if ok
                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
                }
                //Add to child data

                
                        child_node_data.push_back(temp_child);
                    //static_cast<data_node_type*>(children[i]);


                    auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                    for (int i = 0; i < temp_child.size(); i++) {
                        child_values[i].first = temp_child[i];
                        child_values[i].second = temp_child_payload[i];
                    }

                    child->cost_ = child->data_node_type::compute_new_cost(
                    child_values, temp_child.size(),temp_child.size(),count_model,child->level_, data_node_type::kInitDensity_, search_const, 
                    traversal_const,
                    index.params_.expected_insert_frac, &child->model_,
                    index.params_.approximate_cost_computation);


                    *children_cost += child->cost_;
                
            }
            // while(!data_it.is_end()){
            //     //nodes[cur_in] = node_it.current();
                
            //     KEY_TYPE key = data_it.key();
            //     PAYLOAD_TYPE payload = data_it.payload();

            //     // if (std::find((*model_node_data).begin(), (*model_node_data).end(), key) == (*model_node_data).end()) {
            //     //         // If not found, add it to the original vector
            //     //         //changed_data.push_back(element);
            //     //         model_node_data->push_back(key);
            //     //         model_node_payload->push_back(payload);
            //     //     }
            //     // model_node_data->push_back(key);
            //     // model_node_payload->push_back(payload);
            //     //std::cout << model_node_data[i] <<", ";
            //     temp_child.push_back(key);
            //     temp_child_payload.push_back(payload);

            //     //check if 0 if ok
            //     data_it.operator++(0);
            //     // if(node_it.is_end()){
            //     //     break;
            //     // }
            // }
            // //Add to child data

            // if(calculate_childern_costs)
            // {
            //         child_node_data.push_back(temp_child);
            //     //static_cast<data_node_type*>(children[i]);


                 


            //     *children_cost += child->cost_;
            // }
            // else{
            //      *children_cost += child->cost_;
            // }
            // child_node_data.push_back(temp_child);
            // //static_cast<data_node_type*>(children[i]);


            // auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

            // for (int i = 0; i < temp_child.size(); i++) {
            //     child_values[i].first = temp_child[i];
            //     child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
            // }

            // child->cost_ = data_node_type::compute_new_cost(
            // child_values, temp_child.size(),child->level_, data_node_type::kInitDensity_, const_search, const_traversal,
            // index.params_.expected_insert_frac, &child->model_,
            // index.params_.approximate_cost_computation);


        }
        auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[model_data.size()];

            for (int i = 0; i < model_data.size(); i++) {
                child_values[i].first = model_data[i];
                child_values[i].second = model_payload[i];
            }

        auto new_data_node = new (data_node_type::alloc_type(index.get_allocator()).allocate(1))
            data_node_type(model->level_, index.derived_params_.max_data_node_slots,
                            index.key_comp(), index.get_allocator());

            alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

            new_data_node->bulk_load(child_values, model_data.size(), data_node_model,
                                index.params_.approximate_model_computation);


            // model_nodes[m]->cost_ = data_node_type::compute_new_cost(
            //     child_values, child_node_data.size(),child_node_data.size(), model_nodes[m]->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            //     index.params_.expected_insert_frac, data_node_model,
            //     index.params_.approximate_cost_computation);
            model->cost_ = new_data_node->compute_new_cost(
                child_values, model_data.size(),model_data.size(),model_data.size(),model->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                index.params_.expected_insert_frac, data_node_model,
                index.params_.approximate_cost_computation);

                model->num_keys_model_ = model_data.size();

}

void get_useful_model_nodes_by_level( std::vector<model_node_type*> *model_nodes, std::vector<model_node_type*> *model_nodes_useful_by_level, 
int level){

    (*model_nodes_useful_by_level).clear();
    for(model_node_type* node : *model_nodes){
            //get child nodes and check if they are data nodes

            //Also save models by their level (so I can find the parent much easily)
            if(level == node->level_ ){
                model_nodes_useful_by_level->push_back(node);
            }
       
        }
}

void get_data_nodes_by_level( std::vector<data_node_type*> *data_nodes, std::vector<data_node_type*> *data_nodes_by_level, 
int level){

    //(*model_nodes_useful_by_level).clear();
    for(data_node_type* node : *data_nodes){
            //get child nodes and check if they are data nodes

            //Also save models by their level (so I can find the parent much easily)
            if(level == node->level_ ){
                data_nodes_by_level->push_back(node);
            }
       
        }
}

void poison_data_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal, int *poisoned_count,std::vector<KEY_TYPE>* changed_data, double percentile = 0.99){

    // //(*model_nodes_useful_by_level).clear();
    // std::sort((*data_nodes).begin(), (*data_nodes).end(), compare_searches);

    // for(int i =0; i < (*data_nodes).size(); i++){
    //     std::cout << (*data_nodes)[i]->expected_exp_search_iterations_ << " , ";
    // }

    //double percentile = 0.99;
    

    //int idx = std::distance((*data_nodes).begin(), index);

   // *data_nodes[index]->num_exp_search_iterations_;
   // *(data_nodes)[0]
 
    //(*data_nodes)[2];
    // int idx = static_cast<int>((percentile) * ((*data_nodes).size() - 1));
    // double x = (*data_nodes)[idx]->expected_exp_search_iterations_;

    //std::cout << " Searches  " << x << std::endl;

    std::sort((*data_nodes).begin(), (*data_nodes).end(), compare_searches_des);

    int percent_num = static_cast<int>((*data_nodes).size()*(1-percentile));

    // Vector to store elements over x
    std::vector<data_node_type*> data_nodes_percentile((*data_nodes).begin(), (*data_nodes).begin() + percent_num);
    

    // Use std::copy_if with a lambda function
    // std::copy_if((*data_nodes).begin(), (*data_nodes).end(), std::back_inserter(data_nodes_percentile),
    //              [x](const data_node_type* element) { return element->expected_exp_search_iterations_ >= x; });

    std::cout << " NUM  " << data_nodes_percentile.size() << std::endl;

    for(int i =0; i < (data_nodes_percentile).size(); i++){
        std::cout << (data_nodes_percentile)[i]->expected_exp_search_iterations_ << " , ";
    }

    //Now poison
    int idx2 = 0;
    for(data_node_type* node : (data_nodes_percentile)){

        std::cout << "Poisoning node " << idx2 << std::endl;
        idx2++;

            data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(node,0);

            std::vector<KEY_TYPE> temp_child;
            std::vector<PAYLOAD_TYPE> temp_child_payload;

            while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();

                temp_child.push_back(data_it.key());
                temp_child_payload.push_back(data_it.payload());

                //check if 0 if ok
                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
            }
            std::cout <<"Before " << temp_child.size() << std::endl;
            std::cout << "Data capacity : "<< node->data_capacity_  <<std::endl;
                std::cout << "Keys : "<< node->num_keys_  <<std::endl;
            //std::cout <<"Test1 " << temp_child.size() << std::endl;
            // for(int k = 0 ; k < temp_child.size(); k++){
            //     std::cout << temp_child[k] << " ";
            // }
            //std::cout <<"Test2 " << poi_thres << std::endl;

            if(temp_child.size()>2){
                
                std::vector<KEY_TYPE> poisoned_data = perform_poisoning(temp_child, poi_thres);
                //std::cout <<"Poisoned " << temp_child.size() << std::endl;
                //IF I WANTED TO SAMPLE
                //======================
                // std::vector<KEY_TYPE> new_val = sampleRandom(temp_child,0.5);
                // std::sort(new_val.begin(), new_val.end());
                // std::vector<KEY_TYPE> poisoned_data_temp = perform_poisoning(new_val, poi_thres);

                // std::set<KEY_TYPE> uniqueElements(poisoned_data_temp.begin(), poisoned_data_temp.end());
                // uniqueElements.insert(temp_child.begin(), temp_child.end());

                // // Convert the set back to a vector
                // std::vector<KEY_TYPE> poisoned_data(uniqueElements.begin(), uniqueElements.end());
                // std::sort(poisoned_data.begin(), poisoned_data.end());

                //Need a better condition
                
                if(poisoned_data.size() > temp_child.size()){
                

                //std::vector<int> ranks  ;
                //auto poisoned_values =  get_new_rank_values(temp_child, poisoned_data, ranks, temp_child_payload, temp_child.size());
                //  for(int k = 0 ; k < ranks.size();k++){
                //         std::cout << ranks[k] << " ";
                //     }
                    auto poisoned_values = get_poisoned_values(temp_child,poisoned_data,temp_child_payload,poisoned_data.size());
                    poisoned_count += poisoned_data.size() - temp_child.size();
                     alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

                    // for(int k = 0 ; k < poisoned_data.size();k++){
                    //     std::cout << poisoned_data[k] << " ";
                    // }

                    // for(int k = 0 ; k < temp_child.size();k++){
                    //     std::cout << poisoned_values[k].first << " ";
                    // }
                    

                    //Need to change values
                    //bulk_load this node
                    std::cout << "Poisoned Data : "<< poisoned_data.size()  <<std::endl;

                    

                    node->bulk_load(poisoned_values, poisoned_data.size(), data_node_model,
                                        index.params_.approximate_model_computation);

                    // node->bulk_load_new(poisoned_values, temp_child, poisoned_data.size(), data_node_model,
                    //          index.params_.approximate_model_computation);

                    // node->bulk_load_new_rank(poisoned_values, ranks, temp_child.size(), data_node_model,
                    //                     index.params_.approximate_model_computation);

                    std::cout << "Data capacity : "<< node->data_capacity_  <<std::endl;
                    std::cout << "Keys : "<< node->num_keys_  <<std::endl;

                    // node->bulk_load(poisoned_values, poisoned_data.size(), data_node_model,
                    //                     false);

                    // node->cost_ = data_node_type::compute_new_cost(
                    //     poisoned_values, poisoned_data.size(),temp_child.size(), node->level_, data_node_type::kInitDensity_,const_search, const_traversal,
                    //     index.params_.expected_insert_frac, data_node_model,
                    //     index.params_.approximate_cost_computation);

                    //No need to compute the new costs

                    // node->cost_ = node->compute_new_cost(
                    //     poisoned_values, poisoned_data.size(),temp_child.size(),temp_child.size(), node->level_, data_node_type::kInitDensity_,const_search, const_traversal,
                    //     index.params_.expected_insert_frac, data_node_model,
                    //     index.params_.approximate_cost_computation);

                        for (const auto& element : temp_child) {
                    // Check if the element is not in the original vector
                    //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
                        // If not found, add it to the original vector
                        changed_data->push_back(element);
                    //}
                    }
                }
                    // for (const auto& element : temp_child) {
                    // // Check if the element is not in the original vector
                    // //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
                    //     // If not found, add it to the original vector
                    //     changed_data->push_back(element);
                    // //}
                    // }
             

                
            }
    }

             

        // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
        //         data_it2(node,0);

        //     std::vector<KEY_TYPE> temp_child2;
        //     //std::vector<PAYLOAD_TYPE> temp_child_payload;

        //     while(!data_it2.is_end()){
        //         //nodes[cur_in] = node_it.current();

        //         temp_child2.push_back(data_it2.key());
        //         //temp_child_payload.push_back(data_it.payload());

        //         //check if 0 if ok
        //         data_it2.operator++(0);
        //         // if(node_it.is_end()){
        //         //     break;
        //         // }
        //     }
        //     std::cout <<"After " << temp_child2.size() << std::endl;
        // }
}

void poison_all_data_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal, int *poisoned_count, std::vector<KEY_TYPE>* changed_data){
    std::cout << "number of nodes " << (*data_nodes).size() << std::endl;
    //for(int i =0; i < static_cast<int>((*data_nodes).size()*0.015);i++){
    for(int i =0; i < static_cast<int>((*data_nodes).size()*1);i++){
       // for(int i =0; i < 500;i++){
    //for(data_node_type* node : (*data_nodes)){
        data_node_type* node = (*data_nodes)[i];
        //std::cout << "Poisoning node " << std::endl;

            data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(node,0);

            std::vector<KEY_TYPE> temp_child;
            std::vector<PAYLOAD_TYPE> temp_child_payload;

            while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();

                temp_child.push_back(data_it.key());
                temp_child_payload.push_back(data_it.payload());

                //check if 0 if ok
                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
            }
            if(temp_child.size() > 2){
                
                std::cout << "================="<<std::endl;
                std::cout << "Child size: "<< i <<" : " << temp_child.size() <<std::endl;

                std::cout << "Model a : "<< node->model_.a_  <<std::endl;
                std::cout << "Model b : "<< node->model_.b_  <<std::endl;
                std::cout << "Data capacity : "<< node->data_capacity_  <<std::endl;
                std::cout << "Keys : "<< node->num_keys_  <<std::endl;

                std::vector<KEY_TYPE> poisoned_data = perform_poisoning(temp_child, poi_thres);

                std::pair<uint64_t, double>* poisoned_values;

                //Needs a better condition 
                if(poisoned_data.size() > temp_child.size()){

                    
                    poisoned_values = get_poisoned_values(temp_child,poisoned_data,temp_child_payload,poisoned_data.size());

                    *poisoned_count += poisoned_data.size() - temp_child.size();

                    alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

                    //Need to change values
                    //bulk_load this node
                    node->bulk_load(poisoned_values, poisoned_data.size(), data_node_model,
                                        index.params_.approximate_model_computation);
                    //  node->bulk_load_new(poisoned_values, poisoned_data.size(), data_node_model,
                    //                     index.params_.approximate_model_computation);

                    std::cout << "Model a : "<< node->model_.a_  <<std::endl;
                    std::cout << "Model b : "<< node->model_.b_  <<std::endl;
                    std::cout << "Data capacity : "<< node->data_capacity_  <<std::endl;
                    std::cout << "Keys : "<< node->num_keys_  << std::endl;

                    // node->cost_ = data_node_type::compute_new_cost(
                    //     poisoned_values, poisoned_data.size(),temp_child.size(), node->level_, data_node_type::kInitDensity_,const_search, const_traversal,
                    //     index.params_.expected_insert_frac, data_node_model,
                    //     index.params_.approximate_cost_computation);
                    // node->cost_ = node->compute_new_cost(
                    //     poisoned_values, poisoned_data.size(),temp_child.size(), node->level_, data_node_type::kInitDensity_,const_search, const_traversal,
                    //     index.params_.expected_insert_frac, data_node_model,
                    //     index.params_.approximate_cost_computation);
                    // for (const auto& element : temp_child) {
                    // // Check if the element is not in the original vector
                    // //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
                    //     // If not found, add it to the original vector
                    //     changed_data->push_back(element);
                    // //}
                    // }
                }
                    for (const auto& element : temp_child) {
                    // Check if the element is not in the original vector
                    //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
                        // If not found, add it to the original vector
                        changed_data->push_back(element);
                    //}
                    }
                
                // if(i == 6006){
                //     // if(check_duplicates(temp_child)){
                //     //     std::cout <<  "Duplicates" << std::endl;
                //     // }
                //     // else{
                //     //     std::cout <<  "No Duplicates" << std::endl;
                //     // }
                //     for(int k =0; k < temp_child.size();k++){
                //         std::cout <<  temp_child[k]  << " , ";
                //     }
                //     std::cout  <<std::endl;

                //     for(int k =0; k < poisoned_data.size();k++){
                //         std::cout <<  poisoned_values[k].first  << " , ";
                //        // changed_data->push_back(poisoned_data[k]);
                //     }
                    
                // }
            }
        }
}

bool check_parents_removed_correct(model_node_type* model,data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level){

    bool parent_found = false;
    model_node_type* parent;
    //std::vector<model_node_type*> parents;

    //Start from the level above of the best node
    int level = model->level_ - 1;
    int count =0;

        //For all levels I have to replace the parent
        for(int l = level ; level >=0 ; level--){
            //for(int j =0; j < model_nodes_by_level[l].size(); j++){
            //For each model node in the model nodes by level (all models not just useful)
            for(model_node_type* node : model_nodes_by_level[l]){

                //model_node_type* node = model_nodes_by_level[l][j];
                node_type** childrens  = node->children_;
                //node_type* childrens_pt = *childrens;

                int num_child = node->num_children_;

                for(int i =0; i< num_child; i++){
                    //std::cout << "Child "<< childrens[i] << " , ";
                    if(model == childrens[i]){
                        parent = node;
                        //parents.push_back(node);
                        std::cout<<"parent : " << parent << std::endl;
                        //std::cout << "Parent childern size "<< num_child << std::endl;
                        parent_found = true;
                       
                        count++;
                    }
                }
                
            
            }
        }
        std::cout<< "old model found in : " << count << std::endl;
        
        return parent_found;

}

bool check_parents_added_correct(model_node_type* model,data_node_type *new_data_node,std::vector<model_node_type*> *model_nodes_by_level){

    bool parent_found = false;
    model_node_type* parent;
    //std::vector<model_node_type*> parents;

    //Start from the level above of the best node
    int level = model->level_ - 1;
    int count = 0;

        //For all levels I have to replace the parent
        for(int l = level ; level >=0 ; level--){
            //for(int j =0; j < model_nodes_by_level[l].size(); j++){
            //For each model node in the model nodes by level (all models not just useful)
            for(model_node_type* node : model_nodes_by_level[l]){

                //model_node_type* node = model_nodes_by_level[l][j];
                node_type** childrens  = node->children_;
                //node_type* childrens_pt = *childrens;

                int num_child = node->num_children_;

                for(int i =0; i< num_child; i++){
                    //std::cout << "Child "<< childrens[i] << " , ";
                    if(new_data_node == childrens[i]){
                        parent = node;
                        std::cout<<"parent : " << parent << std::endl;
                        //std::cout << "Parent childern size "<< num_child << std::endl;
                        //parents.push_back(node);
                        parent_found = true;
                       
                        count++;
                        
                    }
                }
               
            
            }
        }
        std::cout<<"new data node found in : " << count << std::endl;
        
        return parent_found;

}



void get_data_by_level(std::vector<data_node_type*> *data_nodes_by_level, std::vector<KEY_TYPE> * data){

    
    for(data_node_type* node : (*data_nodes_by_level)){

        data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
            data_it(node,0);

        std::vector<KEY_TYPE> temp_child;
        std::vector<PAYLOAD_TYPE> temp_child_payload;

        while(!data_it.is_end()){

            (*data).push_back(data_it.key());
    
            data_it.operator++(0);
            
        }
    }
}


void calculate_expected_searches_all_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal){


    for(data_node_type* node : (*data_nodes)){

            data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(node,0);

            std::vector<KEY_TYPE> temp_child;
            std::vector<PAYLOAD_TYPE> temp_child_payload;

            while(!data_it.is_end()){
                //nodes[cur_in] = node_it.current();

                temp_child.push_back(data_it.key());
                temp_child_payload.push_back(data_it.payload());

                //check if 0 if ok
                data_it.operator++(0);
                // if(node_it.is_end()){
                //     break;
                // }
            }
            //std::cout <<"Before " << temp_child.size() << std::endl;
            auto poisoned_values = get_poisoned_values(temp_child,temp_child,temp_child_payload,temp_child.size());



            int data_capacity = std::max(static_cast<int>(temp_child.size() / (data_node_type::kInitDensity_)), static_cast<int>(temp_child.size()) + 1);

            alex::ExpectedIterationsAndShiftsAccumulator acc(data_capacity);
            // alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

            // node->cost_ = node->compute_new_cost(
            //             poisoned_values, temp_child.size(),temp_child.size(),temp_child.size(), node->level_, data_node_type::kInitDensity_,const_search, const_traversal,
            //             index.params_.expected_insert_frac, data_node_model,
            //             index.params_.approximate_cost_computation);

            
            //Calculates the error and these expected ones from alex_base.h
            // if(data_node_model == nullptr){
            //     std::cout <<"Before " << temp_child.size() << std::endl;
            // }
            node->build_node_implicit(poisoned_values, temp_child.size(), data_capacity, &acc, &node->model_);
            
            int expected_avg_exp_search_iterations =
                acc.get_expected_num_search_iterations();
            node->expected_exp_search_iterations_ = expected_avg_exp_search_iterations;

            
            
            
    }

             

        // data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
        //         data_it2(node,0);

        //     std::vector<KEY_TYPE> temp_child2;
        //     //std::vector<PAYLOAD_TYPE> temp_child_payload;

        //     while(!data_it2.is_end()){
        //         //nodes[cur_in] = node_it.current();

        //         temp_child2.push_back(data_it2.key());
        //         //temp_child_payload.push_back(data_it.payload());

        //         //check if 0 if ok
        //         data_it2.operator++(0);
        //         // if(node_it.is_end()){
        //         //     break;
        //         // }
        //     }
        //     std::cout <<"After " << temp_child2.size() << std::endl;
        // }
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