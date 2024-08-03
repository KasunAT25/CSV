// RUN using
// g++ alex_poi.cpp -std=c++17 -o alex_poi -march=native -mpopcnt
// ./alex_poi ../../../ext/Data/fb_200M_uint64 1000000

// ./alex_poi test_genome 1000000 1 0.4 0


// I changed regression_benchmark and fast_brute force io helper
// I changed exponentional search to key_type2 from double
// I removed regression_linear and put it in pgm
//I added number_of_searches to alex_nodes.h and then ++ it into the exponential searches (2) also =0 at their top

#define KEY_TYPE uint64_t
#define PAYLOAD_TYPE double


#include <iostream>
#include <fstream>
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
//#include "src/poisoning_new_min_real_new_long_no_outs.h"
//#include "src/poisoning_new_min_real_max_gap.h"
#include "src/poisoning_simple.h"

//#include "src/poisoning_new_min_real_new.h"

//I can delete the unneccessary ones if needed in pgm benchmark

//#include "src/helpers/alex_benchmark_real_new.h"
#include "src/helpers/alex_benchmark_real_new_ind.h"

//I CHANGED NFL BY ADDING MEAN AND VAR FOR THE NUMERICAL_FLOW THERE. BECUASE IT WAS USING IT WITHOUT CALCULATING IT
// AND GIVING ME A SEGMENTATION ERROR WHEN I USE IT BECAUSE THIS RESULTS IN A 0 FOR ALL.

//NFL check conflicts, numercial_flow and nfl
//Because something is wrong with it, it is not finding the correct key (probably the transformation thing)
//#include "src/helpers/nfl_benchmark_real_new2.h"

#include "src/fast_brute_force_real.h"

#include "src/theil_sen.h"

#include <unordered_set>

//#include "src/helpers/nfl_benchmark.h"

// USE THIS TO RUN THE CODE
// source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include test_alex_no_outs.cpp -std=c++17 -o test_alex_no_outs -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

// source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include benchmark.cpp -std=c++17 -o benchmark -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

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
void shiftVector(std::vector<KEY_TYPE>& vec);
void get_data_by_level(std::vector<data_node_type*> *data_nodes_by_level, std::vector<KEY_TYPE> * data);

void calculate_expected_searches_all_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<data_node_type*> *data_nodes, double poi_thres,
double const_search, double const_traversal);

std::vector<KEY_TYPE> sampleRandom(std::vector<KEY_TYPE> data, double percent);

std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_new_rank_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi, std::vector<int>& ranks,
std::vector<PAYLOAD_TYPE> payload, int size);

bool compare_searches_des(const data_node_type* a, const data_node_type* b);

void shift_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
void shift_back_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);

void get_data_node_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,data_node_type* child, std::vector<KEY_TYPE> *model_node_data);
void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1, std::vector<KEY_TYPE> vec2);

void clearMemoryCache() {
      if (system("sync && sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'") != 0) {
        std::cerr << "Failed to clear memory cache." << std::endl;
    }
}

int main(int argc, char *argv[]){

    clearMemoryCache();

    // bool benchmark = true;
    // std::string data_output = "../../../ext/Data/fb.bin";

    // bool poison = false;
    // bool insert = false;
    // double insert_threshold = 0.5;
    // KEY_TYPE orignal_P = 100000000;
    std::string dataset_name = argv[1];
    std::string poi_size_str = argv[4];
    std::string insert_prop_str = argv[5];
    
    bool benchmark = 1;
    // std::string data_output = "../../../ext/Data/";
    // data_output = data_output + dataset_name+".bin";
    std::string data_folder = "../../../mnt2/Data/";
    //std::string data_output = "../../../mnt2/Data/";
    std::string data_output = data_folder + dataset_name+".bin";

    bool poison = (argv[3][0] == '1');

    

    KEY_TYPE orignal_P = static_cast<KEY_TYPE>(std::stoi(poi_size_str));

    
    double insert_threshold = std::stod(insert_prop_str);
    bool insert = (insert_threshold > 0.0);

    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_output);
    //std::vector<KEY_TYPE> legitimate_data = parse_arguments(argc, argv);

    //shiftVector(legitimate_data);

    // for(int i =0; i< legitimate_data.size();i++){
    //     std::cout<< legitimate_data[i] << " , ";
    // }

    //save_data(legitimate_data, "legitimate_data.txt");
    
    


    //PARAMETERS
    //===============
    
    unsigned seed = 123;
    
    
    double original_poi_thres = 0;
    std::string method = "nt3";
    
    KEY_TYPE P = orignal_P;

    if(poison){
        original_poi_thres = std::stod(poi_size_str);
        method = "ntp3";
    }
    
    //double original_poi_thres = 0;
    //std::string method = "nt";
    double poi_thres = original_poi_thres;
    long pois_count = poi_thres*legitimate_data.size();
    bool use_new_cost = true;
    int num_total_poisoning = 0;
    int max_check = 20;
    
    

    //perform_poisoning(legitimate_data, poi_thres);
    
    srand(12345);
   
   //Check if the dataset has any duplicates
    // if (check_duplicates(legitimate_data)) {
    //     std::cout << "The original data has duplicate values." << std::endl;
    // } else {
    //     std::cout << "The original data does not have duplicate values." << std::endl;
    // }

    //perform_poisoning(legitimate_data, poi_thres);

    //Get the dataname, size and the method to save them
    std::string data_name = argv[1];
    std::string data_name2 = "fb_200M_uint64";
    
    // std::string method = "their_cost";
    // std::string method = "thier_cost_size";
    std::string data_size = argv[2];

    //filenames of the benchmark
    std::string model_output = "results/alex_model_output_real";
    std::string original_output = "results/alex_original_index_full";
    std::string original_changed_output = "results/alex_original_index_changed";
    std::string poisoned_output = "results/alex_poisoned_index_full";
    std::string poisoned_changed_output = "results/alex_poisoned_index_changed";

    original_output = original_output+"_"+method+"_"+data_name2+"_"+data_size;
    original_changed_output = original_changed_output+ "_"+method+"_"+data_name2+"_"+data_size;

    poisoned_output = poisoned_output+"_"+method+"_"+data_name2+"_"+data_size;
    poisoned_changed_output = poisoned_changed_output+ "_"+method+"_"+data_name2+"_"+data_size;

    std::string bulk_index_output = "results/bulk_load_indexes.bin";
    std::string insert_index_output = "results/insert_indexes.bin";

    std::string bulk_data_output = "results/bulk_load_data.bin";
    std::string insert_data_output = "results/insert_data.bin";

    std::string performance_output = "results/alex/alex2_original_performance.csv";
    std::string structure_output = "results/alex/alex2_original_structure.csv";

    if(poison){
        performance_output = "results/alex/alex2_poisoned_performance.csv";
        structure_output = "results/alex/alex2_poisoned_structure.csv";
    }
     std::string changed_data_output = "results/changed_data.bin";
     std::string output_model_csv = "results/alex/alex2_model_output_good.csv";

    //To save the changed data
    std::vector<KEY_TYPE> changed_data;
    std::set<KEY_TYPE> full_data_changed;
    std::set<KEY_TYPE> poisoned_data;

    //To save two Alex indexes for testing (original and edited) I only actually need the index
    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index_original;

    //ADDED THIS HERE
    index.params_.approximate_model_computation = false;

    //BULK LOAD IT (std::pair<KEY_TYPE, PAYLOAD_TYPE>[data.size()], data.size())
    //Create the values (append new payloads)
    std::mt19937_64 gen_payload(std::random_device{}());

    int values_size = legitimate_data.size();
    // auto values = create_values(legitimate_data,&values_size);

    //CAN REMOVE FOR FB
    //std::sort(legitimate_data.begin(), legitimate_data.end());

     //legitimate_data.erase(std::unique(legitimate_data.begin(), legitimate_data.end()), legitimate_data.end());
     //save_data_bin(legitimate_data,data_output);

    values_size = legitimate_data.size();
    std::cout << "removed duplicates." << values_size << std::endl;

    //for inserts
    std::vector<int> all_indexes(values_size);
    std::vector<int> bulk_load_indexes;
    std::vector<int> insert_indexes;

    std::iota(all_indexes.begin(), all_indexes.end(), 0);

    

    //Depending on the insertion 
    int numberOfIndexes = 0;

    auto start_build = std::chrono::high_resolution_clock::now();

    if(insert){
        
        //Instead of this, maybe have a file written.
        std::random_device rd;
        std::mt19937 gen(seed);

        std::shuffle(all_indexes.begin(), all_indexes.end(), gen);

        numberOfIndexes = values_size*(1-insert_threshold);
        bulk_load_indexes.resize(numberOfIndexes);
        insert_indexes.resize(values_size - numberOfIndexes);
       
        bulk_load_indexes.assign(all_indexes.begin(), all_indexes.begin() + numberOfIndexes);
        std::sort(bulk_load_indexes.begin(), bulk_load_indexes.end());

        insert_indexes.assign(all_indexes.begin() + numberOfIndexes, all_indexes.end());
        std::sort(insert_indexes.begin(), insert_indexes.end());

        //SAVING THE INDEXES TO FILES
        // saveIndexesToFile(bulk_load_indexes, bulk_index_output);
        // saveIndexesToFile(insert_indexes, insert_index_output);

        // bulk_load_indexes = readIndexesFromFile(bulk_index_output);
        // insert_indexes = readIndexesFromFile(insert_index_output);

        // for(int i : bulk_load_indexes ){
        //      std::cout << i  << " "; 
        // }

        std::vector<KEY_TYPE> subvector(bulk_load_indexes.size());

        //subvector.resize(bulk_load_indexes.size());
        subvector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(subvector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });

        //SAVING THE DATA TO FILES
        // save_data_bin(subvector, bulk_data_output);
        
        // subvector.clear();
        
        // subvector = read_data_bin(bulk_data_output);


        // for(int i : subvector ){
        //      std::cout << i  << " "; 
        // }
        std::cout << "subvector "  <<subvector.size() << std::endl;

        std::cout << "Right"  << std::endl;

        auto values = create_values(subvector,&numberOfIndexes);
       
        //std::pair<KEY_TYPE, PAYLOAD_TYPE> *values2 = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[numberOfIndexes];

        // std::copy_if(values, values + values_size, values2,
        //          [&](const std::pair<uint64_t, double>& value) {
        //              return std::find(bulk_load_indexes.begin(), bulk_load_indexes.end(),
        //                               static_cast<KEY_TYPE>(value.first)) != bulk_load_indexes.end();
        //          });

        //index_original.bulk_load(values, numberOfIndexes);
        index.bulk_load(values, numberOfIndexes);

        // index_original.bulk_load(values, values_size);
        // index.bulk_load(values, values_size);

    }
    else{
        std::cout << "Creating values" << std::endl;
        auto values = create_values(legitimate_data, &values_size);
        std::cout << "Created values" << std::endl;

        //index_original.bulk_load(values, values_size);
        index.bulk_load(values, values_size);
        std::cout << "Created index" << std::endl;
    }
     std::vector<int>().swap(all_indexes);
    std::vector<int>().swap(bulk_load_indexes);

    //auto values = create_values(legitimate_data,legitimate_data.size());

    // std::cout << "Created values" << std::endl;
    std::cout << "============================" << std::endl;

    //bulk_load indexes and measure their times per index build
    

    // index_original.bulk_load(values, legitimate_data.size());
    // index.bulk_load(values, legitimate_data.size());
    

    auto stop_build = std::chrono::high_resolution_clock::now();

    // std::cout << "Created the index" << std::endl;
    // std::cout << "============================" << std::endl;

    long index_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_build - start_build).count();
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;

    //Calculating the new constants by using the created original index
    auto start_tra = std::chrono::high_resolution_clock::now();

    double const_search = 0;
    double const_traversal =  0;

    //Im getting the neccessary constants from here
    get_constants(index,legitimate_data, 1000000, &const_search, &const_traversal);

    std::cout << "const_search " << const_search <<std::endl;
    std::cout << "const_traversal " << const_traversal <<std::endl;

   //Getting all the nodes from the index and the maximum height (default is 0 aka root node only)
    std::vector<node_type*> nodes;
    std::vector<model_node_type*> model_nodes;
    std::vector<data_node_type*> data_nodes;

    int max_model_height =0;
    
    //Getting all nodes and seperates them by model nodes and data nodes also get the max model height
    get_all_nodes( index, &nodes, &model_nodes,&data_nodes,&max_model_height);

    int original_model_nodes_number = model_nodes.size();
    int original_nodes_number = nodes.size();
    int original_data_nodes_number = data_nodes.size();

    std::cout << "Done: Gathered the nodes" << std::endl;
    std::cout << "============================" << std::endl;
    //std::cout << "Max model height " << max_model_height<< std::endl;

    //For all data nodes
    // if(use_new_cost){
    //     //Calculate the new cost
    //     //calculate_data_node_new_cost(data_nodes, index,const_search,const_traversal);    
    // }

    // Add a different condition OR a different method
    //======================
    int num_pois = model_nodes.size() -1;

    //num_pois = 0;
    //num_pois = max_model_height;

    int replaced_models_number = 0;

    double pois_count_model = (pois_count/(1.0*num_pois));

    //Method 3 & 4
    //calculate_data_node_new_cost(data_nodes, index); 
    //Get the model nodes with all the children being data nodes.
    std::vector<model_node_type*> model_nodes_useful;
    //std::vector<double> children_sum;

    // std::vector<model_node_type*> model_nodes_by_level[max_model_height+1];
    std::vector<model_node_type*> *model_nodes_by_level = new std::vector<model_node_type*>[max_model_height + 1];

    //Get the model nodes by level and the useful (all children are data nodes) models
    //I also calculate the children sum for these . (CHILDREN SUM IS  NO LONGER NEEDED)
    get_useful_model_nodes(&model_nodes, &model_nodes_useful, model_nodes_by_level);
    //Get the model nodes with all the children being data nodes.
    std::cout << "DONE: models useful and models by level" << std::endl;
    //std::cout << "models useful " << model_nodes_useful.size() <<std::endl;
    std::cout << "============================" << std::endl;

    if(use_new_cost){
        //For useful models calculate the new cost
        //Calculates the new model costs
        calculate_model_node_new_cost(model_nodes_useful,index,const_search,const_traversal); 
        //Calculates the new children costs as well as the sum of children costs per model
        calculate_data_node_new_cost_new(model_nodes_useful, index,const_search,const_traversal);     
     }

    

    //calculate_model_node_poisoning(model_nodes_useful,index,pois_count_model, const_search,const_traversal);
    bool no_improvement = false;

    //Get the difference between the model node and the children sum (now new cost)
    calculate_cost_difference(model_nodes_useful);

    //std::sort(model_nodes_useful.begin(), model_nodes_useful.end(), compare_models);
    std::vector<model_node_type*> model_nodes_useful_by_level;
    
    //Get the model nodes useful in the max model height
    int current_level = max_model_height;

    get_useful_model_nodes_by_level(&model_nodes_useful, &model_nodes_useful_by_level, current_level);

    
    std::cout << "Current_level " << current_level << std::endl;

    //Do for all the levels except the root (keep the root as it is)
    //Inner index to iterate through the models in the current level (reset in the outter while loop)
    int inner_idx = 0;
    //To keep count of the number of models we have tried without use
    int tried = 0;
    int sucess =  0 ;

    std::vector<model_node_type*> best_nodes;
    std::vector<data_node_type*> replaced_data_nodes;

    std::vector<model_node_type*> new_useful_parents;
    std::vector<KEY_TYPE>  altered_data;
  
    bool continue_poi = true;

    while(current_level > 0 && poison && continue_poi){

        long model_level_size = 0;

        for(int i =0; i< model_nodes_useful_by_level.size();i++){
            model_level_size += model_nodes_useful_by_level[i]->num_keys_model_;
        }

        if(model_level_size > 100000000){
            continue_poi = false;
            break;
        }
    //while(current_level > 0){
        // std::cout << "while outter " << current_level << std::endl;
        // std::cout << "Models Number " << model_nodes_useful_by_level.size() << std::endl;

        

    //for(int poisoing = num_pois ; poisoing > 0; poisoing--){
    //for(int poisoing = 0 ; poisoing < num_pois; poisoing++){

        // if(no_improvement){
        //     break;
        // }

        //Get the best nodes by some condition
        //int cur_in3 = 0;
        model_node_type* best_node;
        //double cost_diff;
        //double temp_cost_diff;

        //std::vector<std::pair<model_node_type*,double>> cost_differences;

        // best_node = calculate_best_node_cost_prev(model_nodes_useful, children_sum);
        //METHOD 1: lowest difference
        //cost_differences = calculate_best_node_cost_prev(model_nodes_useful);
        // std::sort(cost_differences.begin(), cost_differences.end(), 
        // [](const auto& left, const auto& right) {
        //     return left.second < right.second;
        // } );

        //METHOD 2 before and after poisoning cost reduction (does not take children into account)
        // std::vector<double> before = calculate_model_node_new_cost(model_nodes_useful,index);
        // std::vector<double> afrer =calculate_model_node_poisoning(model_nodes_useful,index,poi_thres);

        // cost_differences = calculate_best_node_cost_reduction(model_nodes_useful,before,afrer);

        //METHOD 3 poisoning cost reduction (does take children into account)
        // std::vector<double> afrer =calculate_model_node_poisoning(model_nodes_useful,index,poi_thres);
        
        // cost_differences = calculate_best_node_cost_prev(model_nodes_useful, children_sum);

        //METHOD 4: lowest difference new costs
        // calculate_model_node_new_cost(model_nodes_useful,index);
        // cost_differences = calculate_best_node_cost_prev(model_nodes_useful, children_sum);
        
        //METHOD 5:
        //Calculates the difference
        //Then sort thme

        //calculate_cost_difference(model_nodes_useful);

        //Sort the models in this level by their cost difference
        std::sort(model_nodes_useful_by_level.begin(), model_nodes_useful_by_level.end(), compare_models_new);

        // for(int i =0; i < (model_nodes_useful_by_level).size()/2; i++){
        // std::cout << (model_nodes_useful_by_level)[i]->cost_diff_ << " , ";
        // }

        // for(int i =0; i< model_nodes_useful_by_level.size();i++){
        //     std::cout << model_nodes_useful_by_level[i]->cost_diff_ << " , ";
        // }

        //METHOD 6
        //Calculates the model_node_poisoning (aka the cost after the poisoning)
        // calculate_model_node_poisoning(model_nodes_useful,index,pois_count_model, const_search,const_traversal);
        //Calculates the difference
        // calculate_cost_difference(model_nodes_useful);
        // std::vector<model_node_type*> model_reduce;
        // model_reduce = can_be_done(model_nodes_useful);
        // std::sort(model_nodes_useful.begin(), model_nodes_useful.end(), compare_models);

        //METHOD OLD
        //std::sort(model_nodes_useful.begin(), model_nodes_useful.end(), compare_models_old);
        //sort it by the cost lowest
        

        
        //To store if the model has converted for the current index
        bool model_converted = false;
        inner_idx = 0;

        
        //double check this while loop it may be in a infinite loop

        // do until the model is not converted or more than 10 percent of the models are checked
        //METHOD 1-5
        int size_to_check = model_nodes_useful_by_level.size();
        //METHOD 6
        //int size_to_check = model_reduce.size();
        

        // while(!model_converted && inner_idx < std::max(static_cast<int>(size_to_check*0.1),1)){
        // while(!model_converted && inner_idx < size_to_check){
        tried = 0;
        sucess =  0; 

        //inner while to loop through the inner index until we reach the max models useful in the level
        //Stop if we have reached the max_check limit without changind a model to a data node.
        

        while(inner_idx < size_to_check && tried < max_check){

            if(inner_idx%1000 ==0){
                std::cout << "starting " << inner_idx<<std::endl;
            }

            // if(sucess > 200){
            //     break;
            // }

            // for(int i = 0 ; i < model_nodes_useful_by_level.size() ; i++){
            //     std::cout << model_nodes_useful_by_level[i] << " | " ;
            // }
        //std::cout << std::endl;

            //IF we have tried 
            // if(tried > max_check){
            //     break;
            // }
         //std::cout << "while inner " << inner_idx << std::endl;
        //while(!model_converted && inner_idx < 1){
          //  std::cout << "WHILE " << inner_idx << std::endl;
        //while(inner_idx == 0){

            //METHOD 1-4:
            // best_node = cost_differences[inner_idx].first;
            //METHOD 5:
            //Get the best node as the model node in that level that has the lowest cost difference
            best_node = model_nodes_useful_by_level[inner_idx];

            //METHOD 6
            //best_node = model_reduce[inner_idx];

            //std::cout << "DONE: best model selection " << best_node <<std::endl;
            // std::cout << "============================" << std::endl;


            //std::cout << "Best Model cost: " <<model_nodes_useful.size()<< std::endl;
            //std::cout << "Best Model cost: " <<children_sum.size()<< std::endl;
            //std::cout << "Best Model cost: " <<best_node->cost_<< std::endl;

            //Now get the data for this model node
            double best_children_cost = 0;
            std::vector<KEY_TYPE> model_node_data;
            std::vector<PAYLOAD_TYPE> model_node_payload;

            //Get the children data (aka all the model data) for the best node
            //We don't need to recalculate the children cost because we have already calculated it before (useful models method last para false here)
            //Also get the payloads
            get_children_data(index, best_node, &model_node_data,&model_node_payload, &best_children_cost, const_search, const_traversal, false);

            // std::cout << "DONE: model & children node data" << std::endl;
            // std::cout << "============================" << std::endl;


            //Poison using the model data or till the cost is less than the children nodes sum of costs.
            //Vector to hold the poisoned model data or if no poisoning use this as well. I put these outside the if becasue I need to use it in else as well
            std::vector<KEY_TYPE> poisoned_model_data;
            //poi_thres = (pois_count_model/(1.0*model_node_data.size()));
            poi_thres = original_poi_thres;
            //poi_thres = original_poi_thres*(max_model_height-current_level+1);
            
            int poisoned_model_data_size;
            data_node_type* new_data_node;
            int current_poi =0;

            // if(best_node->cost_diff_ > 0){
            //     inner_idx++;
            //     tried++;
            //     break;
            // }

            
            sucess++;
            int cutoff = -50 ;
            //best_node->cost_diff_ = best_node->cost_diff_/(1.0*model_node_data.size());
            //If the cost differebce is greater than 0 meaning model cost > children cost try poisoning
            if(best_node->cost_diff_ > (cutoff) || true){
                
                
                //std::cout << "while inner POISONING " << inner_idx << std::endl;

                KEY_TYPE first_val = model_node_data[0];
                shift_vector(&model_node_data, first_val);

                //poi_thres = 1/(1.0*model_node_data.size());
                //Do the poisoning for the model data
                poisoned_model_data = perform_poisoning(model_node_data, poi_thres);

                shift_back_vector(&poisoned_model_data, first_val);
                shift_back_vector(&model_node_data, first_val);

                //std::cout << "DONE: Poisoning" << std::endl;
               
                // std::cout << "============================" << std::endl;
                // if(inner_idx == 5043){
                //     std::cout<< model_node_data.size() << std::endl;
                    //  for (int i = 0; i < model_node_data.size(); i++) {
                    // //5562 6249
                    // std::cout<< model_node_data[i] << " , ";
                    // }
               // }

                // if(inner_idx == 5044){
                //      std::cout<< model_node_data.size() << std::endl;
                //      for (int i = 0; i < model_node_data.size(); i++) {
                //     //5562 6249
                //     std::cout<< model_node_data[i] << " , ";
                //     }
                //     exit(0);
                // }
               


                //std::cout << "model data "<< model_node_data.size() << std::endl;

                //int poisoned_model_data_size = poisoned_model_data.size();
                poisoned_model_data_size = poisoned_model_data.size();

                //std::cout << "poisoned model data "<< poisoned_model_data_size << std::endl;

                //get new values for the poisoning data while using the old ones for the existing
                auto poisoned_values = get_poisoned_values(model_node_data, poisoned_model_data, model_node_payload, poisoned_model_data_size);
                //std::cout << "DONE: Getting values" << std::endl;

                // auto poisoned_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[poisoned_model_data.size()];
                //     //std::mt19937_64 gen_payload(std::random_device{}());

                //         for (int i = 0; i < poisoned_model_data.size(); i++) {
                //             poisoned_values[i].first = model_node_data[i];
                //             poisoned_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                //             }

                // std::cout << "DONE: poisoned data values" << std::endl;
                // std::cout << "============================" << std::endl;

                //Create a new data node using these poisoned values

               //ALso check the max data node number here itself
               //Add a constant that if there were no new points also save
                 if(index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size ){
                    inner_idx++;
                    tried++;
                    //replace_parents_data_node(new_data_node,model_nodes_by_level, &best_node_idx, &parent_model_idx);
                    //std::cout << "SKIP" << std::endl;
                    continue;
                }
                // std::cout << "Creating new data" << std::endl;
                // std::cout << "Keys : "<< poisoned_model_data_size  <<std::endl;

                new_data_node = create_new_data_node(index , best_node, poisoned_values,model_node_data, poisoned_model_data_size, model_node_data.size(),
                const_search, const_traversal);
                current_poi = poisoned_model_data_size - model_node_data.size();
                //std::cout << "==========================" << inner_idx << std::endl;
            }
            //Can directly create a new data node instead of a model
            else{
                
                //std::cout << "while inner " << inner_idx << std::endl;
                //poisoned_model_data = perform_poisoning(model_node_data, poi_thres);

                //std::cout << "DONE: No Poisoning" << std::endl;
                // std::cout << "============================" << std::endl;


                //std::cout << "model data "<< model_node_data.size() << std::endl;

                //int poisoned_model_data_size = poisoned_model_data.size();
                //poisoned_model_data_size = poisoned_model_data.size();

                //std::cout << "poisoned model data "<< poisoned_model_data_size << std::endl;

                //get new values for the poisoning data while using the old ones for the existing
                //auto poisoned_values = get_poisoned_values(model_node_data, poisoned_model_data, model_node_payload, poisoned_model_data_size);

                    // std::cout << "DONE: poisoned data values" << std::endl;
                    // std::cout << "============================" << std::endl;

                    auto values_real = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[model_node_data.size()];
                    //std::mt19937_64 gen_payload(std::random_device{}());

                        for (int i = 0; i < model_node_data.size(); i++) {
                            values_real[i].first = model_node_data[i];
                            values_real[i].second = model_node_payload[i];
                            }
                

                //Create a new data node using these values
                //Also check if index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ before creating the new data nodes
                //std::cout << index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ << std::endl;
                //std::cout << model_node_data.size()<< std::endl;
                if(index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= model_node_data.size()){
                    inner_idx++;
                    tried++;
                    //replace_parents_data_node(new_data_node,model_nodes_by_level, &best_node_idx, &parent_model_idx);
                    //std::cout << "SKIP" << std::endl;
                    continue;
                }
                
                new_data_node = create_new_data_node(index , best_node, values_real,model_node_data, model_node_data.size(), model_node_data.size(),
                const_search, const_traversal);
                

                poisoned_model_data_size = model_node_data.size();
            }

            // for (const auto& element : model_node_data) {
            //     // Check if the element is not in the original vector
            //     //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
            //         // If not found, add it to the original vector
            //         changed_data.push_back(element);
            //     //}
            // }
            // poisoned_model_data = perform_poisoning(model_node_data, poi_thres);

             //std::cout << "DONE: Poisoning " << new_data_node << std::endl;
            // std::cout << "============================" << std::endl;


            // //std::cout << "model data "<< model_node_data.size() << std::endl;

            // //int poisoned_model_data_size = poisoned_model_data.size();
            // poisoned_model_data_size = poisoned_model_data.size();

            // //std::cout << "poisoned model data "<< poisoned_model_data_size << std::endl;

            // //get new values for the poisoning data while using the old ones for the existing
            // auto poisoned_values = get_poisoned_values(model_node_data, poisoned_model_data, model_node_payload, poisoned_model_data_size);

            // std::cout << "DONE: poisoned data values" << std::endl;
            // std::cout << "============================" << std::endl;

            // //Create a new data node using these values
            // auto new_data_node = create_new_data_node(index , best_node, poisoned_values, poisoned_model_data_size, model_node_data.size(),
            // const_search, const_traversal);

            // std::cout << "DONE: New data node cost calculated" << std::endl;
            // std::cout << "============================" << std::endl;
            // std::cout << "Before data cost node "<< best_node->cost_ << std::endl;
            //std::cout << "poisoned data cost node "<< new_data_node->cost_ << std::endl;
            //std::cout << "Children cost node "<< best_children_cost << std::endl;

            // std::cout << "num_nodes "<< nodes.size() << "," << index.num_nodes() << std::endl;
            // std::cout << "model nodes number "<< model_nodes.size() << std::endl;
            // std::cout << "useful model nodes number "<< model_nodes_useful.size() << std::endl;
            //std::cout << "data nodes number "<< data_nodes.size() << std::endl;

            //ADD some condition to change the nodes
            int best_level = best_node->level_;
            int best_node_idx = 0;
            bool parent_found = false;
            int parent_model_idx = 0;

            //If the cost of the new data is higher than the children or there is no more space then try the next best node
            //Change the condition here
            // if(new_data_node->cost_ > best_node->children_cost || index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size){
            //std::cout << "Cost : "<< (new_data_node->cost_ - best_node->children_cost)  <<std::endl;

            if((new_data_node->cost_ - best_node->children_cost) > (cutoff) || index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size){
            // if(new_data_node->cost_ > best_node->children_cost || index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size){
                inner_idx++;
                tried++;
                //replace_parents_data_node(new_data_node,model_nodes_by_level, &best_node_idx, &parent_model_idx);
                //std::cout << "SKIP" << std::endl;
                continue;
            }

            double best_children_cost_new = 0;
            std::vector<KEY_TYPE> model_node_data_new;
            std::vector<PAYLOAD_TYPE> model_node_payload_new;

            //Get the children data (aka all the model data) for the best node
            //We don't need to recalculate the children cost because we have already calculated it before (useful models method last para false here)
            //Also get the payloads
            std::vector<KEY_TYPE>  poi_data;
            get_data_node_data(index, new_data_node, &model_node_data_new);
            find_difference(&poi_data,model_node_data, model_node_data_new);

            for(const auto& element : poi_data){
                poisoned_data.insert(element);
            }
                //if(current_level == max_model_height){
                    //std::cout << "Adding"  <<std::endl;
                    for (KEY_TYPE element : model_node_data) {
                    // Check if the element is not in the original vector
                    //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
                        // If not found, add it to the original vector
                        //changed_data.push_back(element);
                        full_data_changed.insert(element);
                        if (poisoned_data.find(element) == poisoned_data.end()) {
                            // Add the element to the result vector
                            altered_data.push_back(element);
                        }
                    //}
                
                    }
                //}
                
            
            
            tried = 0;

            //this means we are going to replace the model node with the new data node 
                
            //std::cout << "Best model height "<< best_level << std::endl;
            // std::cout << "Best model  "<< best_node << std::endl;

            //If the best level is not 0 aka it is not the root model.
            //if(best_level != 0){
                
            // int best_node_idx = 0;
            // bool parent_found = false;
            // int parent_model_idx = 0;

            //std::cout << "model nodes by level " << (*model_nodes_by_level).size() << std::endl;

            //find the parent of this new data node aka the old best node
            std::vector<model_node_type*> parents;
            //parents = find_parent(best_node, model_nodes_by_level, &best_node_idx, &parent_model_idx);
                
            
            // std::cout << "DONE: Finding parents" << std::endl;
            // std::cout << "============================" << std::endl;
            // std::cout << parent << std::endl;

            
            new_data_node = static_cast<data_node_type*>(new_data_node);
            
            // std::cout << "new data node to be replaced "<< new_data_node << std::endl;
            // std::cout << "Parent child by level " << model_nodes_by_level[best_level-1][parent_model_idx]->children_[best_node_idx] << std::endl;
            // std::cout << "Parent child " << parent->children_[best_node_idx] << std::endl;

            //Change the parent's child to this new data node pointer
            //Issue is here
            //parent->children_[best_node_idx] = new_data_node;
            //check_parents_added_correct(best_node,new_data_node,model_nodes_by_level);
            //check_parents_removed_correct(best_node,new_data_node,model_nodes_by_level);
            parents = replace_parents(best_node, new_data_node, model_nodes_by_level, &best_node_idx, &parent_model_idx);
            //check_parents_added_correct(best_node,new_data_node,model_nodes_by_level);
            //check_parents_removed_correct(best_node,new_data_node,model_nodes_by_level);
            //std::cout << "Parent child " << parent->children_[best_node_idx] << std::endl;

            //IF I REALLY NEEDED TO CHECK, THEN I COULD GET THE NODES FROM THE INDEX AGAIN AND CHECK (ALL PREV NODES WILL COME, BUT THERE WILL BE NO LINK)]
            // TO DO
            // ADD A NEW COST MODEL IN alex_nodes.h NEED TO ACCOUNT FOR THE DEPTH AS WELL DEPTH VS.
            // ALSO NEED TO MAKE SURE ALL THE DATA CAN BE FITTED IN THE NEW DATA NODE (SIMPLE IF STATEMENT WITH MAX_NODE SIZE AFTER FINDING BEST NODE)
            // CALCULATE IT INSTEAD OF THE COST FOR BOTH MODEL AND CHILDREN
            // CHECK IF ITS BETTER.
            //  IF SO DO THE REST OF IT 

            //Removing from the vectors 
            // std::cout << "a1 "<< model_nodes.size() << std::endl;
            // std::cout << "b1 "<< model_nodes_by_level[best_level].size() << std::endl;
            // std::cout << "c1 "<< model_nodes_useful.size() << std::endl;

            //find and delete from model_nodes, model_nodes_by_level and model_nodes_useful also check if parent is
            // a useful model node (all children are data nodes) then add that to the model_nodes_useful
            (new_useful_parents).clear();
            update_data_structure(best_node, &model_nodes, model_nodes_by_level, &model_nodes_useful,&model_nodes_useful_by_level, parents,
            index,const_search,const_traversal,&new_useful_parents);

            //std::cout << "data node "<< new_data_node <<std::endl;
            //check_parents_removed_correct(best_node,new_data_node,model_nodes_by_level);
            //check_parents_added_correct(best_node,new_data_node,model_nodes_by_level);

            replaced_models_number++;

            // std::cout << "a2 "<< model_nodes.size() << std::endl;
            // std::cout << "b2 "<< model_nodes_by_level[best_level].size() << std::endl;
            // std::cout << "c2 "<< model_nodes_useful.size() << std::endl;

            // Sort the data vector
            // std::sort(changed_data.begin(), changed_data.end());

            // Iterate through the additional vector
            //Add the changed data here
            // for (const auto& element : model_node_data) {
            //     // Check if the element is not in the original vector
            //     //if (std::find(changed_data.begin(), changed_data.end(), element) == changed_data.end()) {
            //         // If not found, add it to the original vector
            //         changed_data.push_back(element);
            //     //}
            // }
            model_converted = true;

           // }
            
            // std::cout << "model nodes number "<< model_nodes.size() << std::endl;
            // std::cout << "model nodes number "<< model_nodes_by_level[best_level].size() << std::endl;
            // std::cout << "Replaced number " << replaced_models_number <<std::endl;

            // if(inner_idx == std::max(static_cast<int>(size_to_check*0.1),1)-1){
            //     no_improvement = true;
            // }
        
            //std::cout << "Done: Iteration " << poisoing <<std::endl;
            
            num_total_poisoning += current_poi;
            //If model is converted then set tried back to 0 and increase the inner_idx
            if(model_converted){
                inner_idx++;
                tried = 0;
            }
            // if(!model_converted){
            //     tried++;
            // }
            // else{
            //     inner_idx++;
            //     tried = 0;
            // }

            // if(tried > 10){
            //     break;
            // }

        }
        
        // if(!model_converted){
        //     current_level--;
        // }

        //After a level go to the upper level
        std::cout << "=============="  <<std::endl;
        std::cout << "Done: Iteration " << current_level <<std::endl;

        current_level--;
        if(current_level  > 0){
            //Get the new models by level from model nodes useful (which we update in the update_data_structures
            get_useful_model_nodes_by_level(&model_nodes_useful, &model_nodes_useful_by_level, current_level);

             //(model_nodes_useful_by_level).clear();
             //std::cout << "Done: Iteration " << new_useful_parents.size() <<std::endl;
             //model_nodes_useful_by_level = new_useful_parents;
        }
        //std::cout << "==============" << poisoing <<std::endl;
        //std::cout << "Done: Iteration " << poisoing <<std::endl;
        
        //std::cout << "Done: Iteration " << new_useful_parents.size() <<std::endl;
        //std::cout << "New useful " << changed_data.size() <<std::endl;
        

        //Testing

        // std::vector<node_type*> nodes_fin;
        // std::vector<model_node_type*> model_nodes_fin;
        // std::vector<data_node_type*> data_nodes_fin;

        // //Getting all the nodes from the index and the maximum height (default is 0 aka root node only)
        // //Printing information at the end
        // int max_model_height_fin = 0;
        

        // get_all_nodes(index, &nodes_fin, &model_nodes_fin,&data_nodes_fin, &max_model_height_fin);

        // if(model_nodes.size() != model_nodes_fin.size()){
        //     std::cout << "BREAKING" << std::endl;
        //     std::cout << best_node << std::endl;
        //     //break;
        // }
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
        // std::cout << "Replaced number " << replaced_models_number <<std::endl;
            
    }

    if(insert){
   
        for(int i : insert_indexes){
            //index_original.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
            index.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
        }
    }
    //testing
    // get_useful_model_nodes_by_level(&model_nodes, &model_nodes_useful_by_level, 0);

    // for(int i =0; i < model_nodes_useful_by_level.size();i++){
    //     std::cout << model_nodes_useful_by_level[i] << std::endl;
    //     //print all children
    //     node_type** children = model_nodes_useful_by_level[i]->children_;
    //     //std::vector< std::vector<KEY_TYPE> > child_node_data;
        
    //     int num_child = model_nodes_useful_by_level[i]->num_children_;
    //     std:: cout << "Children" << std::endl;
    //     for(int i = 0; i < num_child;i++){
    //         data_node_type* child = static_cast<data_node_type*>(children[i]);
    //        // std:: cout << child << " | ";
    //     }
    //     std:: cout << std::endl;

    // }

    // auto stop_tra = std::chrono::high_resolution_clock::now();
   
    
    // if (check_duplicates(changed_data)) {
    //     std::cout << "The changed data has duplicate values." << std::endl;
    // } else {
    //     std::cout << "The changed data does not have duplicate values." << std::endl;
    // }
    std::sort(changed_data.begin(), changed_data.end());
    std::vector<node_type*> nodes_fin;
    std::vector<model_node_type*> model_nodes_fin;
    std::vector<data_node_type*> data_nodes_fin;

    std::vector<model_node_type*> model_nodes_useful_fin;

    //Getting all the nodes from the index and the maximum height (default is 0 aka root node only)
    //Printing information at the end
    int max_model_height_fin = 0;
    

    get_all_nodes(index, &nodes_fin, &model_nodes_fin,&data_nodes_fin, &max_model_height_fin);
    

    // std::vector<data_node_type*> data_nodes_level;
    // get_data_nodes_by_level(&data_nodes_fin,&data_nodes_level,max_model_height+1);
    //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);

    //poison_data_nodes(index,&data_nodes_fin, poi_thres, const_search, const_traversal,0.9);
    //int *poisoned_count;
    

    //  std::cout << std::endl;
    // std::sort((data_nodes_fin).begin(), (data_nodes_fin).end(), compare_searches);

    // for(int i =0; i < (data_nodes_fin).size(); i++){
    //     std::cout << (data_nodes_fin)[i]->expected_exp_search_iterations_ << " , ";
    // }
    //poison_all_data_nodes(index,&data_nodes_fin, poi_thres, const_search, const_traversal);

    auto stop_tra = std::chrono::high_resolution_clock::now();

    //Original information
    std::cout << "====================" << std::endl;

    std::cout << "original nodes number "<< original_nodes_number << std::endl;
    std::cout << "original model nodes number "<< original_model_nodes_number << std::endl;
    std::cout << "original data nodes number "<< original_data_nodes_number << std::endl;

    std::cout << "====================" << std::endl;

    //afterwards
    std::cout << "nodes number "<< nodes_fin.size() << std::endl;
    std::cout << "model nodes number "<< model_nodes_fin.size() << std::endl;
    std::cout << "data nodes number "<< data_nodes_fin.size() << std::endl;

    std::cout << "====================" << std::endl;

    std::cout << "model_nodes : "<< model_nodes.size() << std::endl;
    std::cout << "model_nodes_useful " << model_nodes_useful.size() << std::endl;
    std::cout << "Replaced number " << replaced_models_number <<std::endl;

    
    //CHECK WHY THERE IS A DIFFERENCE
    // std::cout << "model nodes number "<< model_nodes_fin.size() << std::endl;
    // std::cout << "model nodes number "<< model_nodes.size() << std::endl;

    // for(int i =0; i<model_nodes_fin.size(); i++ ){
    //     std::cout << model_nodes_fin[i]<< " , ";
    // }
    //  std::cout << std::endl;

    //   for(int i =0; i<model_nodes.size(); i++ ){
    //     std::cout << model_nodes[i]<< " , ";
    // }
    //  std::cout << std::endl;
    
    std::cout << "Poisoned number " << num_total_poisoning <<std::endl;
    

    // std::cout << "Constants" <<std::endl;
    // std::cout << const_search <<std::endl;
    // std::cout << const_traversal <<std::endl;

    //Perform the benchmarks
    //std::sort(changed_data.begin(), changed_data.end());
    //For ALEX
    //===========
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

    size_t total_data_count = 0;
    size_t total_node_count = 0;

    std::string structure_results = "ALEX;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning)+";"+ (poi_size_str) +";"
        + std::to_string(altered_data.size())+ ";Data";

    std::vector<model_node_type*> *model_nodes_by_level_fin = new std::vector<model_node_type*>[max_model_height_fin + 1];
    get_useful_model_nodes(&model_nodes_fin, &model_nodes_useful_fin, model_nodes_by_level_fin);

    //TO BENCHMARK LEVELS
    double loss_value = 0.0;
    double mse = 0.0;
    std::string data_node_info = "Data Nodes";
    std::string model_node_info = "Model Nodes";


    for(int i = 1; i <= max_model_height+1 ; i++ ){
        std::cout << "Performance For Level  " << i << std::endl;
        std::vector<data_node_type*> data_nodes_level;
        std::vector<model_node_type*> model_nodes_level;
        get_data_nodes_by_level(&data_nodes_fin,&data_nodes_level,i);
        std::cout << "Data nodes number  " << data_nodes_level.size() << std::endl;

        
        // get_useful_model_nodes_by_level(&model_nodes_fin,&model_nodes_level,i);
        get_useful_model_nodes_by_level(&model_nodes_useful_fin,&model_nodes_level,i);
        std::cout << "model nodes number  " << model_nodes_level.size() << std::endl;
        
        data_node_info  = data_node_info + ";" + std::to_string(data_nodes_level.size());
        model_node_info  = model_node_info + ";" + std::to_string(model_nodes_level.size());
        
        std::vector<KEY_TYPE> level_data;
        get_data_by_level(&data_nodes_level, &level_data);

        std::cout << "data number  " << level_data.size() << std::endl;
        if(i>0){
            total_data_count += level_data.size();
            total_node_count += data_nodes_level.size() + model_nodes_level.size();
        }

        structure_results = structure_results + ";" + std::to_string(level_data.size());
        // benchmark_alex_real_seperate(index_original, level_data, level_data, "ALEX", data_name, poi_thres, original_changed_output);
        //poison_all_data_nodes(index,&data_nodes_level, poi_thres, const_search, const_traversal, &num_total_poisoning, &changed_data);

        for(data_node_type* node : data_nodes_level){
            std::vector<KEY_TYPE> x;
            get_data_node_data(index, node, &x);

            int n = x.size();
            if(n <= 1){
                continue;
            }
            shift_vector(&x, x[0]);

            std::vector<int> y(n);
            std::iota(y.begin(), y.end(), 0);
            long double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXSquare = 0.0;

            for(int i =0; i < n;i++){
                sumX += static_cast<long double>(x[i]);
                sumY += y[i];
                sumXY += static_cast<long double>(x[i]) * y[i];
                sumXSquare += static_cast<long double>(x[i]) * static_cast<long double>(x[i]);
            }

            double numerator = n * sumXY - sumX * sumY;
            double denominator = n * sumXSquare - sumX * sumX;

            double a = numerator / denominator;
            double b = (sumY - a * sumX) / n;

            
            for (int i = 0; i < n; ++i) {
                double predictedY = b + a * x[i];
                double error = predictedY - y[i];
                loss_value += error * error;
            }
            mse = loss_value/(1.0*n);
            }
    }

    structure_results = structure_results + ";losses;" + std::to_string(loss_value)+ ";" + std::to_string(mse);

     structure_results = structure_results + ";" + data_node_info + ";" + model_node_info;
    saveToCSV(structure_results,structure_output);
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

    std::vector<KEY_TYPE> lookup_keys = get_search_keys(legitimate_data, values_size, 1000000, seed);
    // save_data_bin(lookup_keys,lookups_output);
    // lookup_keys.clear();
    // lookup_keys = read_data_bin(lookups_output);

    std::vector<KEY_TYPE> lookup_keys_zipf = get_search_keys_zipf(legitimate_data, values_size, 1000000, seed);


    if(poison){
        save_data_bin(altered_data,changed_data_output);
    }
    if(!poison&benchmark){
        altered_data.clear();
        altered_data = read_data_bin(changed_data_output);
    }

    if(benchmark){
         size_t index_size_after =  index.data_size() + index.model_size() ;
        //benchmark_alex_real(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);
        std::ofstream file;
        file.open(model_output+".txt", std::ios_base::app);
    
        file << "ALEX" << ";" << data_name << ";"  << insert << ";" << insert_threshold << ";"<< legitimate_data.size() << ";" << orignal_P << ";"
        << num_total_poisoning << ";"<< max_model_height << ";" << index_size_after << ";" <<
        total_node_count << ";" << total_data_count << ";" << std::endl;

        file.close();

        //std::cout << "Performance Before " <<std::endl;
        //index_original.bulk_load(values, legitimate_data.size());
       
        // benchmark_alex_real_seperate(index_original, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, original_output);
        // benchmark_alex_real_seperate(index_original, changed_data, changed_data, "ALEX", data_name, poi_thres, original_output);
        // benchmark_alex_real_seperate(index_original, lookup_keys, lookup_keys, "ALEX", data_name+" lookup", poi_thres, original_output);
        // benchmark_alex_real_seperate(index_original, lookup_keys_zipf, lookup_keys_zipf, "ALEX", data_name+" lookup_zipf", poi_thres, original_output);
       
        // benchmark_alex_real(index_original, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, original_output);
        // benchmark_alex_real(index_original, changed_data, changed_data, "ALEX", data_name, poi_thres, original_changed_output);
        
        // //Something is happening when I enter these data to the nfl the data changes
        // //ISSUE is that after the transformation all the data becomes 0 
        // //Some issue with the data type

        std::cout << "Performance " << std::endl;
        // benchmark_alex_real_seperate(index, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, original_output);
        // benchmark_alex_real_seperate(index, changed_data, changed_data, "ALEX", data_name, poi_thres, original_changed_output);

       //benchmark_alex_real_seperate(index, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, poisoned_output);
        // benchmark_alex_real_seperate(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);
        // benchmark_alex_real_seperate(index, lookup_keys, lookup_keys, "ALEX", data_name+" lookup", poi_thres, poisoned_changed_output);
        // benchmark_alex_real_seperate(index, lookup_keys_zipf, lookup_keys_zipf, "ALEX", data_name+" lookup_zipf", poi_thres, poisoned_changed_output);

         benchmark_alex_real_seperate(index,altered_data, altered_data, "ALEX", data_name+ "_changed", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);
         benchmark_alex_real_seperate(index,lookup_keys, lookup_keys, "ALEX", data_name+"_lookup", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);
         benchmark_alex_real_seperate(index,lookup_keys_zipf, lookup_keys_zipf, "ALEX", data_name+ "_lookup_zipf", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);


        std::string model_results = "ALEX;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          (poi_size_str)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" +
        std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        // std::cout << test_output <<std::endl; 
        saveToCSV(model_results,output_model_csv);
         
        // benchmark_alex_real(index, legitimate_data, legitimate_data, "ALEX", data_name, poi_thres, poisoned_output);
        // benchmark_alex_real(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);
        
       
        //std::cout << "Performance Other " <<std::endl;

       

        //std::ofstream file;
        file.open(model_output+".txt", std::ios_base::app);
        //file << regression_name << ";" << data_name << ";" << poisoning_threshold << ";" << data.size() << ";" << lookups.size() << ";" << mean << ";" << median << ";" << log_error << ";" << d_log_error << ";" << mse_error << ";" << build_time << std::endl;
        file << "ALEX " << ";" << data_name << ";" << data_size << ";"<< method <<";" << original_poi_thres << ";" << num_total_poisoning <<";"<< original_nodes_number << ";" << original_model_nodes_number << ";" << original_data_nodes_number << ";" 
        << nodes_fin.size() << ";" << model_nodes_fin.size() << ";" << data_nodes_fin.size() <<  ";"  << replaced_models_number << ";" << tra_time << std::endl;

        file.close();
        
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

        

        // I need to create the index here and then evaluate
        // //For PGM
        // //===========
    //pgm::PGMIndex<KEY_TYPE, 4,4> index_pgm(legitimate_data);
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
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;
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

//The functions required
//==============================================
//==============================================
//==============================================
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
    std::copy_if(vec2.begin(), vec2.end(), std::back_inserter(*existingVector), [&vec1](KEY_TYPE element) {
        return std::find(vec1.begin(), vec1.end(), element) == vec1.end();
    });
}

//Get the data node data
void get_data_node_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,data_node_type* child, std::vector<KEY_TYPE> *model_node_data){
       
 
    data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
    data_it(child,0);


//Check if there are data

    while(!data_it.is_end()){
    model_node_data->push_back(data_it.key());
    
    data_it.operator++(0);
    // if(node_it.is_end()){
    //     break;
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


            model_nodes_useful[i]->cost_diff_ = (model_nodes_useful[i]->cost_ - model_nodes_useful[i]->children_cost);
            
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
    unsigned seed = 123;

    std::random_device rd;
    std::mt19937 generator(123);
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

        levels += static_cast<short>(leaf->level_);
       
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

       // *const_search = *const_search*3;


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
    // std::mt19937_64 gen_payload(std::random_device{}());
    // auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[size];
    // //std::mt19937_64 gen_payload(std::random_device{}());

    //     for (int i = 0; i < size; i++) {
    //         values[i].first = data[i];
    //         values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    //     }

    // return values;

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
        // std:: cout << "model size " << model->num_keys_model_ << std::endl;
        // std:: cout << "cur size " << cur_size << std::endl;
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

void shiftVector(std::vector<KEY_TYPE>& vec) {
    if (vec.empty()) {
        return; // If the vector is empty, no need to shift
    }

    // Find the minimum value in the vector
    KEY_TYPE minVal = vec[0];

    // Shift all elements by subtracting the minimum value
    KEY_TYPE offset = -minVal;

    // Add the offset to each element in the vector
    for (KEY_TYPE& num : vec) {
        num += offset;
    }
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