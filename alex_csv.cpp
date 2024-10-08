// RUN using
// g++ alex_csv.cpp -std=c++17 -o alex_csv -march=native -mpopcnt
// ./alex_csv test 1000000 1 0.2 0

#define KEY_TYPE uint64_t
#define PAYLOAD_TYPE double


#include <iostream>
#include <fstream>
#include "src/log_regression.h"
#include "src/competitor_regression.h"
#include "src/irls.h"


#include "src/helpers/io_handler_real.h"
#include "src/smooth_simple.h"
#include "src/helpers/alex_benchmark_real_new_ind.h"
#include "src/fast_brute_force_real.h"

#include "src/theil_sen.h"

#include <unordered_set>

typedef alex::AlexNode<KEY_TYPE, PAYLOAD_TYPE> node_type;

typedef alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE> model_node_type;
  //This is the data node
typedef alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> data_node_type;

bool check_duplicates(std::vector<KEY_TYPE> data);

std::pair<uint64_t, PAYLOAD_TYPE> * create_values(std::vector<KEY_TYPE> data,int * size);

void get_constants(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<KEY_TYPE> data, int max_data, double *const_search, double *const_traversal);

void get_all_nodes(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index, std::vector<node_type*> *nodes,std::vector<model_node_type*> *model_nodes,
std::vector<data_node_type*> *data_nodes,int *max_model_height);
void get_useful_model_nodes( std::vector<model_node_type*> *model_nodes, std::vector<model_node_type*> *model_nodes_useful, 
std::vector<model_node_type*> *model_nodes_by_level);

std::vector<double>  calculate_model_node_new_cost(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);

void calculate_data_node_new_cost_new(std::vector<model_node_type*>& model_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);
void get_children_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,model_node_type* model, std::vector<KEY_TYPE> *model_node_data, 
std::vector<PAYLOAD_TYPE> *model_node_payload, double *children_cost,
double const_search, double const_traversal, bool calculate_childern_costs);

std::vector<std::pair<model_node_type*,double>> calculate_best_node_cost_prev(std::vector<model_node_type*>& model_nodes_useful);

void calculate_children_new_cost(std::vector<model_node_type*>& model_nodes_useful,alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index);
void calculate_children_new_cost(std::vector<model_node_type*>& model_nodes_useful, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);

void calculate_data_node_new_cost(std::vector<data_node_type*>& data_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);
void calculate_cost_difference(std::vector<model_node_type*>& model_nodes_useful);
void get_useful_model_nodes_by_level( std::vector<model_node_type*> *model_nodes, std::vector<model_node_type*> *model_nodes_useful, int level);

void shift_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
void shift_back_vector(std::vector<KEY_TYPE>* vec, KEY_TYPE key);
std::pair<KEY_TYPE, PAYLOAD_TYPE> * get_poisoned_values(std::vector<KEY_TYPE> data_leg, std::vector<KEY_TYPE> data_poi,
std::vector<PAYLOAD_TYPE> payload, int size);
data_node_type * create_new_data_node(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index , model_node_type* model, std::pair<uint64_t, double> * values, std::vector<KEY_TYPE> leg,int size, int cur_size,
double const_search, double const_traversal);

void get_data_node_data(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,data_node_type* child, std::vector<KEY_TYPE> *model_node_data);
void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1, std::vector<KEY_TYPE> vec2);

std::vector<model_node_type*> replace_parents(model_node_type* model,data_node_type *new_data_node, std::vector<model_node_type*> *model_nodes_by_level, int *model_idx
,int *parent_model_idx);

void update_data_structure(model_node_type* best_node,std::vector<model_node_type*> *model_nodes,
std::vector<model_node_type*> *model_nodes_by_level, std::vector<model_node_type*>* model_nodes_useful,
std::vector<model_node_type*>* model_nodes_useful_by_level, std::vector<model_node_type*> parents,
alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const,std::vector<model_node_type*>* new_useful_parents);
void get_data_nodes_by_level( std::vector<data_node_type*> *data_nodes, std::vector<data_node_type*> *data_nodes_by_level, 
int level);
void get_data_by_level(std::vector<data_node_type*> *data_nodes_by_level, std::vector<KEY_TYPE> * data);
bool check_all_data_nodes(model_node_type* node);
bool compare_models_new(const model_node_type* a, const model_node_type* b);
void get_children_cost(model_node_type* model, double *children_cost,alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const);
bool compare_searches_des(const data_node_type* a, const data_node_type* b);

//====================



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
    std::string cutoff_str = argv[6];
    
    bool benchmark = 1;
    
    std::string data_folder = "data/";
    std::string data_output = data_folder + dataset_name+".bin";

    bool poison = (argv[3][0] == '1');

    

    KEY_TYPE orignal_P = static_cast<KEY_TYPE>(std::stoi(poi_size_str));

    
    double insert_threshold = std::stod(insert_prop_str);
    bool insert = (insert_threshold > 0.0);
    int cutoff_ori = std::stoi(cutoff_str);

    std::vector<KEY_TYPE> legitimate_data = read_data_bin(data_output);
    

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
    
    double poi_thres = original_poi_thres;
    long pois_count = poi_thres*legitimate_data.size();
    bool use_new_cost = true;
    int num_total_poisoning = 0;
    int max_check = 20;
    
        
    srand(12345);

    //Get the dataname, size and the method to save them
    std::string data_name = argv[1];
    std::string data_name2 = "fb_200M_uint64";
    
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

    std::string performance_output = "results/alex/alex_original_performance.csv";
    std::string structure_output = "results/alex/alex_original_structure.csv";

    if(poison){
        performance_output = "results/alex/alex_smooth_performance.csv";
        structure_output = "results/alex/alex_smooth_structure.csv";
    }
     std::string changed_data_output = "results/changed_data.bin";
     std::string output_model_csv = "results/alex/alex_model_output_good.csv";

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
        
        std::cout << "Insert." << values_size << std::endl;

        std::string bulk_index_output = data_folder + "splits/"+  dataset_name +"_bulk.bin";
        bulk_load_indexes = readIndexesFromFile(bulk_index_output);

        std::vector<KEY_TYPE> subvector(bulk_load_indexes.size());
        subvector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(subvector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });
        
       std::cout << "insert_indexes "  << insert_indexes.size() << std::endl;

        std::cout << "subvector "  << subvector.size() << std::endl;
        numberOfIndexes = subvector.size();
        auto values = create_values(subvector,&numberOfIndexes);
        std::cout << "Created values" << std::endl;
       
        index.bulk_load(values, numberOfIndexes);
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
    // std::vector<int>().swap(bulk_load_indexes);

    std::cout << "============================" << std::endl;

    auto stop_build = std::chrono::high_resolution_clock::now();


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
    
    // Add a different condition OR a different method
    //======================
    int num_pois = model_nodes.size() -1;

    int replaced_models_number = 0;

    double pois_count_model = (pois_count/(1.0*num_pois));

    std::vector<model_node_type*> model_nodes_useful;
   
    std::vector<model_node_type*> *model_nodes_by_level = new std::vector<model_node_type*>[max_model_height + 1];

    //Get the model nodes by level and the useful (all children are data nodes) models
    //I also calculate the children sum for these . (CHILDREN SUM IS  NO LONGER NEEDED)
    get_useful_model_nodes(&model_nodes, &model_nodes_useful, model_nodes_by_level);
    //Get the model nodes with all the children being data nodes.
    std::cout << "DONE: models useful and models by level" << std::endl;
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
        model_node_type* best_node;

        //Sort the models in this level by their cost difference
        std::sort(model_nodes_useful_by_level.begin(), model_nodes_useful_by_level.end(), compare_models_new);

        
        //To store if the model has converted for the current index
        bool model_converted = false;
        inner_idx = 0;

        // do until the model is not converted or more than 10 percent of the models are checked
        //METHOD 1-5
        int size_to_check = model_nodes_useful_by_level.size();
        
        tried = 0;
        sucess =  0; 

        //inner while to loop through the inner index until we reach the max models useful in the level
        //Stop if we have reached the max_check limit without changind a model to a data node.
        

        while(inner_idx < size_to_check && tried < max_check){

            if(inner_idx%1000 ==0){
                std::cout << "starting " << inner_idx<<std::endl;
            }

            best_node = model_nodes_useful_by_level[inner_idx];


            //Now get the data for this model node
            double best_children_cost = 0;
            std::vector<KEY_TYPE> model_node_data;
            std::vector<PAYLOAD_TYPE> model_node_payload;

            //Get the children data (aka all the model data) for the best node
            //We don't need to recalculate the children cost because we have already calculated it before (useful models method last para false here)
            //Also get the payloads
            get_children_data(index, best_node, &model_node_data,&model_node_payload, &best_children_cost, const_search, const_traversal, false);

            //Poison using the model data or till the cost is less than the children nodes sum of costs.
            //Vector to hold the poisoned model data or if no poisoning use this as well. I put these outside the if becasue I need to use it in else as well
            std::vector<KEY_TYPE> poisoned_model_data;
            //poi_thres = (pois_count_model/(1.0*model_node_data.size()));
            poi_thres = original_poi_thres;
            //poi_thres = original_poi_thres*(max_model_height-current_level+1);
            
            int poisoned_model_data_size;
            data_node_type* new_data_node;
            int current_poi =0;

            sucess++;
            int cutoff = cutoff_ori ;
            //If the cost differebce is greater than 0 meaning model cost > children cost try poisoning
            if(best_node->cost_diff_ > (cutoff) || true){
                
                KEY_TYPE first_val = model_node_data[0];
                shift_vector(&model_node_data, first_val);

                //Do the poisoning for the model data
                poisoned_model_data = perform_poisoning(model_node_data, poi_thres);

                shift_back_vector(&poisoned_model_data, first_val);
                shift_back_vector(&model_node_data, first_val);


                poisoned_model_data_size = poisoned_model_data.size();

                //get new values for the poisoning data while using the old ones for the existing
                auto poisoned_values = get_poisoned_values(model_node_data, poisoned_model_data, model_node_payload, poisoned_model_data_size);
               
                //Create a new data node using these poisoned values

               //ALso check the max data node number here itself
               //Add a constant that if there were no new points also save
                 if(index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size ){
                    inner_idx++;
                    tried++;
                    continue;
                }

                new_data_node = create_new_data_node(index , best_node, poisoned_values,model_node_data, poisoned_model_data_size, model_node_data.size(),
                const_search, const_traversal);
                current_poi = poisoned_model_data_size - model_node_data.size();
            }
            //Can directly create a new data node instead of a model
            else{
                

                    auto values_real = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[model_node_data.size()];
                    //std::mt19937_64 gen_payload(std::random_device{}());

                        for (int i = 0; i < model_node_data.size(); i++) {
                            values_real[i].first = model_node_data[i];
                            values_real[i].second = model_node_payload[i];
                            }
                

                //Create a new data node using these values
                //Also check if index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ before creating the new data nodes
                if(index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= model_node_data.size()){
                    inner_idx++;
                    tried++;
                    continue;
                }
                
                new_data_node = create_new_data_node(index , best_node, values_real,model_node_data, model_node_data.size(), model_node_data.size(),
                const_search, const_traversal);
                

                poisoned_model_data_size = model_node_data.size();
            }


            // condition to change the nodes
            int best_level = best_node->level_;
            int best_node_idx = 0;
            bool parent_found = false;
            int parent_model_idx = 0;

            //If the cost of the new data is higher than the children or there is no more space then try the next best node
            //Change the condition here
            
            if((new_data_node->cost_ - best_node->children_cost) > (cutoff) || index.derived_params_.max_data_node_slots *data_node_type::kInitDensity_ <= poisoned_model_data_size){
                inner_idx++;
                tried++;
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
                    for (KEY_TYPE element : model_node_data) {
                    // Check if the element is not in the original vector
                        // If not found, add it to the original vector
                        full_data_changed.insert(element);
                        if (poisoned_data.find(element) == poisoned_data.end()) {
                            // Add the element to the result vector
                            altered_data.push_back(element);
                        }
                
                    }
                
            
            
            tried = 0;

            //this means we are going to replace the model node with the new data node 
                
            //find the parent of this new data node aka the old best node
            std::vector<model_node_type*> parents;
            
            new_data_node = static_cast<data_node_type*>(new_data_node);
            
            //Change the parent's child to this new data node pointer
            parents = replace_parents(best_node, new_data_node, model_nodes_by_level, &best_node_idx, &parent_model_idx);
            
            //IF I REALLY NEEDED TO CHECK, THEN I COULD GET THE NODES FROM THE INDEX AGAIN AND CHECK (ALL PREV NODES WILL COME, BUT THERE WILL BE NO LINK)]
            // TO DO
            // ADD A NEW COST MODEL IN alex_nodes.h NEED TO ACCOUNT FOR THE DEPTH AS WELL DEPTH VS.
            // ALSO NEED TO MAKE SURE ALL THE DATA CAN BE FITTED IN THE NEW DATA NODE (SIMPLE IF STATEMENT WITH MAX_NODE SIZE AFTER FINDING BEST NODE)
            // CALCULATE IT INSTEAD OF THE COST FOR BOTH MODEL AND CHILDREN
            // CHECK IF ITS BETTER.
            //  IF SO DO THE REST OF IT 

            //Removing from the vectors 

            //find and delete from model_nodes, model_nodes_by_level and model_nodes_useful also check if parent is
            // a useful model node (all children are data nodes) then add that to the model_nodes_useful
            (new_useful_parents).clear();
            update_data_structure(best_node, &model_nodes, model_nodes_by_level, &model_nodes_useful,&model_nodes_useful_by_level, parents,
            index,const_search,const_traversal,&new_useful_parents);

            replaced_models_number++;

            model_converted = true;
            
            num_total_poisoning += current_poi;
            //If model is converted then set tried back to 0 and increase the inner_idx
            if(model_converted){
                inner_idx++;
                tried = 0;
            }

        }
        
        //After a level go to the upper level
        std::cout << "=============="  <<std::endl;
        std::cout << "Done: Iteration " << current_level <<std::endl;

        current_level--;
        if(current_level  > 0){
            //Get the new models by level from model nodes useful (which we update in the update_data_structures
            get_useful_model_nodes_by_level(&model_nodes_useful, &model_nodes_useful_by_level, current_level);

        }
            
    }

    if(insert){
   
        for(int i : insert_indexes){
            index.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
        }
    }
    
    std::sort(changed_data.begin(), changed_data.end());
    std::vector<node_type*> nodes_fin;
    std::vector<model_node_type*> model_nodes_fin;
    std::vector<data_node_type*> data_nodes_fin;

    std::vector<model_node_type*> model_nodes_useful_fin;

    //Getting all the nodes from the index and the maximum height (default is 0 aka root node only)
    //Printing information at the end
    int max_model_height_fin = 0;
    

    get_all_nodes(index, &nodes_fin, &model_nodes_fin,&data_nodes_fin, &max_model_height_fin);

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
    std::cout << "Poisoned number " << num_total_poisoning <<std::endl;
    

    long tra_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_tra - start_tra).count();

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

    // std::vector<KEY_TYPE> lookup_keys = get_search_keys(legitimate_data, values_size, 1000000, seed);
    // std::vector<KEY_TYPE> lookup_keys_zipf = get_search_keys_zipf(legitimate_data, values_size, 1000000, seed);


    if(poison){
        save_data_bin(altered_data,changed_data_output);
    }
    if(!poison&benchmark){
        altered_data.clear();
        altered_data = read_data_bin(changed_data_output);
    }

    //BENCHMARKING FOR INSET & BULK
//========================
std::vector<KEY_TYPE> lookup_keys;
std::vector<KEY_TYPE> lookup_keys_zipf;

    if(insert){
        std::vector<KEY_TYPE> bulk_vector(bulk_load_indexes.size());
        bulk_vector.clear();

        std::transform(bulk_load_indexes.begin(), bulk_load_indexes.end(), std::back_inserter(bulk_vector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });
        int bulk_size = bulk_load_indexes.size();

        lookup_keys = get_search_keys(bulk_vector, bulk_size, 1000000, seed);

        lookup_keys_zipf = get_search_keys_zipf(bulk_vector, bulk_size, 1000000, seed);
        
        
        benchmark_alex_real_seperate(index,altered_data, altered_data, "ALEX", data_name+"_changed", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);

        benchmark_alex_real_seperate(index,lookup_keys, lookup_keys, "ALEX", data_name+"_lookup", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);
         benchmark_alex_real_seperate(index,lookup_keys_zipf, lookup_keys_zipf, "ALEX", data_name+ "_lookup_zipf", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);

    }
    else{
        lookup_keys = get_search_keys(legitimate_data, values_size, 1000000, seed);
        lookup_keys_zipf = get_search_keys_zipf(legitimate_data, values_size, 1000000, seed);
    }

     size_t index_size_after =  index.data_size() + index.model_size() ;

    std::string model_results = "ALEX;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          (poi_size_str)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" +
        std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        // std::cout << test_output <<std::endl; 
        saveToCSV(model_results,output_model_csv);

    if(insert){
   
        for(int j = 0; j < 5; j++){
            insert_threshold = (j+1)*0.1;
            std::string structure_results2 = "ALEX;" + dataset_name + "_insert_"+std::to_string(j) +";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning);

            std::vector<KEY_TYPE>  inserted_data;
        
            std::string insert_index_output = data_folder+ "splits/"+ data_name +"_insert"+"_"+std::to_string(j)+".bin";
            std::vector<int> insert_indexes2 = readIndexesFromFile(insert_index_output);


            auto start = std::chrono::high_resolution_clock::now();
            for(int i : insert_indexes2){
                //index_original.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
                index.insert(legitimate_data[i],static_cast<PAYLOAD_TYPE>(gen_payload()));
                //inserted_data.push_back(legitimate_data[i]);
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count()/(1.0*insert_indexes2.size());

            std::cout << "Inserted " << insert_indexes2.size() << std::endl;
        //     benchmark_lipp_real(index,inserted_data, inserted_data, "LIPP", dataset_name+ "_changed", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_lipp_real(index,lookup_keys, lookup_keys, "LIPP", dataset_name+"_lookup", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);
        //  benchmark_lipp_real(index,lookup_keys_zipf, lookup_keys_zipf, "LIPP", dataset_name+ "_lookup_zipf", poi_thres, performance_output,
        //  poison,insert,orignal_P,insert_threshold);

    //index_original.scan_and_destory_tree(nodes[0], keys, values,false);
    size_t total_node_count2 = 0;
    size_t total_data_count2 = 0;
    std::string data_node_info2 = "Data Nodes";
    std::string model_node_info2 = "Model Nodes";

    // std::string structure_results = "LIPP;" + dataset_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
    //       std::to_string(orignal_P)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
    //     + std::to_string(num_total_poisoning);

        // std::cout << test_output <<std::endl; 
    
    std::vector<node_type*> nodes_fin2;
    std::vector<model_node_type*> model_nodes_fin2;
    std::vector<data_node_type*> data_nodes_fin2;
    std::vector<model_node_type*> model_nodes_useful_fin2;

    //Getting all the nodes from the index and the maximum height (default is 0 aka root node only)
    //Printing information at the end
    int max_model_height_fin2 = 0;
    get_all_nodes(index, &nodes_fin2, &model_nodes_fin2,&data_nodes_fin2, &max_model_height_fin2);
    std::cout << "Height " << max_model_height_fin2 << std::endl;

    for(int i = 1; i <= max_model_height_fin2+1 ; i++ ){
        
       std::cout << "Performance For Level  " << i << std::endl;
        std::vector<data_node_type*> data_nodes_level;
        std::vector<model_node_type*> model_nodes_level;
        get_data_nodes_by_level(&data_nodes_fin2,&data_nodes_level,i);
        std::cout << "Data nodes number  " << data_nodes_level.size() << std::endl;

        
        // get_useful_model_nodes_by_level(&model_nodes_fin,&model_nodes_level,i);
        get_useful_model_nodes_by_level(&model_nodes_fin2,&model_nodes_level,i);
        std::cout << "model nodes number  " << model_nodes_level.size() << std::endl;
        
        data_node_info2  = data_node_info2 + ";" + std::to_string(data_nodes_level.size());
        model_node_info2  = model_node_info2 + ";" + std::to_string(model_nodes_level.size());
        
        std::vector<KEY_TYPE> level_data;
        get_data_by_level(&data_nodes_level, &level_data);

        std::cout << "data number  " << level_data.size() << std::endl;
        if(i>0){
            total_data_count2 += level_data.size();
            total_node_count2 += data_nodes_level.size() + model_nodes_level.size();
        }

        structure_results2 = structure_results2 + ";" + std::to_string(level_data.size());
    // }
        if(level_data.size() > 0){
            std::vector<KEY_TYPE> lookup_keys_level2 = get_search_keys(level_data, level_data.size(), 1000000, seed);

            // benchmark_alex_real_seperate(index,lookup_keys_level2, lookup_keys_level2, "ALEX", dataset_name+"_lookup_level_"+std::to_string(i)+ "_insert_"+std::to_string(j), 
            //     poi_thres, performance_output,
            // poison,insert,poi_size_str,insert_threshold);
        }
        
        }

        benchmark_alex_real_seperate(index,altered_data, altered_data, "ALEX", dataset_name+"_changed"+ "_insert_"+std::to_string(j),
            poi_thres, performance_output,poison,insert,poi_size_str,insert_threshold);

        benchmark_alex_real_seperate(index,lookup_keys, lookup_keys, "ALEX", dataset_name+"_lookup"+ "_insert_"+std::to_string(j),
            poi_thres, performance_output,poison,insert,poi_size_str,insert_threshold);
        benchmark_alex_real_seperate(index,lookup_keys_zipf, lookup_keys_zipf, "ALEX", dataset_name+ "_lookup_zipf"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);

         std::vector<KEY_TYPE> insert_vector(insert_indexes2.size());

        //subvector.resize(bulk_load_indexes.size());
        insert_vector.clear();

        std::transform(insert_indexes2.begin(), insert_indexes2.end(), std::back_inserter(insert_vector),
                   [&legitimate_data](int index) { return legitimate_data[index]; });

        benchmark_alex_real_seperate(index,insert_vector, insert_vector, "ALEX", dataset_name+ "_insert"+ "_insert_"+std::to_string(j), poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);
        
        structure_results2 = structure_results2 + ";" + data_node_info2 + ";" + model_node_info2+ ";"+ std::to_string(duration);
        saveToCSV(structure_results2,structure_output);

        size_t index_size_after2 =  index.data_size() + index.model_size() ;
    
    std::string model_results = "ALEX;" + data_name + "_insert_"+std::to_string(j)+ ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
          (poi_size_str)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height_fin2) + ";" + std::to_string(index_size_after2) + ";" +
        std::to_string(total_node_count2) + ";" + std::to_string(total_data_count2) + ";" + std::to_string(tra_time)  ;

        // std::cout << test_output <<std::endl; 
        saveToCSV(model_results,output_model_csv);

       // int count_issues = 0;
    //index.link_all_data_nodes();
    //index.validate_structure(true);
    //  for(int k : insert_indexes2){
    //         int ss = index.count(legitimate_data[k]);
    //         if(ss <=0){
    //             count_issues++;
                
    //         }
    //     }
        //std::cout << "Issues " << count_issues << std::endl;

        }
    }

    if(benchmark && !insert){
        //  size_t index_size_after =  index.data_size() + index.model_size() ;
        // //benchmark_alex_real(index, changed_data, changed_data, "ALEX", data_name, poi_thres, poisoned_changed_output);
        // std::ofstream file;
        // file.open(model_output+".txt", std::ios_base::app);
    
        // file << "ALEX" << ";" << data_name << ";"  << insert << ";" << insert_threshold << ";"<< legitimate_data.size() << ";" << orignal_P << ";"
        // << num_total_poisoning << ";"<< max_model_height << ";" << index_size_after << ";" <<
        // total_node_count << ";" << total_data_count << ";" << std::endl;

        // file.close();
        std::cout << "Performance " << std::endl;
       
         benchmark_alex_real_seperate(index,altered_data, altered_data, "ALEX", data_name+ "_changed", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);
         benchmark_alex_real_seperate(index,lookup_keys, lookup_keys, "ALEX", data_name+"_lookup", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);
         benchmark_alex_real_seperate(index,lookup_keys_zipf, lookup_keys_zipf, "ALEX", data_name+ "_lookup_zipf", poi_thres, performance_output,
         poison,insert,poi_size_str,insert_threshold);


        // std::string model_results = "ALEX;" + data_name + ";" + std::to_string(poison) + ";"+ std::to_string(insert) + ";" +
        //   (poi_size_str)  + ";"  + std::to_string(insert_threshold) + ";" +std::to_string(legitimate_data.size()) +";"
        // + std::to_string(num_total_poisoning) + ";" + std::to_string(max_model_height) + ";" + std::to_string(index_size_after) + ";" +
        // std::to_string(total_node_count) + ";" + std::to_string(total_data_count) + ";" + std::to_string(tra_time)  ;

        // saveToCSV(model_results,output_model_csv);
         
        // file.open(model_output+".txt", std::ios_base::app);
        // file << "ALEX " << ";" << data_name << ";" << data_size << ";"<< method <<";" << original_poi_thres << ";" << num_total_poisoning <<";"<< original_nodes_number << ";" << original_model_nodes_number << ";" << original_data_nodes_number << ";" 
        // << nodes_fin.size() << ";" << model_nodes_fin.size() << ";" << data_nodes_fin.size() <<  ";"  << replaced_models_number << ";" << tra_time << std::endl;

        // file.close();

    }

    std::cout << "Max model node height  " << max_model_height << std::endl;
    std::cout << "Time  " << tra_time/1000000000.0 << "s" << std::endl;
    std::cout << "Index Build Time  " << index_time/1000000000.0 << "s" << std::endl;
   
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


bool compare_models_new(const model_node_type* a, const model_node_type* b) {
        return a->cost_diff_ < b->cost_diff_;
}


void find_difference(std::vector<KEY_TYPE>* existingVector, std::vector<KEY_TYPE> vec1, std::vector<KEY_TYPE> vec2) {
    
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
             

}


std::vector<std::pair<model_node_type*,double>> calculate_best_node_cost_prev(std::vector<model_node_type*>& model_nodes_useful)
{
        
        double cost_diff;
        double temp_cost_diff;

        model_node_type* best_node;
        std::vector<std::pair<model_node_type*,double>> costs_diff_models;
        
        //Find the best node to use (the one with the least cost diff)
        //Create a method get best model node
        for(int i =0 ; i< model_nodes_useful.size(); i++){
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
        }
    return  costs_diff_models;
}


void calculate_children_new_cost(std::vector<model_node_type*>& model_nodes_useful, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const)
{

    for(int m =0; m < model_nodes_useful.size(); m++){
        int num_child_best = model_nodes_useful[m]->num_children_;
        node_type** best_children = model_nodes_useful[m]->children_;
        std::mt19937_64 gen_payload(std::random_device{}());
    
        for(int i = 0; i < num_child_best;i++){
            
                data_node_type* child = static_cast<data_node_type*>(best_children[i]);

                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;



            while(!data_it.is_end()){
                temp_child.push_back(data_it.key());

                data_it.operator++(0);
            }
           
            auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

        for (int i = 0; i < temp_child.size(); i++) {
            child_values[i].first = temp_child[i];
            child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        }

            child->cost_ = child->compute_new_cost(
            child_values, temp_child.size(),temp_child.size(),model_nodes_useful[m]->num_keys_model_,child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            index.params_.expected_insert_frac, &child->model_,
            index.params_.approximate_cost_computation);

            
        }
    }
}

void calculate_data_node_new_cost(std::vector<data_node_type*>& data_nodes, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
double search_const, double traversal_const)
{

    std::mt19937_64 gen_payload(std::random_device{}());
    for(int m =0; m < data_nodes.size(); m++){
        
        
                data_node_type* child = data_nodes[m];

                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;



            while(!data_it.is_end()){
               
                temp_child.push_back(data_it.key());

                data_it.operator++(0);
                
            }
            
            auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

        for (int i = 0; i < temp_child.size(); i++) {
            child_values[i].first = temp_child[i];
            child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        }

             child->cost_ = child->compute_new_cost_new(
            child_values, temp_child.size(),temp_child.size(),child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
            index.params_.expected_insert_frac, &child->model_,
            index.params_.approximate_cost_computation);

            
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
            
            int num_child = model_nodes[m]->num_children_;
            node_type** children = model_nodes[m]->children_;
           std::vector<KEY_TYPE> child_node_data;

            for(int i =0; i < num_child; i++){
                 

            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                }
                else{
                    add = false;
                }
            if(add){
                while(!data_it.is_end()){
               
                temp_child.push_back(data_it.key());
                child_node_data.push_back(data_it.key());

                data_it.operator++(0);
                
                }
                auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                for (int i = 0; i < temp_child.size(); i++) {
                    child_values[i].first = temp_child[i];
                    child_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
                }

                    child->cost_ = child->compute_new_cost(
                    child_values, temp_child.size(),temp_child.size(),model_nodes[m]->num_keys_model_ ,child->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                    index.params_.expected_insert_frac, &child->model_,
                    index.params_.approximate_cost_computation);

                    children_cost += child->cost_;

            }
            
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
            
            int num_child = model_nodes[m]->num_children_;
            node_type** children = model_nodes[m]->children_;
           std::vector<KEY_TYPE> child_node_data;

            for(int i =0; i < num_child; i++){
                 

            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                }
                else{
                    add = false;
                }

            if(add){
                while(!data_it.is_end()){
                
                temp_child.push_back(data_it.key());
                child_node_data.push_back(data_it.key());

                data_it.operator++(0);
               
                }

            }
            
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

            model_nodes[m]->cost_ = new_data_node->compute_new_cost(
                child_values, child_node_data.size(),child_node_data.size(),child_node_data.size(), model_nodes[m]->level_, data_node_type::kInitDensity_,search_const, traversal_const,
                index.params_.expected_insert_frac, data_node_model,
                index.params_.approximate_cost_computation);

                model_nodes[m]->num_keys_model_ = child_node_data.size();
            

        new_costs.push_back(model_nodes[m]->cost_);
    }
    return new_costs;
}


void calculate_cost_difference(std::vector<model_node_type*>& model_nodes_useful)
{

     std::vector<double> cost_diff;

    for(int i =0; i < model_nodes_useful.size(); i++){


            model_nodes_useful[i]->cost_diff_ = (model_nodes_useful[i]->cost_ - model_nodes_useful[i]->children_cost);    
    }
}


//Method to get the constants of search per search and traversal per level
void get_constants(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,std::vector<KEY_TYPE> data, int max_data, double *const_search, double *const_traversal){

    int const_size = std::min(static_cast<int>(round(data.size()*1)),max_data);

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
       
        
        for(int i = 0; i < num_child;i++){
            
                data_node_type* child = static_cast<data_node_type*>(children[i]);

                // Get the iterator
                data_node_type::Iterator<data_node_type, PAYLOAD_TYPE, KEY_TYPE>
                data_it(child,0);

                std::vector<KEY_TYPE> temp_child;
                std::vector<PAYLOAD_TYPE> temp_child_payload;


                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                }
                else{
                    add = false;
                }


            //Check if there are data
            if(add){
                while(!data_it.is_end()){
                
                KEY_TYPE key = data_it.key();
                PAYLOAD_TYPE payload = data_it.payload();

                model_node_data->push_back(key);
                model_node_payload->push_back(payload);
                temp_child.push_back(key);
                temp_child_payload.push_back(payload);

                //check if 0 if ok
                data_it.operator++(0);
                }
                //Add to child data

                if(calculate_childern_costs)
                {
                        child_node_data.push_back(temp_child);


                    auto child_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[temp_child.size()];

                    for (int i = 0; i < temp_child.size(); i++) {
                        child_values[i].first = temp_child[i];
                        child_values[i].second = temp_child_payload[i];
                    }
                    child->cost_ = child->compute_new_cost(
                    child_values, temp_child.size(),temp_child.size(),count_model,child->level_, data_node_type::kInitDensity_, const_search, const_traversal,
                    index.params_.expected_insert_frac, &child->model_,
                    index.params_.approximate_cost_computation);


                    *children_cost += child->cost_;
                }
            }
            
        }
        *children_cost = model->children_cost;

}

//To create the values needed for alex key value pairs
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

//create a new data node from the poisoned data
data_node_type * create_new_data_node(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index , model_node_type* model, std::pair<uint64_t, double> * values, std::vector<KEY_TYPE> leg, int size, int cur_size,
double const_search, double const_traversal){

    node_type* model_data = static_cast<node_type*>(model);

    auto new_data_node = new (data_node_type::alloc_type(index.get_allocator()).allocate(1))
            data_node_type(model->level_, index.derived_params_.max_data_node_slots,
                            index.key_comp(), index.get_allocator());

        alex::LinearModel<KEY_TYPE>* data_node_model = nullptr;

        new_data_node->bulk_load(values, size, data_node_model,
                            index.params_.approximate_model_computation);
        
        new_data_node->cost_ = new_data_node->compute_new_cost(
            values, size,cur_size,cur_size, model->level_, data_node_type::kInitDensity_,const_search, const_traversal,
            index.params_.expected_insert_frac, data_node_model,
            index.params_.approximate_cost_computation);
            
        new_data_node->is_leaf_ =true;


    return new_data_node;
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
        
    } else {
    }

    //remove from model nodes by level
    auto it2 = std::find(model_nodes_by_level[best_level].begin(), model_nodes_by_level[best_level].end(), best_node);

    if (it2 != model_nodes_by_level[best_level].end()) {
        // Pointer found, erase it from the vector
        model_nodes_by_level[best_level].erase(it2);
    } else {
    }

    //remove from model nodes useful
    auto it3 = std::find((*model_nodes_useful).begin(), (*model_nodes_useful).end(), best_node);

    if (it3 != (*model_nodes_useful).end()) {
        // Pointer found, erase it from the vector
        (*model_nodes_useful).erase(it3);
    } else {
    }

    //remove from model nodes useful by level
    //Check if the parent is now has all data nodes or not
    // If so, then add to the useful data nodes
    for(int i =0; i< parents.size();i++){
        model_node_type* parent = parents[i];

        if(check_all_data_nodes(parent)){
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
            //For each model node in the model nodes by level (all models not just useful)
            for(model_node_type* node : model_nodes_by_level[l]){
                node_type** childrens  = node->children_;

                int num_child = node->num_children_;

                for(int i =0; i< num_child; i++){
                    if(model == childrens[i]){
                        parent = node;
                        parents.push_back(node);
                        parent_found = true;
                        *model_idx = i;
                        count++;
                        
                        //Setting the new parent pointer
                        parent->children_[i] = new_data_node;
                    }
                }
                *parent_model_idx++;
            
            }
        }        
        return parents;

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
                }
                else{
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

                if (std::find((used_children).begin(), (used_children).end(), child) == (used_children).end()){
                    add = true;
                    used_children.push_back(child);
                }
                else{
                    add = false;
                }


            //Check if there are data
            if(add){
                while(!data_it.is_end()){
                temp_child.push_back(data_it.key());
                temp_child_payload.push_back(data_it.payload());

                model_data.push_back(data_it.key());
                model_payload.push_back(data_it.payload());

                //check if 0 if ok
                data_it.operator++(0);
               
                }
                //Add to child data

                
                        child_node_data.push_back(temp_child);


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
