
// USE THIS TO RUN THE CODE
//  g++ -std=c++17 -arch arm64 -I /opt/homebrew/Cellar/boost/1.82.0_1/include -I /opt/intel/oneapi/mkl/2023.2.0/include -L /opt/homebrew/Cellar/boost/1.82.0_1/lib -L /opt/intel/oneapi/mkl/2023.2.0/lib -lboost_system -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -o nfl_test nfl_test.cpp
#include <iostream>
//#include <boost/optional.hpp>
//#include "afli/afli.h"
///g++ nfl_test.cpp -o nfl_test -lmkl
///g++ nfl_test.cpp -o nfl_test -L/opt/intel/oneapi/mkl/2023.2.0/lib/ -lmkl
//source /opt/intel/oneapi/mkl/2023.2.0/bin/mklvars.sh intel64
//$ source /opt/intel/oneapi/setvars.sh --force intel64
 ///opt/intel/oneapi/mkl
 // g++ -o nfl_test nfl_test.cpp -I/opt/intel/oneapi/mkl/2023.2.0/include -L/opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_intel_ilp64.a -lmkl -Wl,-rpath,/opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_intel_ilp64.a
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include nfl_test.cpp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

// USE THIS TO RUN THE CODE
//source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include nfl_test.cpp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
#include "nfl/nfl.h"
#include "benchmark/workload.h"
#include "util/common.h"

#include <vector>
#include <algorithm>
#include <utility>

void save_to_csv(const std::vector<int64_t>& data, const std::vector<double>& time_total,const std::vector<std::string>& headers, const std::string& filename);
void save_to_csv_insert(const std::vector<int64_t>& data,const std::vector<std::string>& headers, const std::string& filename);

// Get the structure from bechmark.cc
struct NFLConfig {
  int bucket_size;
  int aggregate_size;
  std::string weights_path;

  NFLConfig(std::string path) {
    bucket_size = -1;
    aggregate_size = 0;
    weights_path = "";
    if (path != "") {
      std::ifstream in(path, std::ios::in);
      if (in.is_open()) {
        while (!in.eof()) {
          std::string kv;
          in >> kv;
          std::string::size_type n = kv.find("=");
          if (n != std::string::npos) {
            std::string key = kv.substr(0, n);
            std::string val = kv.substr(n + 1);
            if (key == "bucket_size") {
              bucket_size = std::stoi(val);
            } else if (key == "aggregate_size") {
              aggregate_size = std::stoi(val);
            } else if (key == "weights_path") {
              weights_path = val;
            }
          }
        }
        in.close();
      }
    }
  }
};

enum OperationType {
  kBulkLoad = 0,
  kQuery = 1,
  kInsert = 2,
  kUpdate = 3,
  kDelete = 4
};

template<typename KT, typename VT>
struct Request {
  OperationType op;
  std::pair<KT, VT> kv;
};
// Get the data from workload.h



int main() {
    std::vector<std::string> headers = {"Key","Total Time"};
    //std::string config_path = "../configs/nfl_wiki-ts-200M-0R-zipf.in";
    
    std::string config_path = "../configs/nfl_books-200M-0R-zipf.in";
    //std::string config_path = "../configs/nfl_fb-200M-0R-zipf.in";

    //std::string path = "wiki_ts_200M_uint64";
    std::string path = "books_200M_uint32";
    //std::string path = "fb_200M_uint64";
    
    //KT and VT are int64-t
    //KVT is std::pair<int64_t, int64_t>
    std::cout << "Creating the data" << std::endl;
    
    //IMPORTANT
     // READ FROM A FILE HERE
     //----------------------
     //Check what type of data they have 

     //make it into key value pairs

     //
    // std::vector<std::pair<int64_t, int64_t>> init_data = {
    //     {1, 10},
    //     {1, 10},
    //     {3, 30},
    //     {4, 40},
    //     {5, 50},
    //     {6, 60},
    //     {7, 70},
    //     {8, 80},
    //     {9, 90},
    //     {10, 100}
    // };

    std::vector<std::pair<int64_t, int64_t>> init_data;

    // Add as the requests the init_data. Just need to set the operation types.

    std::vector<Request<int64_t, int64_t>> requests;
    Request<int64_t, int64_t> request;

 // THE ISSUE IS WITH REQUESTS, WHEN THE VALUE IS NOT MET IT GIVES A SEGMENTATION FAULT

    std::cout << "File SEARCHING" << std::endl;
    std::ifstream in(path, std::ios::binary | std::ios::in);
    if (!in.is_open()) {
      std::cout << "File [" << path << "] does not exist HERE" << std::endl;
      exit(-1);
    }
    int num_reqs = 0;
    in.read((char*)&num_reqs, sizeof(int));
    // init_data.reserve(num_reqs);
    // requests.reserve(num_reqs);
    
    for (int i = 0; i < num_reqs; ++ i) {
      Request<int64_t, int64_t> req;
      in.read((char*)&req, sizeof(Request<int64_t, int64_t>));
      if (req.op == kBulkLoad) {
        //init_data.push_back(req.kv);
      } else {
        //std::cout << req.op<< std::endl;
        //init_data.push_back(req.kv);
        //requests.push_back(req);
        if(init_data.size()<5000000){
          init_data.push_back(req.kv);
          requests.push_back(req);
        }
        
        else{

        }
        if(req.op == kQuery){
          //requests.push_back(req);
          //init_data.push_back(req.kv);
        }
        // if(req.op == kQuery && req.op == kInsert){
        //   init_data.push_back(req.kv);
        // }
        
          //requests.push_back(req);
      }
      //When I increase this above a certain number (3300000) it gives a segmentation error. I think it has something to do with tail conflicts
      // if(init_data.size()>3700000)
      //   break;
    }
    in.close();
    
    std::cout << "Init Data Size: " << init_data.size() << " and " << init_data[0].second << std::endl;
    std::cout << "Requests  Size: " << requests.size() << std::endl;

    // //Getting 50% of the dataset
    // std::vector<std::pair<int64_t, int64_t>> train_data;
    // std::vector<std::pair<int64_t, int64_t>> insert_data;
    // //std::random_device rd;
    // std::mt19937 rng(1234);

    // // Shuffle the input vector
    // std::vector<std::pair<int64_t, int64_t>> shuffledInput = init_data;
    // std::shuffle(shuffledInput.begin(), shuffledInput.end(), rng);

    // // Calculate the number of pairs in each half
    // size_t halfSize = shuffledInput.size() / 2;

    // // Split the shuffled input into two halves
    // train_data.assign(shuffledInput.begin(), shuffledInput.begin() + halfSize);
    // insert_data.assign(shuffledInput.begin() + halfSize, shuffledInput.end());

    std::sort(init_data.begin(), init_data.end());
    std::vector<int64_t> keys(init_data.size());

    std::transform(init_data.begin(), init_data.end(), keys.begin(),
                   [](const std::pair<int64_t, int64_t>& pair) { return pair.first; });

    //UNCOMMENT IF I NEEF THE DATASET
    //save_to_csv_insert(keys,{"Key"},"data.csv");

    // std::cout << "Pair values size: " << train_data.size() << " first " << train_data[0].first << " and " << train_data[0].second << std::endl;
    // std::cout << "Pair values size: " << insert_data.size() << " first " << insert_data[1].first << " and " << insert_data[1].second << std::endl;
    std::cout << "Requests values size: " << requests.size() << std::endl;
    std::cout << "Request values 5429366: first " << requests[5429366].kv.first << " and " << requests[5429366].kv.second <<"OperationType" << requests[5429366].op <<std::endl;
    std::cout << "Request values 5429367: first " << requests[5429367].kv.first << " and " << requests[5429367].kv.second <<"OperationType" << requests[5429367].op <<std::endl;
    std::cout << "Request values 5429368: first " << requests[5429368].kv.first << " and " << requests[5429368].kv.second <<"OperationType" << requests[5429368].op <<std::endl;

    OperationType targetOperation = kInsert;
    int insertCount = std::count_if(requests.begin(), requests.end(),
                                    [targetOperation](const auto& request) {
                                        return request.op == targetOperation;
                                    });

    std::cout << "Number of kInsert requests: " << insertCount << std::endl;

    targetOperation = kQuery;
    int queryCount = std::count_if(requests.begin(), requests.end(),
                                    [targetOperation](const auto& request) {
                                        return request.op == targetOperation;
                                    });

    std::cout << "Number of kQuery requests: " << queryCount << std::endl;

    targetOperation = kBulkLoad;
    int bulkCount = std::count_if(requests.begin(), requests.end(),
                                    [targetOperation](const auto& request) {
                                        return request.op == targetOperation;
                                    });

    std::cout << "Number of kBulkLoad requests: " << bulkCount << std::endl;

    targetOperation = kUpdate;
    int updateCount = std::count_if(requests.begin(), requests.end(),
                                    [targetOperation](const auto& request) {
                                        return request.op == targetOperation;
                                    });

    std::cout << "Number of kUpdates requests: " << updateCount << std::endl;

    targetOperation = kDelete;
    int deleteCount = std::count_if(requests.begin(), requests.end(),
                                    [targetOperation](const auto& request) {
                                        return request.op == targetOperation;
                                    });

    std::cout << "Number of kDeletes requests: " << deleteCount << std::endl;



    //CREATE THE REQUESTS HERE
    //----------------------  
    // std::vector<Request<int64_t, int64_t>> requests;
    // Request<int64_t, int64_t> request;
    // std::cout << "Creating the requests" << std::endl;

    // for(int i =0; i<init_data.size();i++){
    //     std::cout << "Pair values: " << init_data[i].first << " and " << init_data[i].second << std::endl;
        
    //     request.kv = init_data[i];
    //     std::cout << "setting the request type" << std::endl;
    //     request.op = kQuery; 
    //     requests.push_back(request);
    // }
  std::cout << "Created the requests" << std::endl;

  const double conflicts_decay = 0.1;
    // Check the order of load data.
    for (int i = 1; i < init_data.size(); ++ i) {
      if (init_data[i].first < init_data[i - 1].first || 
          nfl::compare(init_data[i].first, init_data[i - 1].first)) {
        std::cout << std::fixed 
                  << std::setprecision(std::numeric_limits<int>::digits10) 
                  << "Duplicated data" << std::endl << i - 1 << "th key [" 
                  << init_data[i - 1].first << "]" << std::endl << i 
                  << "th key [" << init_data[i].first << "]" << std::endl;
      }
      nfl::assert_p(init_data[i].first > init_data[i - 1].first, 
              "Unordered case in load data");
    }
    std::cout << "Bulk Loading" << std::endl;

    //---------------------------------------------
    //MODIFY THE GENERATE CONFIGS.SH TO CREATE THE NFL CONFIGS AND THEN USE THEM.
    NFLConfig config(config_path);
    
    // Start to bulk load
    int batch_size = 100;
    nfl::NFL<int64_t, int64_t> nfl(config.weights_path, batch_size);
    std::cout << "FILE CONFIG" << std::endl;
    uint32_t tail_conflicts = nfl.auto_switch(init_data.data(), 
                                              init_data.size());

    std::cout << "Tail conflicts done" << std::endl;

                                              
    nfl.bulk_load(init_data.data(), init_data.size(), tail_conflicts);

    bool show_stat = false;
    bool show_incremental_throughputs =false;
    std::cout << "Created the index" << std::endl;
    int inst_count =0;
    int query_count =0;
    
   
    
    int reps =2;
    std::vector<double> times_total;
    std::vector<double> times_insert;
    std::vector<int64_t> data;
    std::vector<int64_t> data_insert;
    std::cout << data_insert.size()<<std::endl;
    bool is_rep =true;

    bool is_query =true;
    bool is_insert =false;

    int q_num=1;
    int inst_count_limit = 1000000;
    //int inst_count_limit = 10;

    if(is_insert){
      q_num =2;
    }
    std::string filename = "TEST 4 NFL results for books read only";
    //int batch_size = 5429367;
    //ISSUE WITH 5429367
    nfl::ExperimentalResults exp_res(batch_size);
    std::vector<std::pair<int64_t, int64_t>> batch_data;
    batch_data.reserve(batch_size);
    // Perform requests in batch
    //int num_batches = std::ceil(requests.size() * 1. / batch_size);

    //WHEN NOT IN BATCHES
    int num_batches = std::ceil(requests.size());
    //Change according to request type
    //int num_batches = std::ceil(queryCount * 1. / batch_size);
    std::cout << "Number of batches" <<num_batches<< std::endl;
    exp_res.latencies.reserve(num_batches * 3);
    exp_res.need_compute.reserve(num_batches * 3);
    //std::cout << "Number of exp_lat" <<exp_res.latencies.size()<< std::endl;
    int time_transform =0;
    for (int batch_idx = 0; batch_idx < num_batches; ++ batch_idx) {
      is_query =true;
      is_insert =false;

      batch_data.clear();
      int l = batch_idx * batch_size;
      int r = std::min((batch_idx + 1) * batch_size, 
                        static_cast<int>(requests.size()));
        //ADD TIME HERE TOO
        auto start_transform = std::chrono::high_resolution_clock::now();
        nfl.transform(batch_data.data(), batch_data.size());
        auto end_transform = std::chrono::high_resolution_clock::now();
        time_transform += std::chrono::duration_cast<std::chrono::nanoseconds>(end_transform - start_transform).count();

      //std::cout << "R is " <<r<< std::endl;
      for (int i = l; i < r; ++ i) {
        //std::cout << "doing for " <<i<< std::endl;
        batch_data.push_back(requests[i].kv);
      }
      //std::cout << "Created batch " << std::endl;
    
    
      int64_t val_sum = 0;
      int time_t =0;
      int time_i =0;
      // Perform requests
      // ADD THE TIMING CODE HERE IF NEEDED
      auto start = std::chrono::high_resolution_clock::now();
      for(int q =0;q<q_num;q++){
          //std::cout << "Q:" << q<<std::endl;
          if(q > 0){
            //std::cout << data_insert.size()<<std::endl;
            is_query = true;
            is_insert = false;
          }
          //Add the code here
        
        
        
        for (int i = l; i < r; ++ i) {
          time_t = 0;
          time_i = 0;
          int data_idx = i - l;
          if (is_query) {
            data.push_back(batch_data[data_idx].first);
            //Running for the replicates if query
            if(is_rep){
              for(int rep=0; rep<reps; rep++){
              auto start_q = std::chrono::high_resolution_clock::now();
              auto it = nfl.find(data_idx);
              auto end_q = std::chrono::high_resolution_clock::now();
              time_t += std::chrono::duration_cast<std::chrono::nanoseconds>(end_q - start_q).count();
              if (!it.is_end()) {
                val_sum += it.value();
              }
            }
            }
            else{
              auto it = nfl.find(data_idx);
              if (!it.is_end()) {
                val_sum += it.value();
              }
            }
            // auto it = afli.find(batch_data[data_idx].first);
            // if (!it.is_end()) {
            //   val_sum += it.value();
            // }
          } else if (requests[i].op == kUpdate) {
            //bool res = afli.update(batch_data[data_idx]);
          } else if (requests[i].op == kInsert) {
            // if(i>5429366){
            //   continue;
            // }
            if(is_insert && inst_count<inst_count_limit){
              //std::cout << "#############"<<std::endl;
              inst_count++;
              data_insert.push_back(batch_data[data_idx].first);
              auto start_i = std::chrono::high_resolution_clock::now();
              nfl.insert(data_idx);
              auto end_i = std::chrono::high_resolution_clock::now();
              time_i = std::chrono::duration_cast<std::chrono::nanoseconds>(end_i - start_i).count();
              
            }
            
          } else if (requests[i].op == kDelete) {
            //int res = afli.remove(batch_data[data_idx].first);
          }
          //std::cout << "did for " <<i<< std::endl;
          
          if(is_query){
            times_total.push_back(time_t/(double)(reps));
          }
          if(is_insert){
            times_insert.push_back(time_i);
          }
          
        }
      }
      auto end = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      exp_res.sum_indexing_time += time;
      exp_res.num_requests += batch_data.size();
      exp_res.latencies.push_back({0, time});
      exp_res.step();
      //std::cout << "Created batch data" << batch_idx<<std::endl;
    }
    // is_insert =true;
    // is_query = false;
    save_to_csv(data,times_total,headers,filename+".csv");
    if(is_insert){
      save_to_csv_insert(data_insert,{"Key"},filename+"inserts.csv");
    }
    // if(is_query){
    //   save_to_csv(data,times_total,headers,filename+".csv");
    // }
    // else{
    //   save_to_csv(data,times_total,times_insert,headers,filename+".csv");
    // }
    
    std::cout<<"Showing "<< filename<<std::endl;
    exp_res.model_size = nfl.model_size();
    exp_res.index_size = nfl.index_size();

    if (show_stat) {
      nfl.print_stats();
    }

    if (show_incremental_throughputs) {
      exp_res.show_incremental_throughputs();
    } else {
      //Segmentation issue here
      exp_res.show(true);
    }
     // nfl::LinearModel<int64_t>*model2 = new nfl::LinearModel<int64_t>();
    // int64_t pred = model2->predict(2);   

    // std::cout<<pred<<std::endl;

    // uint32_t first_pos = std::min((std::max(static_cast<int64_t>(0), static_cast<int64_t>(0))), 
    //                               static_cast<int64_t>(10 - 1));

    std::cout<<"Number of inserts" << data_insert.size()<<std::endl;
    std::cout<<"Number of queries" << data.size()<<std::endl;
    std::cout<<"Time Transform" << time_transform <<std::endl;
    std::cout<<"Time Transform" << time_transform/(1.0*5000000)<<std::endl;

    
    return 0;
}
void save_to_csv(const std::vector<int64_t>& data, const std::vector<double>& time_total,const std::vector<std::string>& headers, const std::string& filename) {
    
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
    for (size_t i = 0; i < data.size(); i++) {
        if (i < data.size()) {
            outputFile << data[i];
        }
        outputFile << ",";

        if (i < time_total.size()) {
            outputFile << time_total[i];
        }
        outputFile << "\n";
    }

    outputFile.close();
}

void save_to_csv_insert(const std::vector<int64_t>& data,const std::vector<std::string>& headers, const std::string& filename) {
    
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
    for (size_t i = 0; i < data.size(); i++) {
        if (i < data.size()) {
            outputFile << data[i];
        }
        outputFile << "\n";
    }

    outputFile.close();
}

