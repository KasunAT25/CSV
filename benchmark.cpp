// RUN using
// g++ benchmark.cpp -std=c++17 -o benchmark
// ./benchmark test 1000

// ./benchmark ../../../ext/Data/fb_200M_uint64 1000
//sftp://ubuntu@172.26.131.125/home/ubuntu/ext/Data/wiki_ts_200M_uint64
// ../../../ext/data/osm_cellids_800M_uint64 1000

// I changed regression_benchmark and fast_brute force io helper
// I changed exponentional search to key_type2 from double
// I removed regression_linear and put it in pgm

#include <iostream>
#include <fstream>
#include "src/log_regression.h"
#include "src/competitor_regression.h"
#include "src/irls.h"

//#include "src/poisoning.h"
// #include "src/poisoning_new_min_real.h"
//#include "src/poisoning_new_min_decs.h"

//#include "src/poisoning_derv.h"

// #include "src/helpers/io_handler_real.h"

//#include "src/helpers/pgm_benchmark.h"
// #include "src/helpers/pgm_benchmark_real.h"
//#include "src/helpers/pgm_benchmark_alex.h"
//#include "src/helpers/pgm_benchmark_similar.h"
//#include "src/helpers/pgm_benchmark_mse.h"
// #include "src/helpers/pgm_benchmark_grad2.h"

//#include "src/helpers/regression_benchmark.h"
// #include "src/helpers/alex_benchmark_real.h"

// #include "src/helpers/io_handler.h"

//UNCOMMENT FOR NORMAL
//======================
// g++ benchmark.cpp -std=c++17 -o benchmark
// ./benchmark test 1000
// ./benchmark normal 1000

#include "src/poisoning_new_min.h"
#include "src/helpers/io_handler.h"
#include "src/helpers/regression_benchmark.h"
#include "src/helpers/pgm_benchmark.h"
#include "src/helpers/alex_benchmark.h"
#include "src/fast_brute_force.h"

//UNCOMMENT FOR REAL DATA
//======================
// g++ benchmark.cpp -std=c++17 -o benchmark
// ./benchmark ../../../ext/Data/fb_200M_uint64 1000

// #include "src/poisoning_new_min_real.h"
// #include "src/helpers/io_handler_real.h"
// #include "src/helpers/pgm_benchmark_real.h"
// #include "src/helpers/alex_benchmark_real.h"
// #include "src/fast_brute_force_real.h"

#include "src/theil_sen.h"

#include <unordered_set>

//#include "src/helpers/nfl_benchmark.h"

// USE THIS TO RUN THE CODE
// source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include benchmark.cpp -std=c++17 -o benchmark -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

// source /opt/intel/oneapi/setvars.sh --force intel64
// g++ -fopenmp -m64 -I/opt/intel/oneapi/mkl/2023.2.0/include benchmark.cpp -std=c++17 -o benchmark -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl


int main(int argc, char *argv[]){
    //srand(time(NULL));

    std::vector<KEY_TYPE> legitimate_data = parse_arguments(argc, argv);

    save_data(legitimate_data, "legitimate_data.txt");


    int num_lookups = 1000000;
    double poi_thres = 0.2;

    srand(12345);
    std::vector<double> lookups;
    lookups.resize(num_lookups);
    for (int i = 0; i < lookups.size(); i++){
        lookups[i] = legitimate_data[rand() % legitimate_data.size()];
    }
    

    std::unordered_set<int> seen;

    bool hasDuplicates = false;

    for (int num : legitimate_data) {
        if (seen.find(num) != seen.end()) {
            hasDuplicates = true;
            break;
        }
        seen.insert(num);
    }

    if (hasDuplicates) {
        std::cout << "The vector has duplicate values." << std::endl;
    } else {
        std::cout << "The vector does not have duplicate values." << std::endl;
    }


    std::cout << "Generating poisoning keys" << std::endl;

    //bool max = false;

    // Generate vector of poisoning thresholds from 0.01, 0.02, .... , 0.2
    std::vector<double> poisoning_thresholds;
    for (double threshold = poi_thres; threshold <= poi_thres; threshold = threshold + 0.2) {
        poisoning_thresholds.push_back(threshold);
    }
    std::vector<KEY_TYPE> poisoned_data;
    std::vector<KEY_TYPE> pois;
    // Perform poisoning
    auto start_all = std::chrono::high_resolution_clock::now();
    for (auto poisoning_threshold : poisoning_thresholds) {

        // std::vector<double> poisoned_data = perform_poisoning(legitimate_data, poisoning_threshold);
        //poisoned_data = perform_poisoning(legitimate_data, poisoning_threshold);

        std::string data_name = argv[1];
        std::string loutfile = "results/benchmark_fb_100000_higher.csv";

        std::string legitimate_outfile = "results/benchmark_legitimateInv2.csv";
        std::string legitimate_outfile2 = "results/benchmark_legitimate_nov_8_temp";

        std::string poisoned_outfile = "results/benchmark_poisonedInv2.csv";
        std::string poisoned_outfile2 = "results/benchmark_fb_alex_pois";

        std::string inv_poisoned_outfile = "results/benchmark_poisonedInv2.csv";

        

        std::cout << std::endl << std::endl;
        std::cout << "Benchmark results for legitimate data: " << std::endl;
       // benchmark_nfl(legitimate_data,legitimate_data,"NFL",data_name, poisoning_threshold,legitimate_outfile2);

        //benchmark_regression<simple_linear_regression_stable>(legitimate_data,lookups,"SLR",data_name, poisoning_threshold, legitimate_outfile2);
        // // benchmark_regression<create_regression_tournament_selection<LogNorm>>(legitimate_data,lookups,"LogTE",data_name,poisoning_threshold, legitimate_outfile);
        // // benchmark_regression<create_regression_tournament_selection<FastDiscreteLogNorm>>(legitimate_data,lookups,"DLogTE",data_name,poisoning_threshold, legitimate_outfile);
        // // benchmark_regression<build_regression_direct_descent>(legitimate_data,lookups,"2P",data_name, poisoning_threshold, legitimate_outfile);
        // // benchmark_regression<theil_sen>(legitimate_data,lookups,"TheilSen",data_name, poisoning_threshold, legitimate_outfile);
        // //benchmark_regression<create_regression_optimal<L1Norm>>(legitimate_data,lookups,"LAD",data_name,poisoning_threshold, legitimate_outfile);
        //benchmark_alex(legitimate_data, legitimate_data, "ALEX", data_name, poisoning_threshold, loutfile);
        // benchmark_pgm<4>(legitimate_data, lookups, "PGM_4", data_name, poisoning_threshold, legitimate_outfile2);
        // //benchmark_pgm<4>(legitimate_data, lookups, "PGM_4", data_name, poisoning_threshold, legitimate_outfile);
        // //benchmark_pgm<8>(legitimate_data, lookups, "PGM_8", data_name, poisoning_threshold, legitimate_outfile);
        // benchmark_pgm<16>(legitimate_data, lookups, "PGM_16", data_name, poisoning_threshold, legitimate_outfile2);
        // benchmark_pgm<64>(legitimate_data, lookups, "PGM_64", data_name, poisoning_threshold, legitimate_outfile);
        
       
       //pois = benchmark_pgm<4>(legitimate_data, lookups, "PGM_4", data_name, poisoning_threshold, legitimate_outfile2,true);

       //benchmark_alex(legitimate_data, legitimate_data, "ALEX", data_name, poisoning_threshold, legitimate_outfile2);
       benchmark_pgm<4>(legitimate_data, lookups, "PGM_4", data_name, poisoning_threshold, loutfile,true);
       // benchmark_pgm<4>(legitimate_data, lookups, "PGM_4", data_name, poisoning_threshold, legitimate_outfile2,false);


         


        std::cout << std::endl << std::endl;
        std::cout << "Benchmark results for poisoned data: " << std::endl;

       //benchmark_nfl(poisoned_data,legitimate_data,"NFL",data_name, poisoning_threshold,poisoned_outfile2);
        //std::cout << "Benchmark results for poisoned data: " << std::endl;


        //benchmark_regression<simple_linear_regression_stable>(poisoned_data,lookups,"SLR",data_name, poisoning_threshold,poisoned_outfile2);
        // // benchmark_regression<create_regression_tournament_selection<LogNorm>>(poisoned_data,lookups,"LogTE",data_name, poisoning_threshold,poisoned_outfile);
        // // benchmark_regression<create_regression_tournament_selection<FastDiscreteLogNorm>>(poisoned_data,lookups,"DLogTE",data_name,poisoning_threshold, poisoned_outfile);
        // // benchmark_regression<build_regression_direct_descent>(poisoned_data,lookups,"2P",data_name,poisoning_threshold, poisoned_outfile);
        // // benchmark_regression<theil_sen>(poisoned_data,lookups,"TheilSen",data_name, poisoning_threshold,poisoned_outfile);
        // //benchmark_regression<create_regression_optimal<L1Norm>>(poisoned_data,lookups,"LAD",data_name,poisoning_threshold, poisoned_outfile);
        //benchmark_alex(poisoned_data, legitimate_data, "ALEX", data_name, poisoning_threshold, poisoned_outfile2);
         //benchmark_alex(pois, legitimate_data, "ALEX", data_name, poisoning_threshold, poisoned_outfile2);
       // benchmark_pgm<4>(poisoned_data, lookups, "PGM_4", data_name, poisoning_threshold, poisoned_outfile2);
        // //benchmark_pgm<4>(poisoned_data, lookups, "PGM_4", data_name, poisoning_threshold, poisoned_outfile);
        // //benchmark_pgm<8>(poisoned_data, lookups, "PGM_8", data_name, poisoning_threshold, poisoned_outfile);
        // benchmark_pgm<16>(poisoned_data, lookups, "PGM_16", data_name, poisoning_threshold, poisoned_outfile2);
        // //benchmark_pgm<32>(poisoned_data, lookups, "PGM_32", data_name, poisoning_threshold, poisoned_outfile);
        // benchmark_pgm<64>(poisoned_data, lookups, "PGM_64", data_name, poisoning_threshold, poisoned_outfile2);


    
    }
    auto stop_all = std::chrono::high_resolution_clock::now();
    long build_time_all = std::chrono::duration_cast<std::chrono::milliseconds>(stop_all - start_all).count();

    //std::cout << "Total Time Taken for Experiment: " << build_time_all << " ms" << std::endl;


    std::string datafile = "results/data.txt"; // Change this to your desired filename

    // Open a file for writing using an ofstream object
    std::ofstream outputFile(datafile);

    if (outputFile.is_open()) {
        // Use std::copy to write the vector contents to the file
        std::copy(legitimate_data.begin(), legitimate_data.end(), std::ostream_iterator<int>(outputFile, "\n"));

        // The file will be automatically closed when the ofstream object goes out of scope

        std::cout << "Data contents saved to " << datafile << std::endl;
        outputFile.close();

    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }

    //Saving lookups
    //=========

    // std::string lookupfile = "results/looks.txt"; // Change this to your desired filename

    // // Open a file for writing using an ofstream object
    // std::ofstream outputFile2(lookupfile);

    // if (outputFile2.is_open()) {
    //     // Use std::copy to write the vector contents to the file
    //     std::copy(lookups.begin(), lookups.end(), std::ostream_iterator<int>(outputFile2, "\n"));

    //     // The file will be automatically closed when the ofstream object goes out of scope

    //     std::cout << "Lookup contents saved to " << lookupfile << std::endl;
    //     outputFile2.close();

    // } else {
    //     std::cerr << "Failed to open the file for writing." << std::endl;
    // }

    //Saving other data
    //=========

    // save_data(legitimate_data, "legitimate_data.txt");
    // save_data(poisoned_data, "poisoned_data.txt");
    // save_data(pois, "poisoned_data_pgm.txt");
    
}
