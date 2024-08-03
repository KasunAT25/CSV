// cd src/helpers/competitors/SALI/src/examples

// g++ -fopenmp -I /opt/intel/oneapi/tbb/2021.12/include example_multithread.cpp -std=c++17 -o example_multithread -ltbb
// g++ -fopenmp -I /opt/intel/oneapi/tbb/2021.12/include example_multithread.cpp -std=c++17 -o example_multithread -L /opt/intel/oneapi/tbb/2021.12/lib -ltbb

//g++ -fopenmp example_multithread.cpp -o example_multithread -ltbb

//./example_multithread 
// LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH ./example_multithread

// export LD_LIBRARY_PATH=/opt/intel/oneapi/tbb/2021.12/lib:$LD_LIBRARY_PATH

// sudo sh ./l_tbb_oneapi_p_2021.12.0.499_offline.sh --silent --eula accept --components intel.oneapi.lin.tbb.devel

#include <iostream>


#include <omp.h>
#include "../core/sali.h"

using namespace std;

int main() {
  sali::SALI<int, int> sali;

  int key_num = 1000;
  pair<int, int> *keys = new pair<int, int>[key_num];
  for (int i = 0; i < 1000; i++) {
    keys[i]={i,i};
  }
  sali.bulk_load(keys, 1000);

  omp_set_num_threads(12);

#pragma omp parallel for schedule(static, 12)
  for (int i = 1000; i < 2000; i++) {
    sali.insert(i,i);
  }
#pragma omp parallel for schedule(static, 12)
  for (int i = 0; i < 2000; i++) {
    std::cout<<"value at "<<i<<": "<<sali.at(i,i)<<std::endl;
  }

  return 0;
}
