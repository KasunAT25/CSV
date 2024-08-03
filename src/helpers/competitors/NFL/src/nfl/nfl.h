#ifndef NFL_H
#define NFL_H

#include "../afli/afli.h"
#include "../afli/iterator.h"
#include "../benchmark/workload.h"
#include "../models/numerical_flow.h"
#include "../util/common.h"

namespace nfl {

template<typename KT, typename VT>
class NFL {
typedef std::pair<KT, VT> KVT;
typedef std::pair<KT, KVT> KKVT;
private:
  AFLI<KT, VT>* index_;
  uint32_t batch_size_;
  KVT* batch_kvs_;

  bool enable_flow_;
  NumericalFlow<KT, VT>* flow_;
  AFLI<KT, KVT>* tran_index_;
  KKVT* tran_kvs_;

  const float kConflictsDecay = 0.1;
  const uint32_t kMaxBatchSize = 4196;
  const float kSizeAmplification = 1.5;
  const float kTailPercent = 0.99;
public:
  explicit NFL(std::string weights_path, uint32_t batch_size) : batch_size_(batch_size) { 
    enable_flow_ = true;
    flow_ = new NumericalFlow<KT, VT>(weights_path, batch_size);
    index_ = nullptr;
    tran_index_ = nullptr;
    tran_kvs_ = nullptr;
    batch_kvs_ = nullptr;
  }

  ~NFL() {
    if (index_ != nullptr) {
      delete index_;
    }
    if (flow_ != nullptr) {
      delete flow_;
    }
    if (tran_index_ != nullptr) {
      delete tran_index_;
    }
    if (tran_kvs_ != nullptr) {
      delete[] tran_kvs_;
    }
    if (batch_kvs_ != nullptr) {
      delete[] batch_kvs_;
    }
  }

  void set_batch_size(uint32_t batch_size) {
    if (batch_size > batch_size_) {
      if (enable_flow_) {
        delete[] tran_kvs_;
        tran_kvs_ = new KKVT[batch_size];
      } else {
        delete[] batch_kvs_;
        batch_kvs_ = new KVT[batch_size];
      }
    }
    batch_size_ = batch_size;
  }

  uint32_t auto_switch(const KVT* kvs, uint32_t size, uint32_t aggregate_size=0) {
    tran_kvs_ = new KKVT[size];
    std::cout <<  "auto switch in 1" << std::endl;

    uint32_t origin_tail_conflicts = compute_tail_conflicts<KT, VT>(kvs, size, kSizeAmplification, kTailPercent);
    std::cout <<  "auto switch in 2" << std::endl;
    flow_->set_batch_size(kMaxBatchSize);

    
    //I ADDED THIS PART
    //=========================
    //WHY I ADDED
    // BECAUSE IT WAS GIVING 0 AS THE TRANSFORMED VALUE.
    // IT WAS NOT CALCULATING THE MEAN AND VAR BEFORE USING IT FOR THE TRANSFORMATION
    // SO WAS THE TRANSFORMATION CORRECT?????

    const std::pair<KT, VT>& vecRef = *kvs;

    long double sum =0;
    for(int i =0; i<size; i++){
      sum += static_cast<long double>(kvs[i].first);
    }
    // double sum = std::accumulate(vecRef.begin(), vecRef.end(), 0,
    //                           [](int partialSum, const auto& pair) {
    //                               return partialSum + pair.first;
    //                           });

    double mean = static_cast<long double>(sum) / (1.0*size);
     long double sumSquaredDiff = 0 ;
    // for(int i =0; i<size; i++){
    //   long double diff = static_cast<long double>(kvs[i].first) - mean;
    //   sumSquaredDiff += std::pow(diff, 2) ;
    // }
    //  double sumSquaredDiff = std::accumulate(vecRef.begin(), vecRef.end(), 0.0,
    //                                             [mean](long double partialSum, const auto& pair) {
    //                                                 long double diff = static_cast<long double>(pair.first) - mean;
    //                                                 return partialSum + std::pow(diff, 2);
    //                                             });
    
    double var = static_cast<double>(sumSquaredDiff) / (1.0*size);

    flow_->mean_ = mean;
    flow_->var_ = var;

    //std::cout <<  "auto switch in 2" << std::endl;

  //THIS IS WHERE THE ISSUE IS THE TRANSFORMATION GIVES 0
    flow_->transform(kvs, size, tran_kvs_);
    std::sort(tran_kvs_, tran_kvs_ + size, [](const KKVT& a, const KKVT& b) {
      return a.first < b.first;
    });
    
    std::cout <<  "auto switch in 3" << std::endl;

    //const size_t size_1 = 10001;

    // Print all values
    // for (size_t i = 0; i < size; ++i) {
    //     std::cout << "Value: " << tran_kvs_[i].first << " , ";
    // }

    uint32_t tran_tail_conflicts = compute_tail_conflicts<KT, KVT>(tran_kvs_, size, kSizeAmplification, kTailPercent);
    std::cout <<  "auto switch in 3.1" << std::endl;
    if (origin_tail_conflicts <= tran_tail_conflicts
      || origin_tail_conflicts - tran_tail_conflicts 
        < static_cast<uint32_t>(origin_tail_conflicts * kConflictsDecay)) {
      enable_flow_ = false;

      std::cout <<  "auto switch in 3.2" << std::endl;
      delete[] tran_kvs_;
      std::cout <<  "auto switch in 4" << std::endl;
      tran_kvs_ = nullptr;
      return origin_tail_conflicts;
    } else {
      enable_flow_ = true;
      std::cout <<  "auto switch in 4.2" << std::endl;
      return tran_tail_conflicts;
    }
    std::cout <<  "auto switch in 5" << std::endl;
  }

  void bulk_load(const KVT* kvs, uint32_t size, uint32_t tail_conflicts, uint32_t aggregate_size=0) {
    if (enable_flow_) {
      std::cout <<  "bulk flow" << std::endl;
      tran_index_ = new AFLI<KT, KVT>();
      tran_index_->bulk_load(tran_kvs_, size, tail_conflicts, aggregate_size);
      flow_->set_batch_size(batch_size_);
      delete tran_kvs_;
      tran_kvs_ = new KKVT[batch_size_];
      
    } else {
      std::cout <<  "bulk noflow" << std::endl;
      index_ = new AFLI<KT, VT>();
      index_->bulk_load(kvs, size, tail_conflicts, aggregate_size);
      batch_kvs_ = new KVT[batch_size_];   
        
    }
  }

  void transform(const KVT* kvs, uint32_t size) {
    if (enable_flow_) {
      flow_->transform(kvs, size, tran_kvs_);
    } else {
      std::memcpy(batch_kvs_, kvs, sizeof(KVT) * size);
    }
  }

  ResultIterator<KT, VT> find(uint32_t idx_in_batch) {
    if (enable_flow_) {
      auto it = tran_index_->find(tran_kvs_[idx_in_batch].first);
      if (!it.is_end()) {
        return {it.value_addr()};
      } else {
        return {};
      }
    } else {
      return index_->find(batch_kvs_[idx_in_batch].first);
    }
  }

  bool update(uint32_t idx_in_batch) {
    if (enable_flow_) {
      return tran_index_->update(tran_kvs_[idx_in_batch]);
    } else {
      return index_->update(batch_kvs_[idx_in_batch]);
    }
  }

  uint32_t remove(uint32_t idx_in_batch) {
    if (enable_flow_) {
      return tran_index_->remove(tran_kvs_[idx_in_batch].first);
    } else {
      return index_->remove(batch_kvs_[idx_in_batch].first);
    }
  }

  void insert(uint32_t idx_in_batch) {
    if (enable_flow_) {
      tran_index_->insert(tran_kvs_[idx_in_batch]);
    } else {
      index_->insert(batch_kvs_[idx_in_batch]);
    }
  }

  uint64_t model_size() {
    if (enable_flow_) {
      return tran_index_->model_size() + flow_->size();
    } else {
      return index_->model_size();
    }
  }

  uint64_t index_size() {
    if (enable_flow_) {
      return tran_index_->index_size() + flow_->size() 
            + sizeof(NFL<KT, VT>) + sizeof(KKVT) * batch_size_;
    } else {
      return index_->index_size() + sizeof(NFL<KT, VT>) + sizeof(KVT) * batch_size_;
    }
  }

  void print_stats() {
    if (enable_flow_) {
      tran_index_->print_stats();
    } else {
      index_->print_stats();
    }
  }
};

}

#endif