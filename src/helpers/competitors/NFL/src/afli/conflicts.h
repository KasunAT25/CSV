#ifndef CONFLICTS_H
#define CONFLICTS_H

#include "../util/common.h"
#include "../models/linear_model.h"

namespace nfl {

struct ConflictsInfo {
  uint32_t* conflicts_;
  uint32_t* positions_;
  uint32_t num_conflicts_;
  uint32_t max_size_;
  uint32_t size_;

  ConflictsInfo(uint32_t size, uint32_t max_size) 
    : max_size_(max_size), num_conflicts_(0) {
    size_ = size;
    conflicts_ = new uint32_t[size];
    positions_ = new uint32_t[size];
  }

  ~ConflictsInfo() {
    delete[] conflicts_;
    delete[] positions_;
  }

  void add_conflict(uint32_t position, uint32_t conflict) {
    conflicts_[num_conflicts_] = conflict;
    positions_[num_conflicts_] = position;
    num_conflicts_ ++;
  }
};

template<typename KT, typename VT>
ConflictsInfo* build_linear_model(const std::pair<KT, VT>* kvs, uint32_t size,
                                  LinearModel<KT>*& model, 
                                  double size_amp) {

                                     //std::cout <<  "LI 0" << std::endl;
  if (model != nullptr) {
    model->slope_ = model->intercept_ = 0;
  } else {
    model = new LinearModel<KT>();
  }
  // OPT: Find a linear regression model that has the minimum conflict degree
  // The heuristics below is a simple method that scales the positions
  KT min_key = kvs[0].first;
  KT max_key = kvs[size - 1].first;
  KT key_space = max_key - min_key;
  if (compare(min_key, max_key)) {
    delete model;
    model = nullptr;
    //std::cout <<  min_key << " " << max_key<<  std::endl;
    //std::cout <<  "LI 1" << std::endl;
    return nullptr;
  }
  uint32_t max_size = static_cast<uint32_t>(size * size_amp);
  LinearModelBuilder<KT> builder;
  uint32_t j = 0;
  for (uint32_t i = 0; i < size; ++ i) {
    KT key = kvs[i].first;
    // double y = max_size * (key - min_key) / key_space;
    double y = i;
    builder.add(key, y);
  }
  builder.build(model);
  if (compare(model->slope_, 0.)) {
    //std::cout <<  "LI 2" << std::endl;
    // Fail to build a linear model
    delete model;
    model = nullptr;
    return nullptr;
  } else {
    model->intercept_ = -model->slope_ * (min_key) + 0.5;
    int64_t predicted_size = model->predict(max_key) + 1;
    if (predicted_size > 1) {
      max_size = std::min(predicted_size, static_cast<int64_t>(max_size));
    }

    uint32_t first_pos = std::min(std::max(model->predict(min_key), static_cast<int64_t>(0)), 
                                  static_cast<int64_t>(max_size - 1));
    uint32_t last_pos = std::min(std::max(model->predict(max_key), static_cast<int64_t>(0)), 
                                  static_cast<int64_t>(max_size - 1));
    if (last_pos == first_pos) {
      // Model fails to predict since all predicted positions are rounded to the 
      // same one
      model->slope_ = size / (max_key - min_key);
      model->intercept_ = -model->slope_ * (min_key) + 0.5;
    }
    ConflictsInfo* ci = new ConflictsInfo(size, max_size);
    uint32_t p_last = first_pos;
    uint32_t conflict = 1;
    for (uint32_t i = 1; i < size; ++ i) {
      uint32_t p = std::min(std::max(model->predict(kvs[i].first), static_cast<int64_t>(0)), 
                                    static_cast<int64_t>(max_size - 1));
      if (p == p_last) {
        conflict ++;
      } else {
        ci->add_conflict(p_last, conflict);
        p_last = p;
        conflict = 1;
      }
    }
    if (conflict > 0) {
      ci->add_conflict(p_last, conflict);
    }
   // std::cout <<  "LI 1,0" << std::endl;
    
    return ci;
  }
}

template<typename KT, typename VT>
uint32_t compute_tail_conflicts(const std::pair<KT, VT>* kvs, uint32_t size, 
                                double size_amp, float kTailPercent=0.99) {
  // The input keys should be ordered
  LinearModel<KT>* model = new LinearModel<KT>();
  ConflictsInfo* ci = build_linear_model<KT, VT>(kvs, size, model, size_amp);

  std::cout <<  "conflict 1" << std::endl;
  delete model;
  std::cout <<  ci->num_conflicts_ << std::endl;
  std::cout <<  "conflict 2" << std::endl;

  if (ci->num_conflicts_ == 0) {
     std::cout <<  "conflict 3.1" << std::endl;
    delete ci;
    std::cout <<  "conflict 3" << std::endl;
    return 0;
  } else {
    std::cout <<  "conflict 4,0" << std::endl;
    std::sort(ci->conflicts_, ci->conflicts_ + ci->num_conflicts_);
    std::cout <<  "conflict 4,1" << std::endl;
    uint32_t tail_conflicts = ci->conflicts_[std::max(0, int(ci->num_conflicts_ * kTailPercent) - 1)];
    std::cout <<  "conflict 4,2" << std::endl;
    delete ci;
    std::cout <<  "conflict 4" << std::endl;
    return tail_conflicts - 1;
  }
}

}
#endif