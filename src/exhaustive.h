#pragma once

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <cmath>
#include <math.h>
#include <vector>
#include <string>

template<typename T>
std::vector<std::size_t> tag_sort(const std::vector<T>& v)
{
    std::vector<std::size_t> result(v.size());
    std::iota(std::begin(result), std::end(result), 0);
    // std::sort(std::begin(result), std::end(result),
    //           [&v](const auto & lhs, const auto & rhs)
    //           {
    //               return v[lhs] < v[rhs];
    //           }
    // );
    return result;
}
std::vector<KEY_TYPE> get_candidates_virtual(std::vector<KEY_TYPE> & keyset)
{
    std::vector<KEY_TYPE> result;

    KEY_TYPE start = keyset[0]+1;
    KEY_TYPE end = keyset[keyset.size()];
    int pointer = 0;

    for(int i = 0; i< keyset.size()-1 ; i++){
        
        while(start < keyset[i+1] && start > keyset[i]){
            result.push_back(start);
            start++;
        }
        start = keyset[i+1]+1;
    } 
    return result;
}

void get_subsets(const std::vector<KEY_TYPE>& keyset, int n, int index, std::vector<KEY_TYPE>& current, std::vector<std::vector<KEY_TYPE>>& result) {
    // Base case: if current subset size is n, add it to result
    if (current.size() == n) {
        result.push_back(current);
        return;
    }
    
    // If we've processed all elements or remaining elements are not enough, return
    if (index >= keyset.size() || current.size() + (keyset.size() - index) < n) {
        return;
    }
    
    // Include nums[index] in the current subset and recurse
    current.push_back(keyset[index]);
    get_subsets(keyset, n, index + 1, current, result);
    
    // Backtrack: Exclude nums[index] from the current subset and recurse
    current.pop_back();
    get_subsets(keyset, n, index + 1, current, result);
}

std::vector<std::vector<KEY_TYPE>> get_subsets_n(const std::vector<KEY_TYPE>& keyset, int n) {
    std::vector<std::vector<KEY_TYPE>> result;
    std::vector<KEY_TYPE> current;
    
    get_subsets(keyset, n, 0, current, result);
    
    return result;
}

std::vector<std::vector<KEY_TYPE>> get_subsets_n2(const std::vector<KEY_TYPE>& keyset, int k) {
    std::vector<std::vector<KEY_TYPE>> subsets;
    int n = keyset.size();
    if (k > n) return subsets;

    std::vector<int> indices(k);
    for (int i = 0; i < k; ++i) {
        indices[i] = i;
    }

    while (true) {
        // Collect the current subset
        std::vector<KEY_TYPE> subset;
        for (int index : indices) {
            subset.push_back(keyset[index]);
        }
        subsets.push_back(subset);

        // Find the rightmost index that can be incremented
        int i;
        for (i = k - 1; i >= 0; --i) {
            if (indices[i] != i + n - k) {
                break;
            }
        }

        // If no such index exists, we are done
        if (i < 0) break;

        // Increment this index and adjust subsequent indices
        ++indices[i];
        for (int j = i + 1; j < k; ++j) {
            indices[j] = indices[j - 1] + 1;
        }
    }

    return subsets;
}

std::vector<KEY_TYPE> obtain_poisoning_keys(std::vector<KEY_TYPE> & keyset, std::vector<std::vector<KEY_TYPE>> subsets) {
    std::set<KEY_TYPE> poisoning_keys;

    long double min_loss = 0;
    std::vector <KEY_TYPE> best_subset;

     for(int i = 0; i < subsets.size();i++){

        std::vector <KEY_TYPE> keyset2;
        std::vector <KEY_TYPE> subset = subsets[i];
        std::merge(keyset.begin(), keyset.end(), subset.begin(), subset.end(), std::back_inserter(keyset2));

        int size = keyset2.size();

        std::vector<size_t> rankset = tag_sort(keyset2);

        long double sumX2 = 0;
        double sumY2 = 0;
        long double sumX = 0;
        double sumY = 0;
        long double sumXY = 0;

        for(int i =0; i < keyset2.size(); i++){
            sumX2 += static_cast<long double>(keyset2[i]) * keyset2[i];
            sumY2 += rankset[i]*rankset[i];

            sumX += static_cast<long double>(keyset2[i]);
            sumY += rankset[i];

            sumXY += static_cast<long double>(keyset2[i])*rankset[i];
        }

        long double avgX = sumX/(1.0*size);
        double avgY = sumY/(1.0*size);
        long double avgXY = sumXY/(1.0*size);

        long double avgX2 = sumX2/(1.0*size);
        double avgY2 = sumY2/(1.0*size);

        long double vX = avgX2 - pow(avgX,2);
        double vY = avgY2 - pow(avgY,2);
        long double vXY = avgXY - avgX*avgY;


       
        long double loss = (-pow(vXY,2)/vX + vY)*size;
        if(i == 0){
            min_loss = loss;
            best_subset = subset;
        }
        else if(loss < min_loss ){   
            min_loss = loss;
            best_subset = subset;
        }
            //std::cout<< (subsets[i][j]) << " ";
        
        //std::cout<< subset[0] << " loss " << (loss)<< std::endl;
    }
    std::cout<< "loss " << (min_loss)<< std::endl;

    return best_subset;
}

std::vector<KEY_TYPE> perform_poisoning(std::vector<KEY_TYPE> & legit_keys, double poisoning_threshold) {

    int size = legit_keys.size();
    std::vector<size_t> legit_ranks = tag_sort(legit_keys);

    std::vector<KEY_TYPE> candidates = get_candidates_virtual(legit_keys);
    std::cout<< "Candidate keys generated"<< std::endl;
    // for(int i = 0; i < candidates.size();i++){
    //     std::cout<< (candidates[i]) << " ";
    // }
    //  std::cout<< "---------" << std::endl;

    

    int num_virtual = static_cast<int>(poisoning_threshold*size);

    //std::cout << num_virtual << std::endl;

    std::vector<std::vector<KEY_TYPE>> subsets = get_subsets_n2(candidates, num_virtual);
    std::cout<< "Subsets generated"<< std::endl;

    
    // for(int i = 0; i < subsets.size();i++){

    //     for(int j = 0; j < subsets[i].size();j++){
    //         std::cout<< (subsets[i][j]) << " ";
    //     }
    //     std::cout<< "---------" << std::endl;
    // }

    

    std::vector<KEY_TYPE> poisoning_keys = obtain_poisoning_keys(legit_keys, subsets);
    std::cout<< "Poisoned"<< std::endl;

    //std::cout<< "loss " << (poisoning_keys.size())<< std::endl;

    for(int j = 0; j < poisoning_keys.size();j++){
            std::cout<< (poisoning_keys[j]) << " ";
     }

    /*
     * Merge poisoning keys with legitimate keys to generated poisoned keyset
     */
        
    std::vector <KEY_TYPE> poisoned_keyset;
    std::merge(legit_keys.begin(), legit_keys.end(), poisoning_keys.begin(), poisoning_keys.end(), std::back_inserter(poisoned_keyset));
    std::sort(poisoned_keyset.begin(), poisoned_keyset.end());

 

    return poisoned_keyset;

}

