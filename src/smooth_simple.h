#pragma once

#include <set>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <utility>


template<typename T>
std::vector<std::size_t> tag_sort(const std::vector<T>& v)
{
    std::vector<std::size_t> result(v.size());
    std::iota(std::begin(result), std::end(result), 0);
 
    return result;
}


size_t compute_rank_for_endpoint(KEY_TYPE & endpoint, std::vector<KEY_TYPE> & keyset){
   
        //computed_rank_for_endpoint = (rank);
        auto it = std::lower_bound(keyset.begin(), keyset.end(), endpoint);
        int rank = std::distance(keyset.begin(), it);
    //}
    return rank;
}

// Function to find all with the largest gap in a sorted vector of integers
std::vector<std::pair<KEY_TYPE, KEY_TYPE>> find_largest_gaps(std::vector<KEY_TYPE>& keyset) {
    std::vector<std::pair<KEY_TYPE, KEY_TYPE>> largestGapPairs;
    std::vector<KEY_TYPE> gaps;

    std::vector<std::pair<std::pair<KEY_TYPE, KEY_TYPE>, KEY_TYPE>> gaps_pairs;

    KEY_TYPE largestGap = 0;
    KEY_TYPE smallestGap = keyset[keyset.size()-1];
    for (size_t i = 1; i < keyset.size(); ++i) {
        KEY_TYPE currentGap = keyset[i] - keyset[i - 1];
        gaps.push_back(currentGap);


        gaps_pairs.push_back({{keyset[i-1], keyset[i]}, currentGap});
        
    }
    std::sort(gaps_pairs.begin(), gaps_pairs.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    size_t topCount = std::min(static_cast<size_t>(gaps_pairs.size()), static_cast<size_t>(25));
    
    if(topCount > gaps_pairs.size()){
        topCount = gaps_pairs.size();
    }

    for (size_t i = 0; i < topCount; ++i) {
        largestGapPairs.push_back(gaps_pairs[i].first);
    }

    return largestGapPairs;
}

/*
    Extract non-occupied keys for a given sequence of legitimate and poisoning keys
*/
void partition_non_occupied_keys_gaps(std::vector<std::pair<KEY_TYPE, KEY_TYPE>> & K,std::vector<KEY_TYPE>& S,std::vector<std::pair<KEY_TYPE, KEY_TYPE>>& D) {
    
    std::vector<KEY_TYPE> endpoints;

    bool prev = false;


    //Now that I have the keys
    //Find which ones I have to calculate Loss and derivative
        for(int i =0; i< K.size();i++){
            //Get current and next endpoint
            KEY_TYPE start = K[i].first;
            KEY_TYPE end = K[i].second;

            //Change 6 accordinigly
            //Check if the difference is greater than 3
            if(end-start<3){
                
                //Then calculatet the loss for all the points in between (the two endpoints)
                for(KEY_TYPE j = K[i].first+1; j< end; j++){
                    S.push_back(j);

                }
            }
            //Else add to calculate derivative
            else{
                //Get as pairs
                D.push_back(std::make_pair(K[i].first+1,K[i].second-1));
            }
        }

}

void partition_non_occupied_keys(std::vector<KEY_TYPE> & K,std::vector<KEY_TYPE>& S,std::vector<std::pair<KEY_TYPE, KEY_TYPE>>& D) {
    
    std::vector<KEY_TYPE> endpoints;
    bool prev = false;

    //Now that I have the keys
    //Find which ones I have to calculate Loss and derivative
    if(K.size()>2){

        for(int i =0; i< K.size()-1;i++){
            //Get current and next endpoint
            KEY_TYPE start = K[i];
            KEY_TYPE end = K[i+1];

            //Change 6 accordinigly
            //Check if the difference is greater than 3
            if(end-start<3){
                
                //Then calculatet the loss for all the points in between (the two endpoints)
                for(KEY_TYPE j = K[i]+1; j< end; j++){
                    S.push_back(j);

                }
            }
            //Else add to calculate derivative
            else{
                //Get as pairs
                D.push_back(std::make_pair(K[i]+1,K[i+1]-1));
            }
        }

    }

}
void cal_derv(std::vector<KEY_TYPE> & current_keyset2,std::vector<size_t>& rankset,
        std::map<size_t, std::pair< long double,  long double>>& dL, 
        std::vector<KEY_TYPE>& S,
        std::vector<std::pair<KEY_TYPE, KEY_TYPE>>& D,
        std::map<std::pair<KEY_TYPE, KEY_TYPE>, std::pair<  double,   double>>& M , long double &sumXo,
        double &sumYo,
        long double &sumXYo,
        long double &sumX2) {

    //I can probably calculate rankset and T here itself


    //Calculate the basics needed

    int n = current_keyset2.size();

    double dL1 = 0.0;
    double avgXo = sumXo/(1.0*(n));
    int yp;
    KEY_TYPE xp;

    long double avgXp;
    double avgYp;
    double sumY;
    double avgY;
    long double sumXY;
    long double A;
    long double B;
    long double m;
    long double b;
    long double dm;
    long double db;
    long double sign;

    int yp_prev;

    //For all the derivatives
    for(int i = 0; i < D.size(); i++) {
    
    //check if it is the first or second
    //Check if it is first of second
    for(int j =0; j<2;j++){
        //Add both and check
        if(j==0){
               xp = D[i].first;
        }
        else{
            xp = D[i].second;
        }
        yp = compute_rank_for_endpoint(xp, current_keyset2);

        avgXp = (sumXo + static_cast<long double>(xp))/(1.0*(n+1));

        avgYp = (sumYo + n)/(1.0*(n+1));

        //calculate new ranks

        sumY = sumYo + n -yp;

        avgY = sumY/(1.0*n);

        //double check this yp (if it is yp or yp -1 )
        //also check current keyset or current keyset2

        if(i == 0){
            sumXY = sumXYo;
            for(int k = yp; k < n; k++ ){

                sumXY += static_cast<long double>(current_keyset2[k]);
            }

        }
        else{  
            for(int k = yp_prev; k < yp; k++ ){

                sumXY -= static_cast<long double>(current_keyset2[k]);
            }

        }

        A = (n+1)*(sumX2+static_cast<long double>(xp)*static_cast<long double>(xp)) - pow(((n+1)*avgXp),2);

        B = (n+1)*(sumXY+static_cast<long double>(xp)*yp) - (pow((n+1),2)*avgXp*avgYp);

        m = ((sumXY+static_cast<long double>(xp)*yp)-(n+1)*avgXp*avgYp)/(sumX2+static_cast<long double>(xp)*static_cast<long double>(xp)-(n+1)*avgXp*avgXp);

        b = avgYp - m*avgXp;

        dm = (A*(n*(yp-avgY))-B*(2*n*(static_cast<long double>(xp)-avgXo)))/(A*A);

        db = -(m+(n+1)*avgXp*dm)/(1.0*(n+1));

        //Calculate the derivative wiht r xp
        dL1 = 2*(dm*(m*sumX2+n*b*avgXo-sumXY) + n*db*(m*avgXo+b-avgY) + (m*static_cast<long double>(xp)+b-yp)*(dm*static_cast<long double>(xp)+m+db));

        //Add them to the corresponding pair
        //=======================================
        
        if(j==0){
            dL[i].first = dL1;
           

        }
        else{
            dL[i].second = dL1;
        

        }
       
    yp_prev = yp;
        
    }

       

        //Check the sign
        sign = dL[i].first*dL[i].second;

        
        //Maybe can be back inserter
        //Can probablt avoid duplicates if I only save the indexes and use the dL to access them
        //
        //==============================
        //If opposite signs add to find min (M)
         if(sign<0){

            M.insert({std::make_pair(D[i].first, D[i].second), std::make_pair(dL[i].first, dL[i].second)});
        }
        //Calculate L here itself No need to add to S
        //=========================
        else{
            S.push_back(D[i].first);
            S.push_back(D[i].second);
        }
        }
        

}

void cal_loss(std::vector<KEY_TYPE> & current_keyset2,std::vector<size_t>& rankset,
        std::map<size_t,  double>& L, 
        std::vector<KEY_TYPE>& S, bool L_p,   double & L_prev, long double &sumXo,
        double &sumYo,
        long double &sumXYo,
        long double &sumX2,
        double &sumY2o) {




    int n = current_keyset2.size();

    long double avgXo = sumXo/(1.0*(n));


    if(L_p){

        L_prev = (-pow(sumXYo/(1.0*n) - (sumXo/(1.0*n))*(sumYo/(1.0*n)),2)/((sumX2/(1.0*n))-pow(sumXo/(1.0*n),2)) + sumY2o/(1.0*n) - pow(sumYo/(1.0*n),2))*n;
        
        if(L_prev < 0.000000001){
            L_prev = 0;
        }
        if(L_prev < 0){
            std::cout << "MINUS "  << std::endl;
        }
    }

    int yp;
    KEY_TYPE xp;
    long double avgXp;
    double avgYp;
    double sumY;
    double avgY;
    long double sumXY;

    long double m;
    long double b;

    int yp_prev;

    //For all the derivatives
    for(int i = 0; i < S.size(); i++) {

        xp = S[i];
 
        yp = compute_rank_for_endpoint(xp, current_keyset2);

        avgXp = (sumXo + static_cast<long double>(xp))/(1.0*(n+1));

        avgYp = (sumYo + n)/(1.0*(n+1));

        sumY = sumYo + n -yp;

        avgY = sumY/(1.0*n);

        if(i == 0){
            sumXY = sumXYo;
            for(int k = yp; k < n; k++ ){

                sumXY += static_cast<long double>(current_keyset2[k]);
            }
        }
        else{  
            for(int k = yp_prev; k < yp; k++ ){

                sumXY -= static_cast<long double>(current_keyset2[k]);
            }

        }

        long double bottom_m = (static_cast<long double>(sumX2)+static_cast<long double>(xp)*static_cast<long double>(xp)-(n+1)*static_cast<long double>(avgXp)*static_cast<long double>(avgXp));
        m = ((sumXY+static_cast<long double>(xp)*yp)-(n+1)*avgXp*avgYp)/bottom_m;

        b = avgYp - m*avgXp;

        //Calculate the loss 
        
         L[i] = (m*m*sumX2 +2*n*m*b*avgXo - 2*m*sumXY+n*b*b - 2*n*b*avgY + sumY2o + (n*n) - yp*yp + pow(m*static_cast<long double>(xp)+b-yp,2));
        
        if(L[i] < 0.000000001){
            L[i] = 0;
        }

        yp_prev = yp;
        }

}

void minimum_points(std::map<std::pair<KEY_TYPE, KEY_TYPE>, std::pair<  double,   double>>& M, std::vector<KEY_TYPE>& S){
        double val;
        KEY_TYPE xp;
       for (const auto& entry : M) {
            val = (entry.second.second*entry.first.first-entry.second.first*entry.first.second)/(1.0*entry.second.second-entry.second.first);
            xp = static_cast<uint64_t>(std::round(val));
            S.push_back(xp);
         }

}
/*
    Implementation of the greedy poisoning attack on regression models as described in Kornaropoulos et al. ("The Price of Tailoring the Index to Your Data: Poisoning Attacks on Learned Index Structures")
 */
std::set<KEY_TYPE> obtain_poisoning_keys(double poisoning_threshold, std::vector<KEY_TYPE> & keyset, std::vector<size_t> &rankset, bool max = false) {
    
    // Total number of elements
    
    int n = keyset.size();
    double L_prev = 0;
    double original_L_prev = 0;

    // Number of poisoning keys P
    int P = int(poisoning_threshold * n);

    std::set<KEY_TYPE> poisoning_keys;


    long double sumXo = 0;
    
    double sumYo = 0;
    
    long double sumX2 = 0;

    double sumY2o = 0;

    long double sumXYo = 0;

     for(int i =0; i < keyset.size(); i++){
        sumX2 += static_cast<long double>(keyset[i]) * keyset[i];
        sumXo += static_cast<long double>(keyset[i]);
        sumXYo += static_cast<long double>(keyset[i])*rankset[i];

        sumYo += rankset[i];
        sumY2o += rankset[i]*rankset[i];

    }


    KEY_TYPE xp_prev;
    size_t rank_prev;

    while(poisoning_keys.size() < P){

        std::vector <KEY_TYPE> keyset2;
        std::merge(keyset.begin(), keyset.end(), poisoning_keys.begin(), poisoning_keys.end(), std::back_inserter(keyset2));
        n = keyset2.size();

       
        std::vector<size_t> rankset = tag_sort(keyset2);

        
        if(poisoning_keys.size() > 0){
            sumXo = sumXo + static_cast<long double>(xp_prev);
            sumYo = sumYo + (n-1);
            sumXYo += static_cast<long double>(xp_prev)*rank_prev;


            //This maybe a time issue
            
            
            for(int i = (rank_prev+1); i<keyset2.size();i++){
                sumXYo+=keyset2[i];
            }
            

            sumX2 = sumX2 + static_cast<long double>(xp_prev)*static_cast<long double>(xp_prev);
            sumY2o = sumY2o + (n-1)*(n-1) ;
        }
        

        // Partition the non-occupied keys into subsequences such that each subsequence consists of consecutive non-occupied keys;
        // Extract the endpoints of each subsequence and sort them to construct the new sequence of endpoints S(i), where i <= 2(n + j);

        //Calculate the endpoints and which ones to add to caluclate loss or derv

        // To store the calculate loss
        
        //===========================================
        
        

        std::vector<KEY_TYPE> S;
        std::vector<KEY_TYPE> S_temp;

        // To store the calculate derv (stores the two potential endpoints)
        std::vector<std::pair<KEY_TYPE, KEY_TYPE>> D;

        std::map<std::pair<KEY_TYPE, KEY_TYPE>, std::pair<  double,   double>> M;

        std::vector<std::pair<KEY_TYPE, KEY_TYPE>> G;
       

        if(keyset2.size()> 25){
             G = find_largest_gaps(keyset2);
        
            partition_non_occupied_keys_gaps(G,S,D);
        }
        else{
            partition_non_occupied_keys(keyset2,S,D);
        }

        // double temp_dur;
        // auto start_temp = std::chrono::high_resolution_clock::now();
        // auto end_temp = std::chrono::high_resolution_clock::now();
        // temp_dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end_temp - start_temp).count();
        // std::cout << "temp duration: " << temp_dur << std::endl;

        std::map<size_t, std::pair< long double,  long double>> dL;
        std::map<size_t,   double> L;

        //auto derv_start = std::chrono::high_resolution_clock::now();
        //===========================================
         
         
        
        cal_derv(keyset2,rankset,
        dL, 
        S,
        D,
        M,
        sumXo,sumYo,sumXYo,sumX2);

        //auto derv_end = std::chrono::high_resolution_clock::now();

        //auto derv_dur = std::chrono::duration_cast<std::chrono::nanoseconds>(derv_end - derv_start).count();
        //  std::cout << "derv duration: " << derv_dur << std::endl;
        // std::cout << "D: " << D.size() << std::endl;
        // std::cout << "derv duration: " << derv_dur/(2.0*D.size()) << std::endl;
        //===========================================

        //std::cout << "Deriv: "<< std::endl;

        //Now calculate the minimum points
       
        minimum_points(M, S);

        //Now calculate the L for all in S
        std::sort(S.begin(), S.end());

        //===========================================
        


        //Calculate the L_prev for the first one
        
        if(poisoning_keys.size()==0){
            cal_loss(keyset2,rankset,
        L, 
        S,true, L_prev,
        sumXo,sumYo,sumXYo,sumX2,sumY2o);

        original_L_prev = L_prev;
        }

        else{
            cal_loss(keyset2,rankset,
        L, 
        S,false, L_prev,
        sumXo,sumYo,sumXYo,sumX2,sumY2o);
        }

        // auto loss_end = std::chrono::high_resolution_clock::now();

        //auto loss_dur = std::chrono::duration_cast<std::chrono::nanoseconds>(loss_end - loss_start).count();
        // std::cout << "loss duration: " << loss_dur << std::endl;
        // std::cout << "S: " << S.size() << std::endl;
        // std::cout << "loss duration: " << loss_dur/(1.0*S.size()) << std::endl;
        //===========================================
       

        int optimal_key_index;

         auto maxElementIterator = std::min_element(L.begin(), L.end(),[](const std::pair<size_t,  long double>& a, const std::pair<size_t,  long double>& b) {
            return a.second < b.second; // Compare by the second element (the double value)
        });

        optimal_key_index = maxElementIterator->first;
          
        if(maxElementIterator->second >= L_prev || S.size()==0 || L_prev == 0){
       
            break;
        }

    
        poisoning_keys.insert(S[optimal_key_index]);
         L_prev = maxElementIterator->second;
         
        xp_prev = S[optimal_key_index];
        rank_prev = compute_rank_for_endpoint(xp_prev, keyset2);

        
        
        if(poisoning_keys.size()== P){
            
            break;
        }
        
    }

    // auto poisoning_end = std::chrono::high_resolution_clock::now();

    // auto poisoning_dur = std::chrono::duration_cast<std::chrono::microseconds>(poisoning_end - poisoning_start).count();
    // std::cout << "poisoning duration: " << poisoning_dur << std::endl;
    //std::cout << " O_L " << original_L_prev << std::endl;
    //std::cout << " L " << L_prev << std::endl;
    //std::cout << "poisoning duration: " << poisoning_dur/(1.0*poisoning_keys.size()) << std::endl;
    //===========================================
    return poisoning_keys;
}

/*
    Generate poisoning keys based on legitimate keyset (up to specified poisoning threshold)
*/

std::vector<KEY_TYPE> perform_poisoning(std::vector<KEY_TYPE> & legit_keys, double poisoning_threshold,bool max =true) {

    std::vector<size_t> legit_ranks = tag_sort(legit_keys);


    std::set<KEY_TYPE> poisoning_keys = obtain_poisoning_keys(poisoning_threshold, legit_keys, legit_ranks,max);

    /*
     * Merge poisoning keys with legitimate keys to generated poisoned keyset
     */
    
        
    std::vector <KEY_TYPE> poisoned_keyset;
    std::merge(legit_keys.begin(), legit_keys.end(), poisoning_keys.begin(), poisoning_keys.end(), std::back_inserter(poisoned_keyset));
    std::sort(poisoned_keyset.begin(), poisoned_keyset.end());


    return poisoned_keyset;

}
