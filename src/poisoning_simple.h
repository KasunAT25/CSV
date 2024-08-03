#pragma once

#include <set>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <utility>

//typedef u_int64_t key_type;
//PLAN
// Get the gaps
// Get the max gaps
// Then run the derivates and get the S



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
        

        // if (currentGap > largestGap) {
        //     largestGap = currentGap;
        // }
        // if (currentGap < smallestGap) {
        //     smallestGap = currentGap;
        // }
    }
    std::sort(gaps_pairs.begin(), gaps_pairs.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    //size_t topCount = std::max(static_cast<size_t>(0.2 * gaps_pairs.size()), static_cast<size_t>(25));
    size_t topCount = std::min(static_cast<size_t>(gaps_pairs.size()), static_cast<size_t>(25));
    //size_t topCount = 50;
    
    if(topCount > gaps_pairs.size()){
        topCount = gaps_pairs.size();
    }

    for (size_t i = 0; i < topCount; ++i) {
        largestGapPairs.push_back(gaps_pairs[i].first);
    }

    // std::sort(gaps.begin(),gaps.end());
    // int size = gaps.size();
    // KEY_TYPE val = gaps[static_cast<int>(size*0.5)];
    // //std::cout<< largestGap <<std::endl;
    
    // for (size_t i = 1; i < keyset.size(); ++i) {
    //     KEY_TYPE currentGap = keyset[i] - keyset[i - 1];

    //     if (currentGap >= val ) {
    //        largestGapPairs.push_back({keyset[i - 1], keyset[i]}); 
    //        //std::cout<< keyset[i - 1] << " | " << keyset[i] <<std::endl;
    //     }
       
    // }

    return largestGapPairs;
}

//Double check this
// std::vector<std::pair<KEY_TYPE, KEY_TYPE>> find_largest_gaps(std::vector<KEY_TYPE>& values) {
//     int numGaps = 1;

//     std::vector<std::pair<KEY_TYPE, KEY_TYPE>> pairs(values.size() - 1);
//     std::transform(values.begin() + 1, values.end(), values.begin(), pairs.begin(), [](const KEY_TYPE& current, const KEY_TYPE& prev) {
//         return std::make_pair(prev, current);
//     });

//     std::partial_sort(pairs.begin(), pairs.begin() + numGaps, pairs.end(), [](const auto& a, const auto& b) {
//         return a.second - a.first > b.second - b.first;
//     });

//     return std::vector<std::pair<KEY_TYPE, KEY_TYPE>>(pairs.begin(), pairs.begin() + numGaps);
// }
/*
    Extract non-occupied keys for a given sequence of legitimate and poisoning keys
*/
void partition_non_occupied_keys_gaps(std::vector<std::pair<KEY_TYPE, KEY_TYPE>> & K,std::vector<KEY_TYPE>& S,std::vector<std::pair<KEY_TYPE, KEY_TYPE>>& D) {

    //std::cout << "Partitioning " << std::endl;
    
    std::vector<KEY_TYPE> endpoints;

    //endpoints.push_back(K[0]);
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
                //D.push_back(endpoints[i+1]-1);

            }
        }

}

void partition_non_occupied_keys(std::vector<KEY_TYPE> & K,std::vector<KEY_TYPE>& S,std::vector<std::pair<KEY_TYPE, KEY_TYPE>>& D) {

    //std::cout << "Partitioning " << std::endl;
    
    std::vector<KEY_TYPE> endpoints;

    //endpoints.push_back(K[0]);
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
                //D.push_back(endpoints[i+1]-1);

            }
        }

    }

    //return pointsBetween;

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

    //sumXo = std::accumulate(current_keyset2.begin(), current_keyset2.end(), 0.0);
    double avgXo = sumXo/(1.0*(n));
   

    //sumYo = std::accumulate(rankset.begin(), rankset.end(), 0.0);
    //sumXYo = std::inner_product(current_keyset2.begin(), current_keyset2.end(), rankset.begin(), 0.0);
    //sumX2 = std::inner_product(current_keyset2.begin(), current_keyset2.end(), current_keyset2.begin(), 0.0);

    // static_cast<long double>(sumYo);
    // std::cout << "sumXo " << typeid(sumXo).name() << std::endl;
    // std::cout << "sumYo " << typeid(sumYo).name() << std::endl;
    // std::cout << "sumXYo " << typeid( static_cast<long double>(sumXYo)).name() << std::endl;
    // std::cout << "sumY2o " << typeid(sumXo).name() << std::endl;
    
    


    //double sumY;
    

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
    //std::cout << "D: " << D.size() << std::endl;
    for(int i = 0; i < D.size(); i++) {
    
    //check if it is the first or second
   // for(int k =0; k<2;k++){
        

    //Check if it is first of second
    for(int j =0; j<2;j++){

        //std::vector<key_type> current_keyset = current_keyset2;
        //Add both and check
        if(j==0){
               // if(k==0)
               xp = D[i].first;
                //current_keyset.push_back(D[i].first);
                    
                // else
                //     current_keyset.push_back(D[i].first+1);

        }
        else{
            xp = D[i].second;

               // if(k==0)
                    // current_keyset.push_back(D[i].second);
                // else
                //     current_keyset.push_back(D[i].second-1);
        }
        yp = compute_rank_for_endpoint(xp, current_keyset2);
        //std::cout << "yp: " << yp << std::endl;

        avgXp = (sumXo + static_cast<long double>(xp))/(1.0*(n+1));

        avgYp = (sumYo + n)/(1.0*(n+1));

        //calculate new ranks
        //current_keyset.push_back(xp);
        //std::vector<size_t> ranksetp = tag_sort(current_keyset);

        sumY = sumYo + n -yp;

        avgY = sumY/(1.0*n);

        //double check this yp (if it is yp or yp -1 )
        //also check current keyset or current keyset2
        
        // //sumXY = sumXYo + std::accumulate(current_keyset2.begin()+yp, current_keyset2.end(), 0.0);

        // //THIS IS THE ISSUE PART just use that maths arithmatic summation or something
        // sumXY = sumXYo;
        // for(int j = yp; j < n; j++ ){
        //     //std::cout<< "j " ;
        //     //sumXY += static_cast<long double>(current_keyset2[j]);
        //     sumXY += static_cast<long double>(current_keyset2[j]);
        //     //std::cout<< current_keyset2[j] << std::endl;
        // }

        if(i == 0){
            sumXY = sumXYo;
            for(int k = yp; k < n; k++ ){
            //std::cout<< "j " ;
            //sumXY += static_cast<long double>(current_keyset2[j]);
                sumXY += static_cast<long double>(current_keyset2[k]);
            //std::cout<< current_keyset2[j] << std::endl;
            }
            //sumXY_prev = sumXY;
        }
        else{  
            for(int k = yp_prev; k < yp; k++ ){
            //std::cout<< "j " ;
            //sumXY += static_cast<long double>(current_keyset2[j]);
                sumXY -= static_cast<long double>(current_keyset2[k]);
            //std::cout<< current_keyset2[j] << std::endl;
            }

            //sumXY = sumXY_prev;
        }

        //double A = (n+1)*(sumX2+xp*xp) - pow((n*avgXo + xp),2);
        A = (n+1)*(sumX2+static_cast<long double>(xp)*static_cast<long double>(xp)) - pow(((n+1)*avgXp),2);

        B = (n+1)*(sumXY+static_cast<long double>(xp)*yp) - (pow((n+1),2)*avgXp*avgYp);

        m = ((sumXY+static_cast<long double>(xp)*yp)-(n+1)*avgXp*avgYp)/(sumX2+static_cast<long double>(xp)*static_cast<long double>(xp)-(n+1)*avgXp*avgXp);

        b = avgYp - m*avgXp;

        dm = (A*(n*(yp-avgY))-B*(2*n*(static_cast<long double>(xp)-avgXo)))/(A*A);

        db = -(m+(n+1)*avgXp*dm)/(1.0*(n+1));

        //Calculate the derivative wiht r xp
        dL1 = 2*(dm*(m*sumX2+n*b*avgXo-sumXY) + n*db*(m*avgXo+b-avgY) + (m*static_cast<long double>(xp)+b-yp)*(dm*static_cast<long double>(xp)+m+db));


        
        //std::cout << "A " << typeid(A).name() << std::endl;
        // std::cout << "B " << typeid(B).name() << std::endl;
        // std::cout << "sumXY " << typeid((sumXY)).name() << std::endl;
        // std::cout << "A " << A << std::endl;
        // std::cout << "B " << typeid(db).name() << std::endl;
        // std::cout << "sumXY" << sumXY << std::endl;
        
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
        //if(true){
        //if(true){
         if(sign<0){
        //    M[0].first = (dL[0].first);
        //    M[0].second = (dL[0].second);
            M.insert({std::make_pair(D[i].first, D[i].second), std::make_pair(dL[i].first, dL[i].second)});
        }
        //Calculate L here itself No need to add to S
        //=========================
        else{
            // std::cout<< "S" << D[i].first << std::endl;
            // std::cout<< "S" << D[i].second << std::endl;

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
    // std::cout << "n: " << n << std::endl;

    //double dL1 =0.0;
    

   // long double sumXo = std::accumulate(current_keyset2.begin(), current_keyset2.end(), 0.0);
    long double avgXo = sumXo/(1.0*(n));

    // std::cout << "sumXo: " << sumXo << std::endl;
    // std::cout << "avgXo: " << avgXo << std::endl;

    //long double sumYo = std::accumulate(rankset.begin(), rankset.end(), 0.0);
   //long double sumXYo = std::inner_product(current_keyset2.begin(), current_keyset2.end(), rankset.begin(), 0.0);

    // std::cout << "sumYo: " << sumYo << std::endl;
    // std::cout << "sumXYo: " << sumXYo << std::endl;

   // long double sumX2 = std::inner_product(current_keyset2.begin(), current_keyset2.end(), current_keyset2.begin(), 0.0);
    // long double sumY2o = std::inner_product(rankset.begin(), rankset.end(), rankset.begin(), 0.0);

    //DOUBLE CHECK THIS L_PREV
    //  std::cout << "sumXo " << sumXo << std::endl;
    //  std::cout << "sumYo " << sumYo << std::endl;
    // std::cout << "sumXYo " << sumXYo << std::endl;
    // std::cout << "sumY2o " << sumY2o << std::endl;
    // std::cout << "n " << n << std::endl;

    if(L_p){
        // long double top = pow(sumXYo/n - (sumXo/n)*(sumYo/n),2);
        // long double bottom = ((sumX2/n) - pow(sumXo/n,2));
        // double slope = (
        // n * sumXYo - sumXo * sumYo) /
        // (1.0*n * sumX2 - sumXo * sumXo);

        // double intercept = 
        // (sumYo - slope * sumXo) / (1.0*n);

        

        // long double L_prev_temp = 0;
        // for(int l =0; l < n; l++){
        //     // L_prev_temp += pow(static_cast<int>(slope*current_keyset2[l]+intercept)-rankset[l],2);
        //     double pred = slope*current_keyset2[l]+intercept;
        //     //std::cout << "pred " << static_cast<int>(pred) - static_cast<int>(rankset[l])<< std::endl;
        //     //L_prev_temp += pow((static_cast<int>(pred) - static_cast<int>(rankset[l])),2);
        //     L_prev_temp += pow((pred)-rankset[l],2);
        //     //std::cout << "pred " << pow((pred)-rankset[l],2)<< std::endl;
        // }
        // L_prev = L_prev_temp;
        L_prev = (-pow(sumXYo/(1.0*n) - (sumXo/(1.0*n))*(sumYo/(1.0*n)),2)/((sumX2/(1.0*n))-pow(sumXo/(1.0*n),2)) + sumY2o/(1.0*n) - pow(sumYo/(1.0*n),2))*n;
        //L_prev = (-pow(sumXYo/n - (sumXo/n)*(sumYo/n),2)/((sumX2/n)-pow(sumXo/n,2)) + sumY2o/n - pow(sumYo/n,2))*n;
        //long double minus = sumY2o/(1.0*n) - pow(sumYo/(1.0*n),2);

        //L_prev = (-top/bottom + minus)*n;
        
        //L_prev = (-pow(sumXYo/n - (sumXo/n)*(sumYo/n),2)/((sumX2/n)-pow(sumXo/n,2)) + sumY2o/n - pow(sumYo/n,2))*n;
        // std::cout << "L_prev " << L_prev << std::endl;
        // std::cout << "slope " << slope << std::endl;
        // std::cout << "intercept " << intercept << std::endl;
    //     std::cout << "m " << slope  << std::endl;
    //    std::cout << "b " << intercept  << std::endl;
    //    std::cout << "top " << top  << std::endl;
    //    std::cout << "bottom " << bottom  << std::endl;
    //    std::cout << "minus " << minus  << std::endl;
    //    std::cout << "top/bottom " << top/bottom - minus  << std::endl;
        // std::cout << L_prev  << std::endl;
        if(L_prev < 0.000000001){
            L_prev = 0;
            //std::cout << "MINUS "  << std::endl;
        }
        if(L_prev < 0){
            std::cout << "MINUS "  << std::endl;
        }
    }
        

    //  std::cout << "sumX2: " << sumX2 << std::endl;
    //  std::cout << "sumY2o: " << sumY2o << std::endl;



    //double sumY;
    

    int yp;
    KEY_TYPE xp;
    long double avgXp;
    double avgYp;
    double sumY;
    double avgY;
    long double sumXY;
    // long double A;
    // long double B;
    long double m;
    long double b;

    int yp_prev;
    // long double dm;
    // long double db;

    //For all the derivatives
    for(int i = 0; i < S.size(); i++) {
            //std::cout << "=================" << std::endl;


        xp = S[i];
 
        yp = compute_rank_for_endpoint(xp, current_keyset2);
        // std::cout << "xp: " << xp << std::endl;
        // std::cout << "yp: " << yp << std::endl;

        avgXp = (sumXo + static_cast<long double>(xp))/(1.0*(n+1));

        avgYp = (sumYo + n)/(1.0*(n+1));

        // std::cout << "Avgxp: " << avgXp << std::endl;
        // std::cout << "Avgyp: " << avgYp << std::endl;

        // std::cout << "sumX2: " << sumX2 << std::endl;
        // std::cout << "sumY2o: " << sumY2o << std::endl;


        sumY = sumYo + n -yp;

        avgY = sumY/(1.0*n);

        // std::cout << "sumY: " << sumY << std::endl;
        // std::cout << "avgY: " << avgY << std::endl;

        //double check this yp (if it is yp or yp -1 )
        //also check current keyset or current keyset2
        
        // //sumXY = sumXYo + std::accumulate(current_keyset2.begin()+yp, current_keyset2.end(), 0.0);

       

        // //std::cout<< n ;
        // sumXY = sumXYo;
        // for(int j = yp; j < n; j++ ){
        //     //std::cout<< "j " ;
        //     //sumXY += static_cast<long double>(current_keyset2[j]);
        //     sumXY += static_cast<long double>(current_keyset2[j]);
        //     //std::cout<< current_keyset2[j] << std::endl;
        // }

        if(i == 0){
            sumXY = sumXYo;
            for(int k = yp; k < n; k++ ){
            //std::cout<< "j " ;
            //sumXY += static_cast<long double>(current_keyset2[j]);
                sumXY += static_cast<long double>(current_keyset2[k]);
            //std::cout<< current_keyset2[j] << std::endl;
            }
            //sumXY_prev = sumXY;
        }
        else{  
            for(int k = yp_prev; k < yp; k++ ){
            //std::cout<< "j " ;
            //sumXY += static_cast<long double>(current_keyset2[j]);
                sumXY -= static_cast<long double>(current_keyset2[k]);
            //std::cout<< current_keyset2[j] << std::endl;
            }

            //sumXY = sumXY_prev;
        }
         //std::cout << "sumXY: " << sumXY << std::endl;

        //std::cout << "sumXY: " << sumXY << std::endl;

        // A = (n+1)*(sumX2+xp*xp) - pow(n*avgXo + xp,2);
        // B = (n+1)*(sumXY+xp*yp) - (pow(n+1,2)*avgXp*avgYp);

        // std::cout << "A: " << A << std::endl;
        // std::cout << "B: " << B << std::endl;

        long double bottom_m = (static_cast<long double>(sumX2)+static_cast<long double>(xp)*static_cast<long double>(xp)-(n+1)*static_cast<long double>(avgXp)*static_cast<long double>(avgXp));
        // m = ((sumXY+static_cast<long double>(xp)*yp)-(n+1)*avgXp*avgYp)/(sumX2+static_cast<long double>(xp)*static_cast<long double>(xp)-(n+1)*avgXp*avgXp);
        m = ((sumXY+static_cast<long double>(xp)*yp)-(n+1)*avgXp*avgYp)/bottom_m;

        b = avgYp - m*avgXp;

        //  std::cout << "m: " << ((sumXY+static_cast<long double>(xp)*yp)-(n+1)*avgXp*avgYp) << std::endl;
        //  std::cout << "m: " << sumX2+static_cast<long double>(xp)*static_cast<long double>(xp) << std::endl;
        //  std::cout << "m: " << (n+1)*avgXp*avgXp << std::endl;
        //  std::cout << "m: " << bottom_m << std::endl;
        //  std::cout << "m: " << m << std::endl;
         // std::cout << "b: " << b << std::endl;

        // dm = (A*(n*(yp-avgY))-B*(2*n*(xp-avgXo)))/(A*A);

        // db = -(m+(n+1)*avgXp*dm)/(n+1);

        // std::cout << "dm: " << dm << std::endl;
        // std::cout << "db: " << db << std::endl;

        //Calculate the loss 
        //L[i] = (m*m*sumX2 +2*n*m*b*avgXo - 2*m*sumXY+n*b*b - 2*n*b*avgY + sumY2o + (n*n) - yp*yp + pow(m*xp+b-yp,2))/(1.0*(n+1));
        // double L_temp = 0;
        // //std::vector<size_t> rankset_new = tag_sort(current_keyset2);

        // for(int l =0; l < n ; l++){
        //     int rank = static_cast<int>(rankset[l]);
        //     if(rank >= yp){
        //         rank++;
        //     }
        //     double pred = m*current_keyset2[l] + b;
        //    //  L_temp += pow(static_cast<int>(pred) - rank,2);
        //     L_temp += pow((pred) - rank,2);
        // }
        // // L_temp += pow(static_cast<int>(m*xp + b) - yp,2);
        // L_temp += pow((m*xp + b) - yp,2);

        // L[i] = L_temp;
         L[i] = (m*m*sumX2 +2*n*m*b*avgXo - 2*m*sumXY+n*b*b - 2*n*b*avgY + sumY2o + (n*n) - yp*yp + pow(m*static_cast<long double>(xp)+b-yp,2));
        
       // std::cout << "L: " << L[i] << std::endl;
        //std::cout << "L: " << L[i] << std::endl;
        // if(L[i] < 0){
        //     std::cout << "L: " << L[i] << std::endl;
        // }
        if(L[i] < 0.000000001){
            L[i] = 0;
            //std::cout << "MINUS "  << std::endl;
        }

        yp_prev = yp;
        }

}

void minimum_points(std::map<std::pair<KEY_TYPE, KEY_TYPE>, std::pair<  double,   double>>& M, std::vector<KEY_TYPE>& S){
        double val;
        KEY_TYPE xp;
       for (const auto& entry : M) {
            //std::cout << "Key: (" << entry.first.first << ", " << entry.first.second << "), Value: (" << entry.second.first << ", " << entry.second.second << ")" << std::endl;

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

    //std::cout << std::endl << "Number of legitimate keys: " << n << std::endl;

    // Number of poisoning keys P
    int P = int(poisoning_threshold * n);
    //int P = 100;
    //std::cout << "Number of poisoning keys to be generated: " << P << std::endl;

    std::set<KEY_TYPE> poisoning_keys;

    //std::set<KEY_TYPE> poisoning_keys_temp;
    //int count_tries = 0;
    //double loss_value = 0;
    //std::vector<size_t> rankset2 = tag_sort(keyset);

    long double sumXo = 0;
    //long double sumXo = std::accumulate(keyset.begin(), keyset.end(), 0.0L) ;
    
    double sumYo = 0;
    //double sumYo = std::accumulate(rankset.begin(), rankset.end(), 0.0) ;
    
    long double sumX2 = 0;
    //long double sumX2 = std::inner_product(keyset.begin(), keyset.end(), keyset.begin(), 0.0L);

   
    double sumY2o = 0;
    //double sumY2o = std::inner_product(rankset.begin(), rankset.end(), rankset.begin(), 0.0);

    long double sumXYo = 0;
    //long double sumXYo = std::inner_product(keyset.begin(), keyset.end(), rankset.begin(), 0.0L);

     for(int i =0; i < keyset.size(); i++){
        sumX2 += static_cast<long double>(keyset[i]) * keyset[i];
        sumXo += static_cast<long double>(keyset[i]);
        sumXYo += static_cast<long double>(keyset[i])*rankset[i];

        sumYo += rankset[i];
        sumY2o += rankset[i]*rankset[i];

    }

    // long double sumXo_prev = 0;
    // long double sumYo_prev = 0;
    // long double sumXYo_prev = 0;
    // long double sumX2_prev = 0;
    // long double sumY2o_prev = 0;
    //std::vector <KEY_TYPE> keyset_prev;

    KEY_TYPE xp_prev;
    size_t rank_prev;

    
    //auto poisoning_start = std::chrono::high_resolution_clock::now();
     //=====================================
     //int c = 0;
    while(poisoning_keys.size() < P){
        
       // c++;
        //===========================================
        
   

        std::vector <KEY_TYPE> keyset2;
        std::merge(keyset.begin(), keyset.end(), poisoning_keys.begin(), poisoning_keys.end(), std::back_inserter(keyset2));
        n = keyset2.size();
        //std::cout << std::endl << "The current keys number: " << n << std::endl;
        //std::cout << "Current status: " << poisoning_keys.size() << " out of " << P << " poisoning keys generated." << std::endl;

       
        std::vector<size_t> rankset = tag_sort(keyset2);

        //std::vector<size_t> rankset = std::iota(std::begin(keyset2), std::end(keyset2), 0);

        //std::cout << "Ranked: "<< std::endl;
       
        
        if(poisoning_keys.size() > 0){
            sumXo = sumXo + static_cast<long double>(xp_prev);
            sumYo = sumYo + (n-1);
            // sumXYo = sumXYo + std::accumulate(keyset2.begin()+(rank_prev+1), keyset2.end(), 0.0) +rank_prev*static_cast<long double>(xp_prev);
            sumXYo += static_cast<long double>(xp_prev)*rank_prev;


            //This maybe a time issue
            
            
            for(int i = (rank_prev+1); i<keyset2.size();i++){
                sumXYo+=keyset2[i];
            }
            

            sumX2 = sumX2 + static_cast<long double>(xp_prev)*static_cast<long double>(xp_prev);
            sumY2o = sumY2o + (n-1)*(n-1) ;
            //std::inner_product(keyset2.begin(), keyset2.end(), rankset.begin(), 0.0)

            // std::cout << "Rank: " <<  rank_prev << std::endl;
            // std::cout << "sumXo: " <<  sumXo << std::endl;
            // std::cout << "sumYo: " <<  sumYo << std::endl;
            // std::cout << "sumX2: " <<  sumX2 << std::endl;
            // std::cout << "sumY2o: " <<  sumY2o << std::endl;
            
            
        }

        //Time 13.5 ms
        
        
        
        
        //===========================================
        // std::cout << "Rank: " <<  rank_prev << std::endl;
            // std::cout << "sumXo: " <<  sumXo << std::endl;
            // std::cout << "sumYo: " <<  sumYo << std::endl;
            // std::cout << "sumX2: " <<  sumX2 << std::endl;
            // std::cout << "sumY2o: " <<  sumY2o << std::endl;
            // std::cout << "sumXYo: " <<  sumXYo << std::endl;
        

        // Partition the non-occupied keys into subsequences such that each subsequence consists of consecutive non-occupied keys;
        // Extract the endpoints of each subsequence and sort them to construct the new sequence of endpoints S(i), where i <= 2(n + j);

        //Calculate the endpoints and which ones to add to caluclate loss or derv

        // To store the calculate loss
        
        //===========================================
        
        

        std::vector<KEY_TYPE> S;
        std::vector<KEY_TYPE> S_temp;

        // To store the calculate derv (stores the two potential endpoints)
        std::vector<std::pair<KEY_TYPE, KEY_TYPE>> D;
        //std::vector<int> M;

        std::map<std::pair<KEY_TYPE, KEY_TYPE>, std::pair<  double,   double>> M;

        //std::cout << "Before part: "<< std::endl;
        std::vector<std::pair<KEY_TYPE, KEY_TYPE>> G;
        
        //12 (66 s, 1200) ok, 25 (83 s 17) is best, 50 (100 s) is same as normal
       

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
        
        //===========================================

        //std::cout << "Partioned: "<< std::endl;

        //std::cout << S.size() << std::endl;

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
        
        //===========================================
       
        minimum_points(M, S);

        

        //std::cout << "Minimum: "<< std::endl;

        //Now calculate the L for all in S
        std::sort(S.begin(), S.end());

        //===========================================


        //std::cout << "Length of endpoints: " << S.size() << std::endl;

        // std::cout << "Current keys: " << std::endl;
        // for (uint64_t i: keyset2)
        //     std::cout << i << ' ';

        // std::cout << "Current ranks: " << std::endl;
        // for (size_t i: rankset)
        //     std::cout << i << ' ';


        // Endpoint keys are continuously increasing
        
        //std::cout << "Endpoint keys: " << std::endl;
        // for (uint64_t i: S)
        //     std::cout << i << ' ';
        


        //Calculate the L_prev for the first one
        // auto loss_start = std::chrono::high_resolution_clock::now();
        //===========================================

        if(poisoning_keys.size()==0){
            cal_loss(keyset2,rankset,
        L, 
        S,true, L_prev,
        sumXo,sumYo,sumXYo,sumX2,sumY2o);

        original_L_prev = L_prev;
        //std::cout << "L_Prev "<< L_prev << std::endl;

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

        //std::cout << "Loss: "<< std::endl;
        //  cal_loss(keyset2,rankset,
        // L, 
        // S);

       

        int optimal_key_index;
        //auto minElement = L.lower_bound(0);

       
        //===========================================
        
         auto maxElementIterator = std::min_element(L.begin(), L.end(),[](const std::pair<size_t,  long double>& a, const std::pair<size_t,  long double>& b) {
            return a.second < b.second; // Compare by the second element (the double value)
        });

        //std::cout << "Found: "<< std::endl;

        optimal_key_index = maxElementIterator->first;
        //std::cout << "Found: "<< S[optimal_key_index] << " " << maxElementIterator->second << std::endl;
       //if(S.size()==0 || L_prev == 0){    
        if(maxElementIterator->second >= L_prev || S.size()==0 || L_prev == 0){
       // if(poisoning_keys.size()>0 && maxElementIterator->second >= L_prev || S.size()==0){
        //    if(poisoning_keys.size() > 100  || S.size() == 0){
        //if(poisoning_keys.size()>0 && maxElementIterator->second > L_prev || S.size()==0){
            //std::cout << "Stopping " << std::endl;
            //std::cout << "L "<< L_prev << std::endl;
            // if(count_tries > 0){
            //     poisoning_keys = poisoning_keys_temp;
            // }
            break;
        }

        // if(maxElementIterator->second >= L_prev ){
        //     count_tries++;
        //     if(count_tries == 1){
        //         poisoning_keys_temp = poisoning_keys;
        //         loss_value = L_prev;
        //     }
        //     else if(count_tries == 10){
        //         poisoning_keys = poisoning_keys_temp;
        //         break;
        //     }
        // }
        // else if(count_tries > 0 && maxElementIterator->second < loss_value){
        //     count_tries = 0;
        // }
        // else if(count_tries > 0 && maxElementIterator->second >= loss_value){
        //     count_tries++;
        //     if(count_tries == 10){
        //         poisoning_keys = poisoning_keys_temp;
        //         break;
        //     }
        // }
        //std::cout << "count "<< c << std::endl;
        
        poisoning_keys.insert(S[optimal_key_index]);
         L_prev = maxElementIterator->second;
         //std::cout << "L "<< L_prev << std::endl;

        //===========================================
         

        xp_prev = S[optimal_key_index];
        rank_prev = compute_rank_for_endpoint(xp_prev, keyset2);

        
        
        if(poisoning_keys.size()== P){
            // if(count_tries > 0){
            //     poisoning_keys = poisoning_keys_temp;
            // }
            break;
            //std::cout << "L "<< L_prev << std::endl;
        }
        // else if(poisoning_keys.size()== P + 10){
        //     if(count_tries > 0){
        //         poisoning_keys = poisoning_keys_temp;
        //     }
        //     break;
        // }

        
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
// template <class KEY_TYPE, class PAYLOAD_TYPE>
std::vector<KEY_TYPE> perform_poisoning(std::vector<KEY_TYPE> & legit_keys, double poisoning_threshold,bool max =true) {


   //std::cout << "Legitimate keys: " << std::endl;
    // for (double i: legit_keys)
    //     std::cout << i << ' ';


    std::vector<size_t> legit_ranks = tag_sort(legit_keys);


    //std::cout << std::endl << "Legitimate ranks: " << std::endl;

    // for (size_t i: legit_ranks)
    //     std::cout << i << ' ';

    


    std::set<KEY_TYPE> poisoning_keys = obtain_poisoning_keys(poisoning_threshold, legit_keys, legit_ranks,max);

    /*
     * Merge poisoning keys with legitimate keys to generated poisoned keyset
     */
    
        
    std::vector <KEY_TYPE> poisoned_keyset;
    std::merge(legit_keys.begin(), legit_keys.end(), poisoning_keys.begin(), poisoning_keys.end(), std::back_inserter(poisoned_keyset));
    std::sort(poisoned_keyset.begin(), poisoned_keyset.end());

    //std::cout << poisoning_keys.size() << std::endl;


   //std::cout << "Poisoned keyset: " << poisoning_keys.size() << std::endl;
    // for (double i: poisoned_keyset)
    //     std::cout << i << ' ';
    // std::cout << std::endl;

    return poisoned_keyset;

}

/* TO DO 

    find all points that are not the keys that the poisoning can take place
    Do the L thing for all and the same as before
    Add PGM and check performance
    Add time measuring
*/