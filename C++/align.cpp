#include "align.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include <climits>
#include <unordered_map>
#include <stdexcept>

#include "map.h"

using namespace std;

const int MATCH =  1;
const int DIFF  = -1;
const int EMPTY = -1;

int gotovo = 0;

vector<unordered_map<int, unordered_map<short, int>>> AlignAll(string &reference, vector<string> &queries, int w, int k, int eps){
    vector<unordered_map<int, unordered_map<short, int>>> mapInfo {3};

    int n_queries = queries.size();
    vector<unordered_map<int, info>> results[n_queries];
    auto table = Index(&reference, 1, w, k);

    #pragma omp parallel for num_threads(7)
    for (int i = 0; i < n_queries; i++){
        auto &query = queries[i];
        auto regions = Map(table, query, w, k, eps);

        for (auto region : regions) {
            auto result = AlignRegion(reference, query, region);

            #pragma omp critical
            {
                for (int i = 0; i < 3; i++) {
                    auto &map_option = mapInfo[i];
                    for (auto entry : result[i]) {
                        if (map_option.count(entry.first) == 0) {
                            map_option[entry.first] = unordered_map<short, int>();
                        }

                        if(map_option[entry.first].count(entry.second.base) == 0) {
                            map_option[entry.first][entry.second.base] = 1;
                        } else {
                            map_option[entry.first][entry.second.base] += 1;
                        }
                    }
                }

                //cout << gotovo++ << endl;
            }
        }
    }

    return mapInfo;
}

vector<unordered_map<int, info>> AlignRegion(string &reference, string &query, mapping_result region){
    if (region.end < region.start) {
        throw out_of_range("Ref region out of range.");
    }
    if (region.maxQ < region.minQ) {
        throw out_of_range("Query region out of range.");
    }

    alignmnent_info result = align(reference, query, region.start, region.end, region.minQ, region.maxQ);
    return result.infoMap;
}

alignmnent_info align(string &s, string &t, int startPos, int end, int startQ, int endQ){
    int rows = end - startPos + 1;
    int cols = endQ - startQ + 1;
    int elems = rows * cols;

    int *arr = (int *) malloc(sizeof(int) * elems);
    
    for (int j = 1; j < cols; ++j){
        arr[j] = j * EMPTY;
    }
    for (int i = 0; i < rows; ++i) {
        arr[i * cols] = i * EMPTY;
    }

    for (int i = 1; i < rows; ++i){
        for (int j = 1; j < cols; ++j){
            int match = (s[startPos + i-1] == t[startQ + j-1]) ? MATCH : DIFF;
            match += arr[(i-1) * cols + j-1];
            
            int del = arr[(i-1) * cols + j] + EMPTY;
            int insert = arr[i * cols + j-1] + EMPTY;

            arr[i * cols + j] = max(max(match, del), insert);
        }
    }

    int i = rows - 1; 
    int j = cols - 1;
    int maxValue = arr[i * cols + j];

    alignmnent_info al_info {maxValue};
    auto& infoMap = al_info.infoMap;

    int pos = startPos + i - 1;
    while (i > 0 || j > 0){
        if (arr[i*cols + j] == 0) {
            break;
        }

        if ((i > 0) && (j > 0) &&
                (arr[i*cols + j] == (arr[(i-1)*cols + j-1] + ((s[startPos + i-1] == t[startQ + j-1]) ? MATCH : DIFF)))){
            infoMap[0][pos] = info(pos, I_CHA, t[startQ + j - 1]);

            i--; j--; pos--; 
        } else if ((i > 0) && (arr[i*cols + j] == (arr[(i-1)*cols + j] + EMPTY))) {
            infoMap[2][pos] = info(pos, I_DEL, '0');
            i--; pos--;
        } else {
            infoMap[1][pos+1] = info(pos + 1, I_INS, t[startQ + j - 1]);
            j--;
        }
    }

    free(arr);
    return al_info;
}