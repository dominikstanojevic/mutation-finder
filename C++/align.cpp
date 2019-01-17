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

const int MATCH =  5;
const int DIFF  = -4;
const int EMPTY = -8;

int gotovo = 0;

int step;

vector<unordered_map<int, unordered_map<short, int>>> AlignAll(string &reference, vector<string> &queries, int w, int k, int eps){
    vector<unordered_map<int, unordered_map<short, int>>> mapInfo {3};
    step = k * w;

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
            }
        }
    }

    return mapInfo;
}

inline char r(char x) {
    switch (x) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            cout << "Jebem ti život." << endl;
            exit(-1);
    }
}

vector<unordered_map<int, info>> AlignRegion(string &reference, string &query, mapping_result region){
    if (region.end < region.start) {
        throw out_of_range("Ref region out of range.");
    }
    if (region.maxQ < region.minQ) {
        throw out_of_range("Query region out of range.");
    }

    int startPos = (region.start - step) <= 0 ? 0 : (region.start - step);
    int end = (region.end + step) > reference.size() ? reference.size() : (region.end + step);

    int startQ = (region.minQ - step) <= 0 ? 0 : (region.minQ - step);
    int endQ = (region.maxQ + step) > query.size() ? query.size() : (region.maxQ + step);

    if(region.reversed) {
        string q = "";
        for(int i = query.size() - 1; i >= 0; i--) {
            q += r(query[i]);
        }

        return align(reference, q, startPos, end, startQ, endQ).infoMap;
    } else {
        return align(reference, query, startPos, end, startQ, endQ).infoMap;
    }
}

inline void emit_change(unordered_map<int, info> &map, int &pos, string &t, int startQ, int &i, int &j) {
    map[pos] = info(pos, I_CHA, t[startQ + j - 1]);
    i--; j--; pos--; 
}
inline void emit_ins(unordered_map<int, info> &map, int pos, string &t, int startQ, int &j) {
    map[pos+1] = info(pos + 1, I_INS, t[startQ + j - 1]);
    j--;
}
inline void emit_del(unordered_map<int, info> &map, int &pos, string &t, int &i) {
    map[pos] = info(pos, I_DEL, '0');
    i--; pos--;
}

inline bool check_diag(int *arr, string &s, string &t, int i, int j, int cols, int startPos, int startQ) {
    if ((i > 0) && (j > 0) &&
        (arr[i*cols + j] == (arr[(i-1)*cols + j-1] + ((s[startPos + i-1] == t[startQ + j-1]) ? MATCH : DIFF)))){
            return true;
        } else return false;
}

inline bool check_up(int *arr, int i, int j, int cols) {
    if ((i > 0) && (arr[i*cols + j] == (arr[(i-1)*cols + j] + EMPTY))) return true;
    return false;
}

inline bool check_left(int *arr, int i, int j, int cols) {
    if ((j > 0) && (arr[i*cols + j] == (arr[i*cols + j - 1] + EMPTY))) return true;
    return false;
}

alignmnent_info align(string &s, string &t, int startPos, int end, int startQ, int endQ){
    int rows = end - startPos + 1;
    int cols = endQ - startQ + 1;
    int elems = rows * cols;

    int *arr = (int *) malloc(sizeof(int) * elems);
    
    for (int j = 0; j < cols; ++j){
        arr[j] = j * EMPTY;
    }

    for (int i = 1; i < rows; ++i){
        arr[i * cols] = 0;
        for (int j = 1; j < cols; ++j){
            int match = (s[startPos + i-1] == t[startQ + j-1]) ? MATCH : DIFF;
            match += arr[(i-1) * cols + j-1];
            
            int del = arr[(i-1) * cols + j] + EMPTY;
            int insert = arr[i * cols + j-1] + EMPTY;

            arr[i * cols + j] = max(max(match, del), insert);
        }
    }

    int i = 0; 
    int j = cols - 1;
    int maxValue = arr[j];
    for (int e = 1; e < rows; e++) {
        if (arr[e * cols + j] >= maxValue) {
            i = e;
            maxValue = arr[e * cols + j];
        }
    }

    alignmnent_info al_info {maxValue};
    auto& infoMap = al_info.infoMap;
    
    int pos = startPos + i - 1;
    while (j > 0){ 
        bool diag = check_diag(arr, s, t, i, j, cols, startPos, startQ);
        bool up = check_up(arr, i, j, cols);
        bool left = check_left(arr, i, j, cols);

        if (diag && up && left) {
            switch (rand() % 3) {
                case 0:
                    emit_change(infoMap[0], pos, t, startQ, i, j);
                    break;
                case 1:
                    emit_ins(infoMap[1], pos, t, startQ, j);
                    break;
                case 2:
                    emit_del(infoMap[2], pos, t, i);
                    break;
            }
        } else if (diag && up) {
            switch (rand() % 2) {
                case 0:
                    emit_change(infoMap[0], pos, t, startQ, i, j);
                    break;
                case 1:
                    emit_del(infoMap[2], pos, t, i);
                    break;
            }
        } else if (diag && left) {
            switch (rand() % 2) {
                case 0:
                    emit_change(infoMap[0], pos, t, startQ, i, j);
                    break;
                case 1:
                    emit_ins(infoMap[1], pos, t, startQ, j);
                    break;
            }
        } else if (up && left) {
            switch (rand() % 2) {
                case 0:
                    emit_ins(infoMap[1], pos, t, startQ, j);
                    break;
                case 1:
                    emit_del(infoMap[2], pos, t, i);
                    break;
            }
        } else if (diag) emit_change(infoMap[0], pos, t, startQ, i, j);
        else if (left)  emit_ins(infoMap[1], pos, t, startQ, j);
        else emit_del(infoMap[2], pos, t, i);

        /* if ((i > 0) && (j > 0) &&
                (arr[i*cols + j] == (arr[(i-1)*cols + j-1] + ((s[startPos + i-1] == t[startQ + j-1]) ? MATCH : DIFF)))){
            infoMap[0][pos] = info(pos, I_CHA, t[startQ + j - 1]);

            i--; j--; pos--; 
        } else if ((i > 0) && (arr[i*cols + j] == (arr[(i-1)*cols + j] + EMPTY))) {
            if ((j > 0) && (arr[i*cols + j] == (arr[i*cols + j - 1] + EMPTY))) {
                if (rand() % 2 == 0) {
                    infoMap[2][pos] = info(pos, I_DEL, '0');
                    i--; pos--;
                } else {
                    infoMap[1][pos+1] = info(pos + 1, I_INS, t[startQ + j - 1]);
                    j--;
                }
            } else {
                infoMap[2][pos] = info(pos, I_DEL, '0');
                i--; pos--;
            }
        } else {
            infoMap[1][pos+1] = info(pos + 1, I_INS, t[startQ + j - 1]);
            j--;
        } */
    }

    free(arr);
    return al_info;
}