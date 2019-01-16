#include "align.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include <climits>
#include <unordered_map>

#include "map.h"

using namespace std;

const int MATCH =  4;
const int DIFF  = -1;
const int EMPTY = -2;

int gotovo = 0;

vector<int> arr;
int maxArr = 0;

vector<unordered_map<int, vector<info>>> AlignAll(string &reference, vector<string> &queries, int w, int k, int eps){
    vector<unordered_map<int, vector<info>>> mapInfo {3};

    auto table = Index(&reference, 1, w, k);
    for (auto query : queries){
        auto regions = Map(table, query, w, k, eps);

        auto result = AlignQuery(reference, query, regions);

        for (int i = 0; i < 3; ++i){
            for (auto entry : result[i]){
                if (mapInfo[i].count(entry.first) == 0){
                    mapInfo[i][entry.first] = vector<info>();
                }
                mapInfo[i][entry.first].push_back(entry.second);
            }
        }

        cout << gotovo++ << endl;
    }

    return mapInfo;
}

vector<unordered_map<int, info>> AlignQuery(string &reference, string &query, vector<mapping_result> &regions){
    alignmnent_info best {INT_MIN};
    int bestStart = -1;

    for (auto region : regions) {
        //cout << region.start << " " << region.end << " ";

        int len = region.end - region.start;
        int diff = (query.size() - len) / 2;
        diff = 0;

        int start = diff > 0 ? (region.start - diff) : region.start;
        if (start < 0) start = 0;

        int end = diff > 0 ? (region.end + diff) : region.end;
        if (end > reference.size()) end = reference.size();

        alignmnent_info result = align(reference, query, start, end);
        if (best.value < result.value){
            best = result;
            bestStart = start;
        }

        //cout << gotovo << endl;
    }
    //cout << endl;

    cout << best.infoMap[0].size() << " " << best.infoMap[1].size() << " " << best.infoMap[2].size() << endl;

    return best.infoMap;
}

alignmnent_info align(string &s, string &t, int startPos, int end){
    int rows = end - startPos + 1;
    int cols = t.size() + 1;
    int elems = rows * cols;

    if (elems > maxArr) {
        arr.resize(elems);
        maxArr = elems;
    }
    
    for (int j = 1; j < cols; ++j){
        arr[j] = j * EMPTY;
    }

    for (int i = 1; i < rows; ++i){
        arr[i * cols] = 0;
        for (int j = 1; j < cols; ++j){
            int match = (s[startPos + i-1] == t[j-1]) ? MATCH : DIFF;
            match += arr[(i-1) * cols + j-1];
            
            int del = arr[(i-1) * cols + j] + EMPTY;
            int insert = arr[i * cols + j-1] + EMPTY;

            arr[i * cols + j] = max(max(match, del), insert);
        }
    }

    int i = 0;
    int j = cols - 1;
    int maxValue = arr[j];

    for (int e = 1; e < rows; ++e){
        if (arr[e*cols + j] > maxValue){
            i = e;
            maxValue = arr[e*cols + j];
        }
    }

    alignmnent_info al_info {maxValue};
    auto& infoMap = al_info.infoMap;

    unordered_map<int, info> temp;
    bool last = false;

    int pos = startPos + i - 1;
    bool first = true;
    while (i > 0 || j > 0){
        if ((i > 0) && (j > 0) &&
                (arr[i*cols + j] == (arr[(i-1)*cols + j-1] + ((s[startPos + i-1] == t[j-1]) ? MATCH : DIFF)))){
            infoMap[0][pos] = info(pos, I_CHA, t[j - 1]);
            if (last) {
                infoMap[2].insert(temp.begin(), temp.end());
                temp.clear();
            }

            i--; j--; pos--; 
            first = last = false;
        } else if ((j > 0) && (arr[i*cols + j] == (arr[i*cols + j-1] + EMPTY))) {
            infoMap[1][pos] = info(pos, I_INS, t[j - 1]);
            if (last) {
                infoMap[2].insert(temp.begin(), temp.end());
                temp.clear();
            }

            j--;
            first = last = false;
        } else {
            if (first) {
                i--;
                pos--;
                continue;
            }

            temp[pos] = info(pos, I_DEL, '0');

            i--;
            pos--;
            last = true;
        }
    }

    return al_info;
}