#include "align.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "map.h"

using namespace std;

const int MATCH =  4;
const int DIFF  = -1;
const int EMPTY = -2;

map<int, vector<info> > AlignAll(string &reference, vector<string> &queries, int w, int k, int eps){
    map<int, vector<info> > mapInfo;

    auto table = Index(&reference, 1, w, k);
    for (auto query : queries){
        auto regions = Map(table, query, w, k, eps);
        auto result = AlignQuery(reference, query, regions);

        for (auto entry : result){
            if (mapInfo.count(entry.first) == 0){
                mapInfo[entry.first] = vector<info>();
            }
            mapInfo[entry.first].push_back(entry.second);
        }
    }

    return mapInfo;
}

map<int, info> AlignQuery(string &reference, string &query, vector<mapping_result> regions){
    map<int, info> infoMap;
    alignmnent_info best;
    int bestStart = -1;

    for (auto region : regions) {
        int len = region.end - region.start;
        int diff = query.size() - len;

        int start = diff > 0 ? (region.start - diff) : region.start;
        if (start < 0) start = 0;

        int end = diff > 0 ? (region.end + diff) : region.end;
        if (end > reference.size()) end = reference.size();

        string str = reference.substr(start, end-start);
        alignmnent_info result = align(str, query);
        if (best.alX == "" || best.value < result.value){
            best = result;
            bestStart = start;
        }
    }

    if (best.alX == "") {
        return infoMap;
    }

    int i = 0;
    while (best.alX[i] == '-') ++i;
    int j = query.size() - 1;
    while (best.alY[j] == '-') --j;

    for (; i < j; ++i){
        int pos = bestStart + i;

        if (query[i] == '-'){
            infoMap[pos] = info(pos, I_DEL, '0');
        } else if (reference[pos] == '-') {
            infoMap[pos] = info(pos, I_INS, query[i]);
        } else {
            infoMap[pos] = info(pos, I_CHA, query[i]);
        }
    }

    return infoMap;
}

alignmnent_info align(string &s, string &t){
    int rows = s.size() + 1;
    int cols = t.size() + 1;

    int *array = (int*)malloc(sizeof(int) * rows * cols);
    
    for (int i = 0; i < rows; ++i){
        array[i * cols] = 0;
    }
    for (int j = 0; j < cols; ++j){
        array[j] = 0;
    }
    for (int i = 1; i < rows; ++i){
        for (int j = 1; j < cols; ++j){
            int match = (s[i-1] == t[j-1]) ? MATCH : DIFF;
            match += array[(i-1) * cols + j-1];
            int del = array[(i-1) * cols + j] + EMPTY;
            int insert = array[i * cols + j-1] + EMPTY;

            array[i * cols + j] = max(max(match, del), insert);
        }
    }

    int i = rows-1;
    int j = 0;
    int maxValue = array[i * cols + j];

    for (int e = 1; e < cols; ++e){
        if (array[i*cols + e] > maxValue){
            j = e;
            maxValue = array[i*cols + e];
        }
    }

    for (int e = 0; e < rows; ++e){
        if (array[e*cols + cols-1] > maxValue){
            i = e;
            maxValue = array[e*cols + cols-1];
        }
    }

    string aS = "", aT = "";
    if (i == rows-1 && j != cols-1) {
        for (int e = 0; e < cols-1-j; ++e){
            aS.append("-");
        }
        aT = t.substr(j, cols-j-1);
    } else if (i != rows-1 && j == cols-1) {
        aS = s.substr(i, rows-1-i);
        for (int e = 0; e < rows-1-i; ++e){
            aT.append("-");
        }
    }


    string raS = "", raT = "";

    while (i > 0 || j > 0){
        if ((i > 0) && 
                (j > 0) &&
                (array[i*cols + j] == (array[(i-1)*cols + j-1] + ((s[i-1] == t[j-1]) ? MATCH : DIFF)))){
            raS += s[i-1];
            raT += t[j-1];

            i--; j--;
        } else if ((j > 0) && (array[i*cols + j] == (array[i*cols + j-1] + EMPTY))) {
            raS.append("-");
            raT += t[j-1];

            j--;
        } else {
            raS += s[i-1];
            raT.append("-");

            i--;
        }
    }

    reverse(raS.begin(), raS.end());
    reverse(raT.begin(), raT.end());

    aS = raS + aS;
    aT = raT + aT;

    free(array);
    return alignmnent_info(maxValue, aS, aT);
}