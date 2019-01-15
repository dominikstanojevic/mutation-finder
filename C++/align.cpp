#include "align.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>

#include "map.h"

using namespace std;

const int MATCH =  4;
const int DIFF  = -1;
const int EMPTY = -2;

int gotovo = 0;

void AlignAll(string &reference, vector<string> &queries, int w, int k, int eps,
        vector< map<int, vector<info> > > &mapInfo){
    mapInfo.resize(3);

    auto table = Index(&reference, 1, w, k);
    for (auto query : queries){
        auto regions = Map(table, query, w, k, eps);

        vector< map<int, info> > result; 
        AlignQuery(reference, query, regions, result);

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
}

void AlignQuery(string &reference, string &query, vector<mapping_result> &regions,
        vector< map<int, info> > &infoMap){
    infoMap.resize(3);
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
        return;
    }

    int i = 0;
    while (best.alY[i] == '-') ++i;
    int j = best.alY.size() - 1;
    while (best.alY[j] == '-') --j;

    int pos = bestStart + i;
    for (; i < j; ++i){
        if (best.alY[i] == '-'){
            infoMap[2][pos] = info(pos, I_DEL, '0');
            pos++;
        } else if (best.alX[i] == '-') {
            infoMap[1][pos] = info(pos, I_INS, best.alY[i]);
        } else {
            infoMap[0][pos] = info(pos, I_CHA, best.alY[i]);
            pos++;
        }
    }
}

alignmnent_info align(string &s, string &t){
    int rows = s.size() + 1;
    int cols = t.size() + 1;

    int *array = (int*)malloc(sizeof(int) * rows * cols);
    
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

    int i = 0;
    int j = cols - 1;
    int maxValue = array[j];

    for (int e = 1; e < rows; ++e){
        if (array[e*cols + j] > maxValue){
            i = e;
            maxValue = array[e*cols + j];
        }
    }

    string aS = "", aT = "";
    aS = s.substr(i, rows-1-i);
    for (int e = 0, end = rows - 1 - i; e < end; ++e){
        aT.append("-");
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