#ifndef MAP_H
#define MAP_H

#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

struct index_data{
    int t;
    int i;
    int r;

    index_data(int _t, int _i, int _r): t(_t), i(_i), r(_r) {

    }

    friend bool operator < (const index_data &a, const index_data &b){
        if (a.t != b.t) return a.t < b.t;
        if (a.i != b.i) return a.i < b.i;
        return a.r < b.r;
    }
};

struct map_data{
    int t;
    int r;
    int c;
    int i;

    map_data(int _t, int _r, int _c, int _i): t(_t), r(_r), c(_c), i(_i) {

    }

    friend bool operator < (const map_data &a, const map_data &b){
        if (a.r != b.r) return a.r < b.r;
        return a.c < b.c;
    }
};

struct mapping_result{
    int start;
    int end;

    mapping_result(int _start, int _end): start(_start), end(_end) {

    }

    friend bool operator < (const mapping_result &a, const mapping_result &b){
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    }
};

map<long, set<index_data> > Index(string sequences[], int n_seq, int w, int k);
vector<mapping_result> Map(map<long, set<index_data> > &table, string &q, int w, int k, int eps);

#endif // MAP_H