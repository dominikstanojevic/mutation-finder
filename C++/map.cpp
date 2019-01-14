#include "map.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include "sketch.h"

using namespace std;

map<long, set<index_data> > Index(string sequences[], int n_seq, int w, int k){
    map<long, set<index_data> > table;

    for (int i = 0; i < n_seq; ++i){
        vector<minimizer> minimizers = MinimizerSketch(sequences[i], w, k);

        for (auto m : minimizers){
            if (table.count(m.hash) == 0){
                table[m.hash] = set<index_data>();
            }
            set<index_data> *entry = &table[m.hash];
            entry->insert(index_data(i, m.end_position, m.strand));
        }
    }

    return table;
}

void Map(map<long, set<index_data> > &table, string &q, int w, int k, int eps){
    vector<map_data> arr;
    vector<minimizer> minimizers = MinimizerSketch(q, w, k);

    for (auto m : minimizers){
        if (table.count(m.hash) == 0){
            continue;
        }
        for (auto ind : table[m.hash]){
            if (m.strand == ind.r){
                arr.push_back(map_data(ind.t, 0, m.end_position - ind.i, ind.i));
            } else {
                arr.push_back(map_data(ind.t, 1, m.end_position + ind.i, ind.i));
            }
        }
    }

    sort(arr.begin(), arr.end());

    for (auto x : arr){
        cout << x.t << " " << x.r << " " << x.c << " " << x.i << endl;
    }
}
