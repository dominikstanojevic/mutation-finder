#ifndef SKETCH_H
#define SKETCH_H

#include <iostream>
#include <set>

using namespace std;

struct minimizer {
    uint64_t hash;
    int end_position;
    int strand;

    friend bool operator < (const minimizer &x, const minimizer &y){
        if (x.hash != y.hash) return x.hash < y.hash;
        if (x.end_position != y.end_position) return x.end_position < y.end_position;
        return x.strand < y.strand;
    }
};

set<minimizer> MinimizerSketch(string s, int w, int k);

#endif // SKETCH_H