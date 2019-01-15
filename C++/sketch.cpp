#include "sketch.h"

#include <assert.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <set>
#include <vector>

#include "utils.h"

using namespace std;

int CharToBase(char c){
    if (c == 'A' || c == 'a') return 0;
    if (c == 'C' || c == 'c') return 1;
    if (c == 'G' || c == 'g') return 2;
    if (c == 'T' || c == 't') return 3;

    cerr << "Unsupported allele: " << c << endl;
    exit(1);
}

int ComplementBase(int base){
    assert(base >= 0 && base < 4);
    return 3 - base;
}

set<minimizer> MinimizerSketch(string &s, int w, int k){
    set<minimizer> minimizers;
    int len = s.size();

    uint64_t mask = (1ULL << (2*k)) - 1;
    uint64_t kmer[2] = {0, 0};
    minimizer min;
    min.hash = UINT64_MAX;
    min.end_position = 0;

    minimizer *window; 
    window = (minimizer *) malloc(w * sizeof(minimizer));
    memset(window, 0, w * sizeof(minimizer));
    int window_index = 0;

    int min_index = 0;
    minimizer tmp_min;
    for (int i = 0; i < len; ++i){
        char base = CharToBase(s[i]);
        kmer[0] = (kmer[0] << 2 | base) & mask;
        kmer[1] = (kmer[1] >> 2) | (ComplementBase(base) << (2 * (k - 1)));

        if (kmer[0] == kmer[1]) continue;
        int strand = kmer[0] < kmer[1] ? 0 : 1;

        tmp_min.hash = Hash(kmer[strand], mask);
        tmp_min.end_position = i - k + 1;
        tmp_min.strand = strand;
        window[window_index] = tmp_min;

        if (tmp_min.hash <= min.hash){
            if (min.end_position >= 1){
                minimizers.insert(min);
            }
            min = tmp_min;
            min_index = window_index;
        } else if (window_index == min_index){
            if (min.end_position >= 1){
                minimizers.insert(min);
            }

            min.hash = UINT64_MAX;
            for (int j = window_index + 1; j < window_index + w + 1; ++j){
                if (window[j % w].hash < min.hash){
                    min = window[j % w];
                    min_index = j % w;
                }
            }

            if (min.end_position >= 1){
                for (int j = window_index + 1; j < window_index + w + 1; ++j){
                    if (window[j % w].hash == min.hash && window[j % w].end_position != min.end_position){
                        minimizers.insert(window[j % w]);
                    }
                }
            }
        }

        window_index++;
        if (window_index >= w){
            window_index = 0;
        }
    }

    if (min.hash < UINT64_MAX){
        minimizers.insert(min);
    }

    free(window);
    return minimizers;
}