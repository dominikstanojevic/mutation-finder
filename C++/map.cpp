#include "map.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "sketch.h"

using namespace std;

const int MIN_COUNT = 4;
const int MIN_READS = 1000;
const int GAP = 10000;

struct region_data{
    vector<map_data> lis;
    int start;
    int end;

    region_data(vector<map_data> _lis, int _start, int _end): lis(_lis), start(_start), end(_end) {

    }
};

map<long, set<index_data> > Index(string sequences[], int n_seq, int w, int k){
    map<long, set<index_data> > table;

    for (int i = 0; i < n_seq; ++i){
        auto minimizers = MinimizerSketch(sequences[i], w, k);

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

inline int GetQPosition(const map_data &x){
    if (x.r == 0){
        return x.i - x.c;
    } else {
        return x.c - x.i;
    }
}

bool CompareLt(int a, int b){
    return a < b;
}
bool CompareGt(int a, int b){
    return a > b;
}



vector<map_data> LongestIncreasingSubsequence(vector<map_data> &arr, int start, int len, bool (*Compare)(int, int) ){
    vector<int> P;
    vector<int> M;
    P.resize(len);
    M.resize(len + 1);
    
    int l = 0;
    for (int i = 0; i < len; ++i){
        int lo = 1, hi = l;

        while (lo <= hi){
            int mid = (int)ceil((lo + hi) / 2.0);
            if (Compare(arr[start + M[mid]].i, arr[start + i].i)) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }

        int newL = lo;
        P[i] = M[newL - 1];
        M[newL] = i;

        if (newL > l){
            l = newL;
        }
    }

    vector<map_data> m;
    m.resize(l);
    int k = M[l];
    for (int i = l - 1; i >= 0; --i){
        m[i] = arr[start + k];
        k = P[k];
    }
    return m;
}

vector<mapping_result> Map(map<long, set<index_data> > &table, string &q, int w, int k, int eps){
    vector<map_data> arr;
    auto minimizers = MinimizerSketch(q, w, k);

    for (auto m : minimizers){
        if (table.count(m.hash) == 0){
            continue;
        }
        for (auto ind : table[m.hash]){
            if (m.strand == ind.r){
                arr.push_back(map_data(ind.t, 0, ind.i - m.end_position, ind.i));
            } else {
                arr.push_back(map_data(ind.t, 1, m.end_position + ind.i, ind.i));
            }
        }
    }
    
    sort(arr.begin(), arr.end());

    vector<mapping_result> maps;

    int b = 0, l = arr.size();
    for (int e = 0; e < l; ++e){
        vector<region_data> regions;
        
        if (e == l - 1 ||
                arr[e + 1].t != arr[e].t ||
                arr[e + 1].r != arr[e].r ||
                arr[e + 1].c - arr[e].c > eps) {
            
            //sublist = arr.subList(b, e+1);
            struct CustomCompare {
                bool operator ()(const map_data &a, const map_data &b)
                {
                    int q1 = GetQPosition(a);
                    int q2 = GetQPosition(b);

                    if (q1 != q2) return q1 < q2;
                    return a.i < b.i;
                }
            };
            sort(arr.begin() + b, arr.begin() + e + 1, CustomCompare());
            vector<map_data> C;
            
            if (arr[e].r == 0) {
                C = LongestIncreasingSubsequence(arr, b, e + 1 - b, CompareLt);
            } else {
                C = LongestIncreasingSubsequence(arr, b, e + 1 - b, CompareGt);                
            }

            int start = 0;
            for (int end = 1; end < C.size(); ++end){
                if (abs(C[end].i - C[end - 1].i) >= GAP){
                    regions.push_back(region_data(C, start, end));
                    start = end;
                }
            }
            regions.push_back(region_data(C, start, C.size()));
            b = e + 1;
        }

        for (auto region : regions){
            int region_len = region.end - region.start;
            if (region_len < MIN_COUNT){
                continue;
            }
            int min, max;

            int minQ = GetQPosition(region.lis[region.start]);
            int maxQ = GetQPosition(region.lis[region.end - 1]);

            if (arr[e].r == 0){
                min = region.lis[region.start].i;
                max = region.lis[region.end - 1].i;
            } else {
                min = region.lis[region.end - 1].i;
                max = region.lis[region.start].i;
            }
            if (max - min >= MIN_READS) {
                maps.push_back(mapping_result(min, max + k, minQ, maxQ + k));
            }
        }
    }
    return maps;
}