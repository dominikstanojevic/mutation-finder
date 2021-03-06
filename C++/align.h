#ifndef ALIGN_H
#define ALIGN_H

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <unordered_map>

#include "map.h"

using namespace std;

#define I_NULL (short)-1
#define I_DEL (short)1
#define I_INS (short)2
#define I_CHA (short)3

#define I_A (short)11
#define I_C (short)12
#define I_G (short)13
#define I_T (short)14

struct info {
    int pos;
    short option;
    short base;

    info(int _pos, short _option, short _base):
        pos(_pos),
        option(_option),
        base(_base) {
    }

    info(int _pos, short _option, char _base): 
       pos(_pos), option(_option) {

        switch(_base) {
            case 'a':
            case 'A':
                base = I_A;
                break;
            case 'c':
            case 'C':
                base = I_C;
                break;
            case 'g':
            case 'G':
                base = I_G;
                break;
            case 't':
            case 'T':
                base = I_T;
                break;
            default:
                base = I_NULL;
        }
    }

    info(){
        pos = I_NULL;
        option = I_NULL;
        base = I_NULL;
    }
};

struct alignmnent_info {
    int value;
    vector<unordered_map<int, info>> infoMap;

    alignmnent_info(){
        value = 0;
    }

    alignmnent_info(int _value):
        value(_value) {
            infoMap.resize(3);
    }
};

vector<unordered_map<int, unordered_map<short, int>>> AlignAll(string &reference, vector<string> &queries, int w, int k, int eps);
vector<unordered_map<int, info>> AlignRegion(string &reference, string &query, mapping_result region);
alignmnent_info align(string &s, string &t, int startPos, int end, int startQ, int endQ);


#endif // ALIGN_H