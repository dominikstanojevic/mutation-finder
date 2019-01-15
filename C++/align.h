#ifndef ALIGN_H
#define ALIGN_H

#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "map.h"

using namespace std;

struct info {
    int pos;
    char option;
    char base;

    info(int _pos, char _option, char _base): 
        pos(_pos),
        option(_option),
        base(_base) {
    }

    info(){
        pos = 0;
        option = '\0';
        base = '\0';
    }
};

struct alignmnent_info {
    int value;
    string alX;
    string alY;

    alignmnent_info(){
        value = 0;
        alX = "";
        alY = "";
    }

    alignmnent_info(int _value, string _alX, string _alY):
        value(_value),
        alX(_alX),
        alY(_alY) {

    }
};

map<int, vector<info> > AlignAll(string &reference, vector<string> &queries, int w, int k, int eps);
map<int, info> AlignQuery(string &reference, string &query, vector<mapping_result> regions);
alignmnent_info align(string &s, string &t);

#endif // ALIGN_H