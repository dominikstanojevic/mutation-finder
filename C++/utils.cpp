#include "utils.h"

#include <iostream>
#include <fstream>

using namespace std;

void ReadFASTAFile(string &file_name, vector<string> &data){
    ifstream file(file_name);
    if (file.is_open()) {
        string line;
        string genome;
        while (!file.eof()){
            getline(file, line);
            if (line[0] == '>') {   // Description line
                if (genome.size() > 0) {
                    data.push_back(genome);
                    genome = "";
                }
                continue;
            }

            genome += line;
        }
        data.push_back(genome);
    } else {
        cerr << "File not found!" << endl;
        exit(1);
    }
    file.close();
}

uint64_t Hash(uint64_t x, uint64_t mask){
    x = (~x + (x << 21)) & mask;
    x = x ^ x >> 24;
    x = ((x + (x << 3)) + (x << 8)) & mask;
    x = x ^ x >> 14;
    x = ((x + (x << 2)) + (x << 4)) & mask;
    x = x ^ x >> 28;
    x = (x + (x << 31)) & mask;
    return x;
}