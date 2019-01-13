#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <iostream>

using namespace std;

void ReadFASTAFile(string &file_name, vector<string> &data);
uint64_t Hash(uint64_t x, uint64_t mask);

#endif // UTILS_H