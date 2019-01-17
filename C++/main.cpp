#include <iostream>

#include <algorithm>
#include <cstdlib>
#include <ctime>  
#include <fstream>
#include "align.h"
#include "utils.h"
#include "sketch.h"
#include "map.h"

using namespace std;

int MIN_COVERAGE = 8;

inline short CharToBaseShort(char x){
    switch(x) {
            case 'a':
            case 'A':
                return I_A;
            case 'c':
            case 'C':
                return I_C;
            case 'g':
            case 'G':
                return I_G;
            case 't':
            case 'T':
                return I_T;
            default:
                return I_NULL;
        }
}

inline char BaseToChar(short x){
    switch (x) {
        case I_A:
            return 'A';
        case I_C:
            return 'C';
        case I_G:
            return 'G';
        case I_T:
            return 'T';
        default:
            return 'x';
    }
}

inline char ChangeToChar(short x){
    switch (x) {
        case I_CHA:
            return 'X';
        case I_DEL:
            return 'D';
        case I_INS:
            return 'I';
        default:
            return 'x';
    }
}

pair<short, int> FindMostFrequent(unordered_map<short, int> &occChange){
    int counts[] {occChange[I_A], occChange[I_C], occChange[I_G], occChange[I_T]};

    vector<short> bases;
    bases.push_back(I_A);
    int max = counts[0];

    if (counts[1] > max){
        bases.clear();
        bases.push_back(I_C);
        max = counts[1];
    } else if (counts[1] == max){
        bases.push_back(I_C);
    }

    if (counts[2] > max){
        bases.clear();
        bases.push_back(I_G);
        max = counts[2];
    } else if (counts[2] == max){
        bases.push_back(I_G);
    }

    if (counts[3] > max){
        bases.clear();
        bases.push_back(I_T);
        max = counts[3];
    } else if (counts[3] == max){
        bases.push_back(I_T);
    }

    return make_pair(bases[rand()%bases.size()], max);
}

void FindChanges(string &ref, vector<unordered_map<int, unordered_map<short, int>>> &a, 
        vector<info> &changes){
    
    for(int i = 0; i < ref.size(); ++i) {
        //   base   pos
        pair<short, int> mostFrequentChange;
        if(a[0].count(i) == 0) {
            mostFrequentChange = make_pair(I_NULL, I_NULL);
        } else {
            mostFrequentChange = FindMostFrequent(a[0][i]);
        }

        pair<short, int> mostFrequentIns;
        if (a[1].count(i) == 0) {
            mostFrequentIns = make_pair(I_NULL, I_NULL);
        } else {
            mostFrequentIns = FindMostFrequent(a[1][i]);
            mostFrequentIns.second *= 2;
            //cout << mostFrequentIns.first << ":" << mostFrequentIns.second << endl;
        }

        int delSize = 0;
        if (a[2].count(i) != 0) {
            delSize = a[2][i][-1];
        }

        cout << i << ref[i] << " " << mostFrequentChange.second << ChangeToChar(mostFrequentChange.first) << " ";
        cout << mostFrequentIns.second << " " << delSize << endl;

        if (mostFrequentChange.second <= 0 && mostFrequentIns.second <= 0 && delSize == 0) {
            continue;
        }

        int threshold = mostFrequentIns.second + delSize + (ref[i] == ChangeToChar(mostFrequentChange.first) ? 0 : mostFrequentChange.second);
        if (threshold < MIN_COVERAGE){
            continue;
        }

        if (mostFrequentChange.second >= mostFrequentIns.second && mostFrequentChange.second >= delSize) {
            if (CharToBaseShort(ref[i]) != mostFrequentChange.first) {
                changes.push_back(info(i, I_CHA, mostFrequentChange.first));
            }
        } else if (mostFrequentIns.second > mostFrequentChange.second && mostFrequentIns.second > delSize) {
            changes.push_back(info(i, I_INS, mostFrequentIns.first));
        } else if (delSize > mostFrequentChange.second && delSize > mostFrequentIns.second) {
            if (changes.size() > 0 &&
                (changes.back().pos == i-1) && changes.back().option == I_DEL) {
                changes.pop_back();
            }
            
            changes.push_back(info(i, I_DEL, I_NULL));
        } else if (delSize == mostFrequentIns.second) {
            if (rand() % 2 == 0) {
                changes.push_back(info(i, I_INS, mostFrequentIns.first));
            } else {
                if (changes.size() > 0 && 
                (changes.back().pos == i-1) && changes.back().option == I_DEL) {
                    changes.pop_back();
                }

                changes.push_back(info(i, I_DEL, I_NULL));
            }
        } 
    }
}


/**
 * Input arguments:
 *      - path to the fasta file containing the whole genome
 *      - path to the fasta file containing the chunks of the genome
 */
int main(int argc, char *argv[]){
    srand (time(NULL));
    if (argc != 3){
        cerr << "Input not specified" << endl;
        exit(1);
    }

    int w   =  10;//10;
    int k   =  15;//15;
    int eps = 500;

    string f_genome = argv[1];
    vector<string> v_genome;
    ReadFASTAFile(f_genome, v_genome);
    string genome = v_genome[0];

    string f_reads = argv[2];
    vector<string> v_reads;
    ReadFASTAFile(f_reads, v_reads);

    auto a = AlignAll(genome, v_reads, w, k, eps);

    vector<info> changes;
    FindChanges(genome, a, changes);

    ofstream output ("result.csv");
    if (output.is_open()) {
        char letter;
        for (auto change : changes){
            if (change.option == I_DEL){
                letter = '-';
            } else {
                letter = BaseToChar(change.base);
            }
            output << ChangeToChar(change.option) << "," << change.pos << "," << letter << "\n";
        }
        output.close();
    }
    
    return 0;
}