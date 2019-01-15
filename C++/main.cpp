#include <iostream>

#include <algorithm>
#include "align.h"
#include "utils.h"
#include "sketch.h"
#include "map.h"

using namespace std;

/**
 * Input arguments:
 *      - path to the fasta file containing the whole genome
 *      - path to the fasta file containing the chunks of the genome
 */
int main(int argc, char *argv[]){
    if (argc != 3){
        cerr << "Input not specified" << endl;
        exit(1);
    }

    int w   =  3;//10;
    int k   =  5;//15;
    int eps = 500;

    string f_genome = argv[1];
    vector<string> v_genome;
    ReadFASTAFile(f_genome, v_genome);
    string genome = v_genome[0];

    string f_reads = argv[2];
    vector<string> v_reads;
    ReadFASTAFile(f_reads, v_reads);

    auto minimizers = MinimizerSketch(genome, w, k);

    auto table = Index(&genome, 1, w, k);
    Map(table, v_reads[0], w, k, eps);
    return 0;
}