#include <iostream>

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

    string f_genome = argv[1];
    vector<string> v_genome;
    ReadFASTAFile(f_genome, v_genome);
    string genome = v_genome[0];

    string f_reads = argv[2];
    vector<string> v_reads;
    ReadFASTAFile(f_reads, v_reads);

    vector<minimizer> minimizers = MinimizerSketch(genome, 10, 15);

    cout << minimizers.size() << endl;

    auto table = Index(&genome, 1, 10, 15);
    Map(table, v_reads[0], 10, 15, 500);
    return 0;
}