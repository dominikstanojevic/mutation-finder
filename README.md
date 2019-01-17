# Mutation Finding With Third Generation Sequencing

This project is part of the (Bioinformatics course)[https://www.fer.unizg.hr/en/course/bio] at the Faculty of Electrical Engineering and Computing, University of Zagreb.

## Project setup
In order to compile the project, clone the repository and position yourself in the `C++` folder. Then simply run the `make` command, which should create `mutation_finder`.

## How to run
The program expect 2 arguments:
- path to fasta file containing the referent genome
- path to fasta file containing reads of genome sequences

The program will create a file names `result.csv` that contains all the mutations found in the referent genome.
An example of input files is given in the `data` folder.
