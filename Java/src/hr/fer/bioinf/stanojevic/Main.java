package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.mapping.Mapping;
import hr.fer.bioinf.stanojevic.mapping.Minimizer;

public class Main {
    public static void main(String[] args) {
        String s = "GTCATGCACGTTCAC";

        var minimizers = Mapping.minimizerSketch(s, 3, 3);

        for (Minimizer min : minimizers) {
            System.out.println(min.getHash() + " " + min.getPosition() + " " + min.getReversed());
        }
    }
}
