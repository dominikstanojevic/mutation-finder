package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.mapping.Mapping;
import hr.fer.bioinf.stanojevic.mapping.Minimizer;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static void main(String[] args) throws IOException {
        String ref = String.join("", Files.readAllLines(Paths.get(REF_PATH)));
        String query = String.join("", Files.readAllLines(Paths.get(Q_PATH)));

        var minimizers = Mapping.minimizerSketch(ref, 10, 15);


        Mapping.map(Mapping.index(new String[]{ref}, 10, 15), query, 10, 15, 500);

        for (Minimizer min : minimizers) {
            System.out.println(min.getHash() + " " + min.getPosition() + " " + min.getReversed());
        }

    }
}
