package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.mapping.Mapping;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static final String S = "GTCATGCACGTTCAC";

    public static void main(String[] args) throws IOException {
        String ref = String.join("", Files.readAllLines(Paths.get(REF_PATH)));
        String query = String.join("", Files.readAllLines(Paths.get(Q_PATH)));

        var table = Mapping.index(new String[]{ref}, 3,5);
        Mapping.map(table, query, 3, 5, 500);
    }
}
