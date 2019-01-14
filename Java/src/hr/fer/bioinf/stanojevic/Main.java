package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.mapping.Mapping;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static final int W = 3;
    public static final int K = 5;
    public static final int EPS = 500;

    public static void main(String[] args) throws IOException {
        String ref = String.join("", Files.readAllLines(Paths.get(REF_PATH)));
        String query = String.join("", Files.readAllLines(Paths.get(Q_PATH)));

        var table = Mapping.index(new String[]{ref}, W,K);
        Mapping.map(table, query, W, K, EPS);
    }
}
