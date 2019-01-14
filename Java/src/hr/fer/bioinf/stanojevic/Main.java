package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.alignment.Alignment;

import java.io.IOException;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static final int W = 3;
    public static final int K = 5;
    public static final int EPS = 500;

    public static void main(String[] args) throws IOException {
        String s = "ACTCCCCCATTT";
        String t = "GGGCCCGAT";

        var pair = Alignment.align(s, t);
        System.out.println(pair.first + " " + pair.second);
    }
}
