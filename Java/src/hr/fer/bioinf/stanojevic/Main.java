package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.alignment.Alignment;
import hr.fer.bioinf.stanojevic.alignment.Info;
import hr.fer.bioinf.stanojevic.alignment.Option;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static final int W = 3;
    public static final int K = 5;
    public static final int EPS = 500;

    public static void main(String[] args) throws IOException {
        String ref = String.join("", Files.readAllLines(Paths.get(REF_PATH)));
        String query = String.join("", Files.readAllLines(Paths.get(Q_PATH)));


        Map<Integer, List<Info>> a = Alignment.alignAll(ref, new String[]{query}, W, K, EPS);

        int other = 0;
        for(int i = 0; i < ref.length(); i++, other++) {
            var list = a.get(i);
            if (list == null) {
                System.out.println(i + " " + ref.charAt(other) + " -");
            } else {
                var info = list.get(0);
                if (info.option == Option.INSERTION) {
                    other--;
                    System.out.println(i + " - " + info.base.getLetter());
                } else if (info.option == Option.DELETION) {
                    System.out.println(i + " " + ref.charAt(other) + " -");
                    throw new RuntimeException(Integer.toString(i));
                } else {
                    System.out.println(i + " " + ref.charAt(other) + " " + a.get(i).get(0).base.getLetter());
                }
            }
        }
    }
}
