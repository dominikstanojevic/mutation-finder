package hr.fer.bioinf.stanojevic;

import hr.fer.bioinf.stanojevic.alignment.Alignment;
import hr.fer.bioinf.stanojevic.alignment.Info;
import hr.fer.bioinf.stanojevic.alignment.Option;
import hr.fer.bioinf.stanojevic.mapping.Nucleobase;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static final int W = 10;
    public static final int K = 15;
    public static final int EPS = 500;

    public static void main(String[] args) throws IOException {
        var ref = readFile(Paths.get("data/ecoli.fasta")).get(0);

        var queries = readFile(Paths.get("data/ecoli_simulated_reads.fasta"));
        queries = queries.subList(73, queries.size());
        System.out.println(queries.size());

        Map<Integer, List<Info>> a = Alignment.alignAll(ref, queries.toArray(new String[queries.size()]), W, K, EPS);

        int other = 0;
        for(int i = 0; i < ref.length(); i++, other++) {
            var list = a.get(i);
            if (list == null) {
                System.out.println(i + " " + ref.charAt(other) + " -");
            } else {
                Info info = vote(list);
                if (info == null) {
                } else if (info.option == Option.INSERTION) {
                    other--;
                    System.out.println(i + " - " + info.base.getLetter());
                } else if (info.option == Option.DELETION) {
                    System.out.println(i + " " + ref.charAt(other) + " -");
                } else {
                    if (!(ref.charAt(other) == info.base.getLetter())) {
                        System.out.println(i + " " + ref.charAt(other) + " " + info.base.getLetter());
                    }
                }
            }
        }
    }

    private static Info vote(List<Info> list) {
        if(list.size() == 0) {
            return null;
        }

        int[] bases = new int[4];
        for(int i = 0, e = list.size(); i < e; i++) {
            bases[list.get(i).base.getNumber()] += 1;
        }

        List<Integer> mostCommon = new ArrayList<>();
        mostCommon.add(0);
        int count = bases[0];
        for(int i = 1; i < 4; i++) {
            if (count < bases[i]) {
                mostCommon.clear();
                mostCommon.add(i);
                count = bases[i];
            } else if (count == bases[i]) {
                mostCommon.add(i);
            }
        }

        int base = new Random().nextInt(mostCommon.size());
        base = mostCommon.get(base);
        Nucleobase nb;
        if (base == 0) {
            nb = Nucleobase.ADENINE;
        } else if (base == 1) {
            nb = Nucleobase.CYTOSINE;
        } else if (base == 2) {
            nb = Nucleobase.GUANINE;
        } else {
            nb = Nucleobase.THYMINE;
        }

        Map<Option, Integer> options = new HashMap<>();
        for (Info aList : list) {
            Option o = aList.option;
            if (options.containsKey(o)) {
                options.put(o, options.get(o) + 1);
            } else {
                options.put(o, 1);
            }
        }

        List<Option> mc = new ArrayList<>();
        int c = -1;
        for(var entry : options.entrySet()) {
            if (c == -1) {
                mc.add(entry.getKey());
                c = entry.getValue();
            } else if (c < entry.getValue()) {
                mc.clear();
                mc.add(entry.getKey());
                c = entry.getValue();
            } else if (c == entry.getValue()) {
                mc.add(entry.getKey());
            }
        }

        base = new Random().nextInt(mc.size());

        return new Info(0, mc.get(base), nb);
    }

    public static List<String> readFile(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path);

        List<String> strings = new ArrayList<>();
        StringBuilder sb = new StringBuilder();

        for (String line : lines) {
            line = line.trim();
            if (line.startsWith(">")) {
                String s = sb.toString();
                if (!s.equals("")) {
                    strings.add(s);
                }
                sb = new StringBuilder();
            } else {
                sb.append(line);
            }
        }

        String s = sb.toString();
        if (!s.equals("")) {
            strings.add(s);
        }

        return strings;
    }
}
