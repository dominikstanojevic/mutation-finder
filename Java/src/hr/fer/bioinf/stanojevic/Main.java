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
import java.util.stream.Collectors;

public class Main {

    public static final String REF_PATH = "data/test.fa";
    public static final String Q_PATH = "data/query.fq";

    public static final int W = 10;
    public static final int K = 15;
    public static final int EPS = 500;

    public static final Random RANDOM = new Random();

    public static int max = -1;

    public static void main(String[] args) throws IOException {
        /*var ref = readFile(Paths.get("data/ecoli.fasta")).get(0);

        var queries = readFile(Paths.get("data/ecoli_simulated_reads.fasta"));
        max = queries.stream().mapToInt(q -> q.length()).max().getAsInt();
        System.out.println(queries.size());
*/
        //var a = Alignment.alignAll(ref, queries.toArray(new String[queries.size()]), W, K, EPS);
        //var a = Alignment.alignAll(ref, new String[]{queries.get(1847)}, W, K, EPS);

        String ref = String.join("", Files.readAllLines(Paths.get(REF_PATH)));
        String query = String.join("", Files.readAllLines(Paths.get(Q_PATH)));

        var a = Alignment.alignAll(ref, new String[]{query}, W, K, EPS);

        var changes = findChanges(ref, a);
        StringJoiner sj = new StringJoiner("\n");

        for(Info change : changes) {
            char letter = change.option == Option.DELETION ? '-' : change.base.getLetter();
            sj.add(change.option.getLetter() + "," + change.pos + "," + letter);
        }

        Files.write(Paths.get("data/result.csv"), sj.toString().getBytes());
    }

    private static List<Info> findChanges(String ref, Map<Integer, List<Info>>[] a) {
        List<Info> changes = new ArrayList<>();
        for(int i = 0, end = ref.length(); i < end; i++) {
            List<Info> change = a[0].get(i);
            List<Info> ins = a[1].get(i);
            List<Info> del = a[2].get(i);

            Utils.Pair<Nucleobase, Integer> mostFrequentChange;
            if(change == null) {
                mostFrequentChange = new Utils.Pair<>(null, -1);
            } else {
                var occChange = change.stream().collect(Collectors.groupingBy(b -> b.base, Collectors.counting()));
                mostFrequentChange = findMostFrequent(occChange);
            }

            Utils.Pair<Nucleobase, Integer> mostFrequentIns;
            if (ins == null) {
                mostFrequentIns = new Utils.Pair<>(null, -1);
            } else {
                var occIns = ins.stream().collect(Collectors.groupingBy(b -> b.base, Collectors.counting()));
                mostFrequentIns = findMostFrequent(occIns);
            }

            int delSize = 0;
            if (del != null) {
                delSize = del.size();
            }

            if (mostFrequentChange.second <= 0 && mostFrequentIns.second <= 0 && delSize == 0) {
                continue;
            }

            if (mostFrequentChange.second >= mostFrequentIns.second && mostFrequentChange.second >= delSize) {
                if (ref.charAt(i) != mostFrequentChange.first.getLetter()) {
                    changes.add(new Info(i, Option.CHANGE, mostFrequentChange.first));
                }
            } else if (mostFrequentIns.second > mostFrequentChange.second && mostFrequentIns.second > delSize) {
                changes.add(new Info(i, Option.INSERTION, mostFrequentIns.first));
            } else if (delSize > mostFrequentChange.second && del.size() > mostFrequentIns.second) {
                changes.add(new Info(i, Option.DELETION, null));
            } else if (delSize == mostFrequentIns.second) {
                if (RANDOM.nextBoolean()) {
                    changes.add(new Info(i, Option.INSERTION, mostFrequentIns.first));
                } else {
                    changes.add(new Info(i, Option.DELETION, null));
                }
            } else {
                throw new RuntimeException("Da li se može ovo dogoditi?");
            }
        }

        return changes;
    }

    private static Utils.Pair<Nucleobase, Integer> findMostFrequent(Map<Nucleobase, Long> occ) {
        long[] counts = new long[4];
        counts[0] = occ.containsKey(Nucleobase.ADENINE) ? occ.get(Nucleobase.ADENINE) : 0;
        counts[1] = occ.containsKey(Nucleobase.CYTOSINE) ? occ.get(Nucleobase.CYTOSINE) : 0;
        counts[2] = occ.containsKey(Nucleobase.GUANINE) ? occ.get(Nucleobase.GUANINE) : 0;
        counts[3] = occ.containsKey(Nucleobase.THYMINE) ? occ.get(Nucleobase.THYMINE) : 0;

        List<Nucleobase> bases = new ArrayList<>();
        bases.add(Nucleobase.ADENINE);
        long max = counts[0];

        for(int i = 1; i < 4; i++) {
            if (counts[i] > max) {
                bases.clear();
                bases.add(Nucleobase.getNucleobase(i));
                max = counts[i];
            } else if (counts[i] == max) {
                bases.add(Nucleobase.getNucleobase(i));
            }
        }

        if (max == 0) {
            return new Utils.Pair<>(null, -1);
        }

        int bound = bases.size();
        int pos = RANDOM.nextInt(bound);

        return new Utils.Pair<>(bases.get(pos), (int) max);
    }

    public static void main2(String[] args) {
        String ref = "AATCGACTAGGACGGTAACGCATTGAGAGTT";
        String query = "AACGTTAGAGTT";

        var result = Alignment.align(ref, query);
        System.out.println(result.aligned.first + " " + " " + result.aligned.second);
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
