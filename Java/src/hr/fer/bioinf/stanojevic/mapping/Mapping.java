package hr.fer.bioinf.stanojevic.mapping;

import hr.fer.bioinf.stanojevic.Utils;

import java.util.*;

public class Mapping {
    public static Set<Minimizer> minimizerSketch(String s, int w, int k) {
        Set<Minimizer> minimizers = new HashSet<>();

        Map<Integer, Utils.Pair<Long, Long>> cache = new HashMap<>();

        for (int i = 0, end = s.length() - w - k + 1; i <= end; i++) {
            long m = Long.MAX_VALUE;


            for (int j = 0; j < w; j++) {
                int start = i + j;
                var pair = kmerToLongs(s.substring(start, start + k));

                long u = invertibleHash(pair.first);
                long v = invertibleHash(pair.second);

                cache.put(start, new Utils.Pair<>(u, v));

                if (u != v) {
                    long tempMin = Long.compareUnsigned(u, v) < 0 ? u : v;
                    m = Long.compareUnsigned(m, tempMin) < 0 ? m : tempMin;
                }
            }

            for (int j = 0; j < w; j++) {
                int start = i + j;
                var pair = cache.get(start);

                long u = pair.first;
                long v = pair.second;

                if (Long.compareUnsigned(u, v) < 0 && u == m) {
                    minimizers.add(new Minimizer(m, start, 0));
                } else if (Long.compareUnsigned(v, u) < 0 && v == m) {
                    minimizers.add(new Minimizer(m, start, 1));
                }
            }
        }

        return minimizers;
    }

    public static HashMap<Long, Set<IndexData>> index(String[] sequences, int w, int k) {
        HashMap<Long, Set<IndexData>> table = new HashMap<>();

        for (int i = 0; i < sequences.length; i++) {
            var minimizers = minimizerSketch(sequences[i], w, k);
            for (Minimizer m : minimizers) {
                var set = table.computeIfAbsent(m.getHash(), k1 -> new HashSet<>());

                set.add(new IndexData(i, m.getPosition(), m.getReversed()));
            }
        }

        return table;
    }

    public static void map(HashMap<Long, Set<IndexData>> table,
                           String q, int w, int k, int eps) {
        List<MapData> arr = new ArrayList<>();
        Set<Minimizer> minimizers = minimizerSketch(q, w, k);

        for (Minimizer m : minimizers) {
            if (table.containsKey(m.getHash())) {
                for (IndexData ind : table.get(m.getHash())) {
                    if (m.getReversed() == ind.r) {
                        arr.add(new MapData(ind.t, 0, m.getPosition() - ind.i, ind.i));
                    } else {
                        arr.add(new MapData(ind.t, 1, m.getPosition() + ind.i, ind.i));
                    }
                }
            }
        }

        arr.sort(Comparator.naturalOrder());

        int b = 1;
        for (int e = 1, l = arr.size(); e <= l; e++) {
            if (e == l || arr.get(e + 1).t != arr.get(e).t || arr.get(e + 1).r != arr.get(e).r) {
                //TODO
            }
        }
    }



    public static long invertibleHash(long x) {
        x = (~x) + (x << 21);
        x = x ^ (x >> 24);
        x = (x + (x << 3)) + (x << 8);
        x = x ^ (x >> 14);
        x = (x + (x << 2)) + (x << 4);
        x = x ^ (x >> 28);
        x = x + (x << 31);

        return x;
    }

    private static Utils.Pair<Long, Long> kmerToLongs(String s) {
        long kmer = 0;
        long rev = 0;

        for (int i = 0, l = s.length() - 1; i <= l; i++) {
            kmer = (kmer << 2) + Nucleobase.getNucleobase(s.charAt(i)).getNumber();
            rev = (rev << 2) + Nucleobase.getNucleobase(s.charAt(l - i)).getReverseNumber();
        }

        return new Utils.Pair<>(kmer, rev);
    }

    public static class IndexData {
        public int t;
        public int i;
        public int r;

        public IndexData(int t, int i, int r) {
            this.t = t;
            this.i = i;
            this.r = r;
        }
    }

    public static class MapData implements Comparable<MapData>  {
        public int t;
        public int r;
        public int c;
        public int i;

        public MapData(int t, int r, int c, int i) {
            this.t = t;
            this.r = r;
            this.c = c;
            this.i = i;
        }


        @Override
        public int compareTo(MapData other) {
            if (this.t < other.t) return -1;
            if (this.t > other.t) return 1;

            if (this.r < other.r) return -1;
            if (this.r > other.r) return 1;

            if (this.c < other.c) return -1;
            if (this.c > other.c) return 1;

            return Integer.compare(this.i, other.i);
        }
    }
}
