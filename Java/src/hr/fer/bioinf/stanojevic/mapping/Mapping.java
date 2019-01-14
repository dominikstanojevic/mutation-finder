package hr.fer.bioinf.stanojevic.mapping;

import hr.fer.bioinf.stanojevic.Utils;

import java.util.*;

public class Mapping {
    public static Set<Minimizer> minimizerSketch(String s, int w, int k) {
        Set<Minimizer> minimizers = new LinkedHashSet<>();
        long mask = (1L << (2*k)) - 1;

        for (int i = 0, end = s.length() - w - k + 1; i <= end; i++) {
            long m = Long.MAX_VALUE;

            Map<Integer, Utils.Pair<Long, Long>> cache = new HashMap<>();
            for (int j = 0; j < w; j++) {
                int start = i + j;
                var pair = kmerToLongs(s.substring(start, start + k));

                long u = invertibleHash(pair.first, mask);
                long v = invertibleHash(pair.second, mask);

                cache.put(start, new Utils.Pair<>(u, v));

                if (!pair.first.equals(pair.second)) {
                    if (Long.compareUnsigned(pair.first, pair.second) < 0)
                        m = Long.compareUnsigned(m, u) < 0 ? m : u;
                    else
                        m = Long.compareUnsigned(m, v) < 0 ? m : v;
                }
            }

            for (int j = 0; j < w; j++) {
                int start = i + j;
                var pair = cache.get(start);

                long u = pair.first;
                long v = pair.second;

                if (u == m) {
                    minimizers.add(new Minimizer(m, start, 0));
                } else if (v == m) {
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
                long bucket = m.getHash();
                var set = table.computeIfAbsent(bucket, k1 -> new HashSet<>());

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
            long bucket = m.getHash();
            if (table.containsKey(bucket)) {
                for (IndexData ind : table.get(bucket)) {
                    if (m.getReversed() == ind.r) {
                        arr.add(new MapData(ind.t, 0, ind.i - m.getPosition(), ind.i));
                        //System.out.println(ind.r + " " + (ind.i - m.getPosition()) + " " + (k - 1 + ind.i));
                    } else {
                        arr.add(new MapData(ind.t, 1, m.getPosition() + ind.i, ind.i));
                        //System.out.println(ind.r + " " + (ind.i + m.getPosition()) + " " + (k - 1 + ind.i));
                    }
                }
            }
        }


        arr.sort((o1, o2) -> {
            int r = Integer.compare(o1.r, o2.r);
            if (r != 0) {
                return r;
            }

            int i = Integer.compare(o1.c, o2.c);
            if (i != 0) {
                return i;
            }

            return Integer.compare(o1.i, o2.i);
        });


        for (MapData md : arr) {
            System.out.println(md.r + " " + md.c + " " + md.i);
        }

        int b = 0;
        //System.out.println("Veličina: " + arr.size());
        for (int e = 0, l = arr.size(); e < l; e++) {
            if (e == l - 1 ||
                    arr.get(e + 1).t != arr.get(e).t ||
                    arr.get(e + 1).r != arr.get(e).r ||
                    arr.get(e + 1).c - arr.get(e).c > eps) {
                //TODO Ovo ide nakon svega
                if (e + 1 - b < 4) {
                    continue;
                }

                var subList = arr.subList(b, e + 1);
                //System.out.println(b + " " + (e + 1));


                MapData[] C = longestIncreasingSubsequence(subList);
                //System.out.println(b + "    " + (e + 1));

                int min = subList.stream().mapToInt(c -> c.i).min().getAsInt();
                int max = subList.stream().mapToInt(c -> c.i).max().getAsInt();
                //System.out.println(min + " " + max);

                b = e + 1;
            }
        }
    }

    public static MapData[] longestIncreasingSubsequence(List<MapData> eAway) {
        int[] P = new int[eAway.size()];
        int[] M = new int[eAway.size() + 1];

        int L = 0;
        for (int i = 0; i < eAway.size(); i++) {
            int lo = 1, hi = L;

            while (lo <= hi) {
                int mid = (int) Math.ceil((lo + hi) / 2.);
                if (eAway.get(M[mid]).c < eAway.get(i).c) {
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }

            int newL = lo;
            P[i] = M[newL - 1];
            M[newL] = i;

            if (newL > L) {
                L = newL;
            }
        }

        MapData[] S = new MapData[L];
        int k = M[L];
        for (int i = L - 1; i >= 0; i--) {
            S[i] = eAway.get(k);
            k = P[k];
        }

        return S;
    }

    public static long invertibleHash(long x, long mask) {
        x = ((~x) + (x << 21)) & mask;
        x = x ^ (x >>> 24);
        x = ((x + (x << 3)) + (x << 8)) & mask;
        x = x ^ (x >>> 14);
        x = ((x + (x << 2)) + (x << 4)) & mask;
        x = x ^ (x >>> 28);
        x = (x + (x << 31)) & mask;

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
            if (this.c < other.c) return -1;
            if (this.c > other.c) return 1;

            return Integer.compare(this.i, other.i);
        }

        @Override
        public String toString() {
            return "(" + t + ", " + r + ", " + c + ", " + i + ")";
        }
    }
}
