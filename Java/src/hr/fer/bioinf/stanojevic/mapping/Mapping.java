package hr.fer.bioinf.stanojevic.mapping;

import hr.fer.bioinf.stanojevic.Utils;

import java.util.*;

public class Mapping {
    public static final int MIN_COUNT = 4;
    public static final int MIN_READS = 100;

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

            return Integer.compare(o1.c, o2.c);
        });

        int b = 0;
        //System.out.println("Veličina: " + arr.size());
        for (int e = 0, l = arr.size(); e < l; e++) {
            if (e == l - 1 ||
                    arr.get(e + 1).t != arr.get(e).t ||
                    arr.get(e + 1).r != arr.get(e).r ||
                    arr.get(e + 1).c - arr.get(e).c > eps) {

                var subList = arr.subList(b, e + 1);
                //System.out.println(b + " " + (e + 1));

                subList.sort((c1, c2) -> {
                    int q1 = getQPosition(c1.c, c1.i, c1.r);
                    int q2 = getQPosition(c2.c, c2.i, c2.r);

                    int i = Integer.compare(q1, q2);
                    if (i < 0) return -1;
                    if (i > 0) return 1;

                    return Integer.compare(c1.i, c2.i);
                });

                Comparator<Integer> comp = arr.get(e).r == 0 ? Comparator.naturalOrder() : Comparator.naturalOrder();
                MapData[] C = longestIncreasingSubsequence(subList, comp);
                //System.out.println(Arrays.toString(C));
                //System.out.println(b + "    " + (e + 1));

                if(C.length >= MIN_COUNT) {
                    int min = C[0].i;
                    int max = C[C.length - 1].i;

                    if (max - min >= MIN_READS) {
                        System.out.println("Strand: " + arr.get(e).r);
                        System.out.println(min + " " + max);
                        System.out.println();
                    }
                }

                b = e + 1;
            }
        }
    }

    public static int getQPosition(int c, int i, int r) {
        if (r == 0) {
            return i - c;
        } else {

            return c + i;
        }
    }

    public static MapData[] longestIncreasingSubsequence(List<MapData> eAway, Comparator<Integer> comparator) {
        int[] P = new int[eAway.size()];
        int[] M = new int[eAway.size() + 1];

        int L = 0;
        for (int i = 0; i < eAway.size(); i++) {
            int lo = 1, hi = L;

            while (lo <= hi) {
                int mid = (int) Math.ceil((lo + hi) / 2.);
                if (comparator.compare(eAway.get(M[mid]).i, eAway.get(i).i) < 0) {
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
