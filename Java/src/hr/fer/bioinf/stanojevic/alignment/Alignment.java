package hr.fer.bioinf.stanojevic.alignment;

import hr.fer.bioinf.stanojevic.Utils;
import hr.fer.bioinf.stanojevic.mapping.Mapping;
import hr.fer.bioinf.stanojevic.mapping.MappingResult;
import hr.fer.bioinf.stanojevic.mapping.Nucleobase;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Alignment {
    public static final int MATCH = 4;
    public static final int DIFF = -1;
    public static final int EMPTY = -2;

    public static Map<Integer, List<Info>> alignAll(String reference, String[] queries, int w, int k, int eps) {
        Map<Integer, List<Info>> info = new HashMap<>();

        var table = Mapping.index(new String[]{reference}, w, k);
        for (String query : queries) {
            var regions = Mapping.map(table, query, w, k, eps);
            var result = Alignment.alignQuery(reference, query, regions);

            for (var entry : result.entrySet()) {
                var list = info.computeIfAbsent(entry.getKey(), s -> new ArrayList<>());
                list.add(entry.getValue());
            }
        }

        return info;
    }

    public static Map<Integer, Info> alignQuery(String reference, String query, List<MappingResult> regions) {
        Map<Integer, Info> info = new HashMap<>();

        AlignmentInfo best = null;
        int bestStart = -1;
        for (MappingResult region : regions) {
            int len = region.end - region.start;

            int start = region.start - (query.length() - len);
            start = start > 0 ? start : 0;

            int end = region.end + (query.length() - len);
            end = end < reference.length() ? end : reference.length();

            AlignmentInfo result = align(reference.substring(start, end), query);
            if (best == null || best.value < result.value) {
                best = result;
                bestStart = start;
            }
        }

        var aligned = best.aligned;

        int i = 0;
        while (aligned.second.charAt(i) == '-') {
            i++;
        }

        int j = query.length() - 1;
        while (aligned.second.charAt(j) == '-') {
            j--;
        }

        for(; i < j; i++) {
            int pos = bestStart + i;
            if (reference.charAt(pos) == '-' && query.charAt(i) == '-') {
                System.out.println("Pička ti materina");
            } else if (reference.charAt(pos) == '-') {
                info.put(pos, new Info(pos, Option.INSERTION, Nucleobase.getNucleobase(query.charAt(i))));
            } else if (query.charAt(i) == '-') {
                info.put(pos, new Info(pos, Option.DELETION, null));
            } else {
                info.put(pos, new Info(pos, Option.CHANGE, Nucleobase.getNucleobase(query.charAt(i))));
            }
        }

        return info;
    }


    private static AlignmentInfo align(String s, String t) {
        int rows = s.length();
        int cols = t.length();

        int[][] array = new int[rows + 1][cols + 1];

        for(int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                int match = array[i - 1][j - 1] + (s.charAt(i - 1) == t.charAt(j - 1) ? MATCH : DIFF);
                int delete = array[i - 1][j] + EMPTY;
                int insert = array[i][j - 1] + EMPTY;

                array[i][j] = Math.max(Math.max(match, delete), insert);
            }
        }

        return traceback(array, s, t);
    }

    private static AlignmentInfo traceback(int[][] arr, String s, String t) {
        int rows = arr.length - 1;
        int cols = arr[0].length - 1;

        int i = rows;
        int j = 0;
        int maxValue = arr[rows][0];

        for(int e = 1; e <= cols; e++) {
            if (arr[rows][e] > maxValue) {
                j = e;
                maxValue = arr[rows][j];
            }
        }

        for(int e = 0; e <= rows; e++) {
            if (arr[e][cols] > maxValue) {
                i = e;
                maxValue = arr[i][cols];
            }
        }

        StringBuilder aS, aT;
        if (i == rows && j == cols) {
            aS = new StringBuilder();
            aT = new StringBuilder();
        } else if (i == rows) {
            aS = new StringBuilder();
            for(int e = 0, end = cols - j; e < end; e++) {
                aS.append('-');
            }

            aT = new StringBuilder(t.substring(j, cols));
        } else {
            aS = new StringBuilder(s.substring(i, rows));

            aT = new StringBuilder();
            for(int e = 0, end = rows - i; e < end; e++) {
                aT.append('-');
            }
        }

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 &&
                    arr[i][j] == arr[i - 1][j - 1] + (s.charAt(i - 1) == t.charAt(j - 1) ? MATCH : DIFF)) {
                aS.insert(0, s.charAt(i - 1));
                aT.insert(0, t.charAt(j - 1));

                i--; j--;
            } else if (i > 0) {
                aS.insert(0, s.charAt(i - 1));
                aT.insert(0, '-');

                i--;
            } else {
                aS.insert(0, '-');
                aT.insert(0, t.charAt(j - 1));

                j--;
            }
        }

        return new AlignmentInfo(maxValue, new Utils.Pair<>(aS.toString(), aT.toString()));
    }

    private static class AlignmentInfo {
        public int value;
        public Utils.Pair<String, String> aligned;

        public AlignmentInfo(int value, Utils.Pair<String, String> aligned) {
            this.value = value;
            this.aligned = aligned;
        }
    }
}
