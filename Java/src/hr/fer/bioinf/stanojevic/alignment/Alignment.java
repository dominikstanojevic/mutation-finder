package hr.fer.bioinf.stanojevic.alignment;

import hr.fer.bioinf.stanojevic.Main;
import hr.fer.bioinf.stanojevic.Utils;
import hr.fer.bioinf.stanojevic.mapping.Mapping;
import hr.fer.bioinf.stanojevic.mapping.MappingResult;
import hr.fer.bioinf.stanojevic.mapping.Nucleobase;

import java.util.*;

public class Alignment {
    public static final int MATCH = 4;
    public static final int DIFF = -1;
    public static final int EMPTY = -2;

    public static int counter = 0;

    public static int[][] array;

    public static Map<Integer, List<Info>>[] alignAll(String reference, String[] queries, int w, int k, int eps) {
        System.out.println(2 * Main.max);
        array = new int[30_000][30_000];

        //0 CHANGE 1 INSERTION 2 DELETION
        Map<Integer, List<Info>>[] info = new HashMap[3];
        info[0] = new HashMap<>();
        info[1] = new HashMap<>();
        info[2] = new HashMap<>();

        var table = Mapping.index(new String[]{reference}, w, k);
        for (String query : queries) {
            var regions = Mapping.map(table, query, w, k, eps);
            var result = Alignment.alignQuery(reference, query, regions);

            for(int i = 0; i < 3; i++) {
                for (var entry : result[i].entrySet()) {
                    var list = info[i].computeIfAbsent(entry.getKey(), s -> new ArrayList<>());
                    list.add(entry.getValue());
                }
            }
        }

        return info;
    }

    public static Map<Integer, Info>[] alignQuery(String reference, String query, List<MappingResult> regions) {
        //0 CHANGE 1 INSERTION 2 DELETION
        Map<Integer, Info>[] info = new HashMap[3];
        info[0] = new HashMap<>();
        info[1] = new HashMap<>();
        info[2] = new HashMap<>();

        AlignmentInfo best = null;
        int bestStart = -1;
        for (MappingResult region : regions) {
            int len = region.end - region.start;
            if (len < 0) {
                throw new RuntimeException("Puče");
            }

            int diff = query.length() - len;

            int start = diff > 0 ? (region.start - diff) : region.start;
            start = start > 0 ? start : 0;

            int end = diff > 0 ? (region.end + diff) : region.end;
            end = end < reference.length() ? end : reference.length();

            AlignmentInfo result = align(reference.substring(start, end), query);
            if (best == null || best.value < result.value) {
                best = result;
                bestStart = start;
            }
        }

        /*for(int i = 0; i < best.aligned.first.length(); i++) {
            System.out.println(i + " " + best.aligned.first.charAt(i) + " " + best.aligned.second.charAt(i));
        }
*/

        if (best == null) {
            return info;
        }
        var aligned = best.aligned;

        int i = 0;
        while (aligned.second.charAt(i) == '-') {
            i++;
        }

        int j = aligned.second.length() - 1;
        while (aligned.second.charAt(j) == '-') {
            j--;
        }

        int pos = bestStart + i;
        for(; i < j; i++) {
            if (aligned.second.charAt(i) == '-') {
                info[2].put(pos, new Info(pos, Option.DELETION, null));
                pos++;
            } else if (aligned.first.charAt(i) == '-') {
                info[1].put(pos, new Info(pos, Option.INSERTION, Nucleobase.getNucleobase(aligned.second.charAt(i))));
            } else {
                info[0].put(pos, new Info(pos, Option.CHANGE, Nucleobase.getNucleobase(aligned.second.charAt(i))));
                pos++;
            }
        }

        System.out.println("Gotov " + counter++);
        return info;
    }


    public static AlignmentInfo align(String s, String t) {
        int rows = s.length();
        int cols = t.length();

        for (int j = 1; j <= cols; j++) {
            array[0][j] = j * EMPTY;
        }

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
        int rows = s.length();
        int cols = t.length();

        int i = 0;
        int j = cols;
        int maxValue = arr[0][cols];

        for(int e = 1; e <= rows; e++) {
            if (arr[e][cols] > maxValue) {
                i = e;
                maxValue = arr[e][cols];
            }
        }

        StringBuilder aS, aT;
        aS = new StringBuilder(s.substring(i, rows));
        aT = new StringBuilder();
        for(int e = 0, end = rows - i; e < end; e++) {
                aT.append('-');
        }

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 &&
                    arr[i][j] == arr[i - 1][j - 1] + (s.charAt(i - 1) == t.charAt(j - 1) ? MATCH : DIFF)) {
                i--; j--;

                aS.insert(0, s.charAt(i));
                aT.insert(0, t.charAt(j));
            } else if (j > 0 && arr[i][j] == arr[i][j - 1] + EMPTY) {
                aS.insert(0, '-');
                aT.insert(0, t.charAt(j - 1));

                j--;
            } else {
                aS.insert(0, s.charAt(i - 1));
                aT.insert(0, '-');

                i--;
            }

           /* if (i > 0 && j > 0 &&
                    arr[i][j] == arr[i - 1][j - 1] + (s.charAt(i - 1) == t.charAt(j - 1) ? MATCH : DIFF)) {
                i--; j--;

                aS.insert(0, s.charAt(i));
                aT.insert(0, t.charAt(j));
            } else if (i > 0) {
                aS.insert(0, s.charAt(i - 1));
                aT.insert(0, '-');

                i--;
            } else {
                aS.insert(0, '-');
                aT.insert(0, t.charAt(j - 1));

                j--;
            }*/
        }

        return new AlignmentInfo(maxValue, new Utils.Pair<>(aS.toString(), aT.toString()));
    }

    public static class AlignmentInfo {
        public int value;
        public Utils.Pair<String, String> aligned;

        public AlignmentInfo(int value, Utils.Pair<String, String> aligned) {
            this.value = value;
            this.aligned = aligned;
        }
    }
}
