package hr.fer.bioinf.stanojevic.alignment;

import hr.fer.bioinf.stanojevic.Utils;
import hr.fer.bioinf.stanojevic.mapping.MappingResult;

import java.util.List;

public class Alignment {
    public static final int MATCH = 4;
    public static final int DIFF = -1;
    public static final int EMPTY = -2;

    public static void alignQuery(String reference, String query, List<MappingResult> regions) {
        for (MappingResult region : regions) {
            int len = region.end - region.start;
            int start = region.start - (query.length() - len);
            int end = region.end + (query.length() - len);

            var resultS = align(reference.substring(start, end), query);

        }
    }


    public static Utils.Pair<String, String> align(String s, String t) {
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

    private static Utils.Pair<String, String> traceback(int[][] arr, String s, String t) {
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

        return new Utils.Pair<>(aS.toString(), aT.toString());
    }
}