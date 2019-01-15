package hr.fer.bioinf.stanojevic.mapping;

public class MappingResult {
    public int start;
    public int end;

    public MappingResult(int start, int end) {
        if (start >= end) {
            throw new RuntimeException("JEBIGA");
        }

        this.start = start;
        this.end = end;
    }
}
