package hr.fer.bioinf.stanojevic.mapping;

import java.util.Objects;

public class Minimizer implements Comparable<Minimizer> {
    private long hash;
    private int position;
    private int reversed;

    public Minimizer(long hash, int position, int reversed) {
        this.hash = hash;
        this.position = position;
        this.reversed = reversed;
    }

    @Override
    public int compareTo(Minimizer minimizer) {
        if (this.hash < minimizer.hash) return -1;
        if (this.hash > minimizer.hash) return 1;

        if (this.position < minimizer.position) return -1;
        if (this.position > minimizer.position) return 1;

        return Integer.compare(this.reversed, minimizer.reversed);
    }

    public long getHash() {
        return hash;
    }

    public int getPosition() {
        return position;
    }

    public int getReversed() {
        return reversed;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Minimizer minimizer = (Minimizer) o;
        return hash == minimizer.hash &&
                position == minimizer.position &&
                reversed == minimizer.reversed;
    }

    @Override
    public int hashCode() {
        return Objects.hash(hash, position, reversed);
    }
}
