package hr.fer.bioinf.stanojevic;

public class Utils {
    public static class Pair<T, U> {
        public T first;
        public U second;

        public Pair(T first, U second) {
            this.first = first;
            this.second = second;
        }
    }

    public static class Triplet<T, U, V> {
        public T first;
        public U second;
        public V third;

        public Triplet(T first, U second, V third) {
            this.first = first;
            this.second = second;
            this.third = third;
        }
    }
}
