package hr.fer.bioinf.stanojevic.mapping;

public enum Nucleobase {
    ADENINE('A',  0, 'T'),
    CYTOSINE('C', 1, 'G'),
    GUANINE('G', 2, 'C'),
    THYMINE('T', 3, 'A');

    private char letter;
    private int num;
    private char reverse;

    Nucleobase(char letter, int num, char reverse) {
        this.letter = letter;
        this.num = num;
        this.reverse = reverse;
    }

    public int getNumber() {
        return this.num;
    }

    public int getReverseNumber() {
        return 3 - this.num;
    }

    public char getReverse() {
        return this.reverse;
    }

    public char getLetter() {
        return letter;
    }

    public static Nucleobase getNucleobase(char c) {
        c = Character.toUpperCase(c);
        switch (c) {
            case 'A':
                return ADENINE;
            case 'C':
                return CYTOSINE;
            case 'G':
                return GUANINE;
            case 'T':
                return THYMINE;
            default:
                throw new RuntimeException("Invalid character");
        }
    }
}
