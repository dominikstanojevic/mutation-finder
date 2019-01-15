package hr.fer.bioinf.stanojevic.alignment;

public enum Option {
    CHANGE("X"),
    INSERTION("I"),
    DELETION("D");

    private String letter;

    Option(String letter) {
        this.letter = letter;
    }

    public String getLetter() {
        return letter;
    }
}
