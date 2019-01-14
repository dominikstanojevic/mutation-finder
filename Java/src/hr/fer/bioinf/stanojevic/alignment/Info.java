package hr.fer.bioinf.stanojevic.alignment;

import hr.fer.bioinf.stanojevic.mapping.Nucleobase;

public class Info {
    public int pos;
    public Option option;
    public Nucleobase base;

    public Info(int pos, Option option, Nucleobase base) {
        this.pos = pos;
        this.option = option;
        this.base = base;
    }
}
