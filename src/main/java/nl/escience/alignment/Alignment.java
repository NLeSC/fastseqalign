package nl.escience.alignment;


/**
 * Smith-Waterman alignment
 *
 * @author Daniel Cameron
 *
 */
public class Alignment {
    /**
     * score1	the best alignment score
     */
    public final short score;
    /**
     * sub-optimal alignment score
     */
    /**
     * 0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning position is not available
     */
    public final int ref_begin;
    /**
     * 0-based best alignment ending position on reference
     */
    public final int ref_end;
    /**
     * 0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
     position is not available
     */
    public final int read_begin;
    /**
     * 0-based best alignment ending position on read
     */
    public final int read_end;
    /**
     * 0-based sub-optimal alignment ending position on read
     */
    /**
     * best alignment cigar; stored the same as that in BAM format
     */
    //public final String cigar;

    //public final String alignment;

    public Alignment(
            short score,
            int ref_begin,
            int ref_end,
            int	read_begin,
            int read_end
            /*String cigar,
            String alignment*/) {
        this.score = score;
        this.ref_begin = ref_begin;
        this.ref_end = ref_end;
        this.read_begin = read_begin;
        this.read_end = read_end;
        //this.cigar = cigar;
        //this.alignment = alignment;
    }
    @Override
    public String toString() {
            return String.format("score=%d,ref_begin=%d,ref_end=%d,read_begin=%d,read_end=%d",
                score, ref_begin, ref_end, read_begin, read_end);
    }
}