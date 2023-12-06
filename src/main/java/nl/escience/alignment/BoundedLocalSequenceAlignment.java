package nl.escience.alignment;

import java.util.Stack;

/**
 * Implements required functionalities for bounded pseudo-global alignment of
 * two nucleotide/peptide sequences
 *
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, the Netherlands
 */
public class BoundedLocalSequenceAlignment {

    public static int match[][];
    private int matrix[][];
    private char direction[][];
    private int up[][];
    private int left[][];
    private StringBuilder seq1;
    private StringBuilder seq2;
    private StringBuilder cigar;
    private Stack<Character> operation_stack;
    private Stack<Integer> count_stack;
    private int similarity;
    private double identity;
    private int GAP_OPEN;
    private int GAP_EXT;
    private int MAX_LENGTH;
    private int[] score_array;
    private int max_i;
    private int max_j;
    private int deletions;
    private int insertions;
    private int BOUND;
    private char TYPE;
    private int offset;
    private int range_len;
    private int mismatch_penalty, insertion_penalty;
    private int CLIPPING_STRIGENCY;


    /**
     * Initializes the alignment object.
     *
     * @param gap_open The gap opening penalty
     * @param gap_ext The gap extension penalty
     * @param max_length The maximum possible length of the alignment
     * @param band Length of the he band
     * @param clip The stringency of soft-clipping in the range [0..3]
     * @param type Type of the input sequences(N for nucleotide P for peptide).
     */
    public BoundedLocalSequenceAlignment(int gap_open, int gap_ext, int max_length, int band, int clip, char type) {
        int i, j;
        seq1 = new StringBuilder();
        seq2 = new StringBuilder();
        GAP_OPEN = gap_open;
        GAP_EXT = gap_ext;
        MAX_LENGTH = max_length;
        BOUND = band;
        TYPE = type;
        CLIPPING_STRIGENCY = clip;
        cigar = new StringBuilder();
        // initialize matrixes
        matrix = new int[MAX_LENGTH + 1][2 * BOUND + 3];
        direction = new char[MAX_LENGTH + 1][2 * BOUND + 3];
        up = new int[MAX_LENGTH + 1][2 * BOUND + 3];
        left = new int[MAX_LENGTH + 1][2 * BOUND + 3];
        score_array = new int[MAX_LENGTH];
        operation_stack = new Stack();
        count_stack = new Stack();
        direction[0][0] = 'M';
        matrix[0][0] = 0;
        up[0][0] = left[0][0] = -1000;
        for (i = 1; i <= MAX_LENGTH; i++) {
            direction[i][0] = 'M';
            direction[i][2 * BOUND + 2] = 'M';
            // below the bound
            up[i][0] = -1000;
            left[i][0] = -1000;
            matrix[i][0] = 0;
            // above the bound
            up[i][2 * BOUND + 2] = -1000;
            left[i][2 * BOUND + 2] = -1000;
            matrix[i][2 * BOUND + 2] = 0;
        }
        for (j = 1; j <= 2 * BOUND + 2; j++) {
            direction[0][j] = 'D';
            up[0][j] = -1000;
            left[0][j] = -1000;
            matrix[0][j] = 0;
        }
        if (TYPE == 'N')
            initialize_NUCC_matrix();
        else if (TYPE == 'P')
            initialize_BLOSUM_matrix(62);
        else
            System.out.println("Aligner type should be N or P");
        switch (CLIPPING_STRIGENCY) {
            case 1:
                mismatch_penalty = -1;
                insertion_penalty = -1;
                break;
            case 2:
                mismatch_penalty = -4;
                insertion_penalty = -2;
                break;
            case 3:
                mismatch_penalty = -9;
                insertion_penalty = -3;
                break;
        }
    }

    /**
     * Initializes the NUCC.1 scoring matrix.
     */
    public final void initialize_NUCC_matrix() {
        match = new int[256][256];

        match['A']['A'] = 5;
        match['A']['T'] = -4;
        match['A']['G'] = -4;
        match['A']['C'] = -4;
        match['A']['S'] = -4;
        match['A']['W'] = 1;
        match['A']['R'] = 1;
        match['A']['Y'] = -4;
        match['A']['K'] = -4;
        match['A']['M'] = 1;
        match['A']['B'] = -4;
        match['A']['V'] = -1;
        match['A']['H'] = -1;
        match['A']['D'] = -1;
        match['A']['N'] = -2;
        match['T']['A'] = -4;
        match['T']['T'] = 5;
        match['T']['G'] = -4;
        match['T']['C'] = -4;
        match['T']['S'] = -4;
        match['T']['W'] = 1;
        match['T']['R'] = -4;
        match['T']['Y'] = 1;
        match['T']['K'] = 1;
        match['T']['M'] = -4;
        match['T']['B'] = -1;
        match['T']['V'] = -4;
        match['T']['H'] = -1;
        match['T']['D'] = -1;
        match['T']['N'] = -2;
        match['G']['A'] = -4;
        match['G']['T'] = -4;
        match['G']['G'] = 5;
        match['G']['C'] = -4;
        match['G']['S'] = 1;
        match['G']['W'] = -4;
        match['G']['R'] = 1;
        match['G']['Y'] = -4;
        match['G']['K'] = 1;
        match['G']['M'] = -4;
        match['G']['B'] = -1;
        match['G']['V'] = -1;
        match['G']['H'] = -4;
        match['G']['D'] = -1;
        match['G']['N'] = -2;
        match['C']['A'] = -4;
        match['C']['T'] = -4;
        match['C']['G'] = -4;
        match['C']['C'] = 5;
        match['C']['S'] = 1;
        match['C']['W'] = -4;
        match['C']['R'] = -4;
        match['C']['Y'] = 1;
        match['C']['K'] = -4;
        match['C']['M'] = 1;
        match['C']['B'] = -1;
        match['C']['V'] = -1;
        match['C']['H'] = -1;
        match['C']['D'] = -4;
        match['C']['N'] = -2;
        match['S']['A'] = -4;
        match['S']['T'] = -4;
        match['S']['G'] = 1;
        match['S']['C'] = 1;
        match['S']['S'] = -1;
        match['S']['W'] = -4;
        match['S']['R'] = -2;
        match['S']['Y'] = -2;
        match['S']['K'] = -2;
        match['S']['M'] = -2;
        match['S']['B'] = -1;
        match['S']['V'] = -1;
        match['S']['H'] = -3;
        match['S']['D'] = -3;
        match['S']['N'] = -1;
        match['W']['A'] = 1;
        match['W']['T'] = 1;
        match['W']['G'] = -4;
        match['W']['C'] = -4;
        match['W']['S'] = -4;
        match['W']['W'] = -1;
        match['W']['R'] = -2;
        match['W']['Y'] = -2;
        match['W']['K'] = -2;
        match['W']['M'] = -2;
        match['W']['B'] = -3;
        match['W']['V'] = -3;
        match['W']['H'] = -1;
        match['W']['D'] = -1;
        match['W']['N'] = -1;
        match['R']['A'] = 1;
        match['R']['T'] = -4;
        match['R']['G'] = 1;
        match['R']['C'] = -4;
        match['R']['S'] = -2;
        match['R']['W'] = -2;
        match['R']['R'] = -1;
        match['R']['Y'] = -4;
        match['R']['K'] = -2;
        match['R']['M'] = -2;
        match['R']['B'] = -3;
        match['R']['V'] = -1;
        match['R']['H'] = -3;
        match['R']['D'] = -1;
        match['R']['N'] = -1;
        match['Y']['A'] = -4;
        match['Y']['T'] = 1;
        match['Y']['G'] = -4;
        match['Y']['C'] = 1;
        match['Y']['S'] = -2;
        match['Y']['W'] = -2;
        match['Y']['R'] = -4;
        match['Y']['Y'] = -1;
        match['Y']['K'] = -2;
        match['Y']['M'] = -2;
        match['Y']['B'] = -1;
        match['Y']['V'] = -3;
        match['Y']['H'] = -1;
        match['Y']['D'] = -3;
        match['Y']['N'] = -1;
        match['K']['A'] = -4;
        match['K']['T'] = 1;
        match['K']['G'] = 1;
        match['K']['C'] = -4;
        match['K']['S'] = -2;
        match['K']['W'] = -2;
        match['K']['R'] = -2;
        match['K']['Y'] = -2;
        match['K']['K'] = -1;
        match['K']['M'] = -4;
        match['K']['B'] = -1;
        match['K']['V'] = -3;
        match['K']['H'] = -3;
        match['K']['D'] = -1;
        match['K']['N'] = -1;
        match['M']['A'] = 1;
        match['M']['T'] = -4;
        match['M']['G'] = -4;
        match['M']['C'] = 1;
        match['M']['S'] = -2;
        match['M']['W'] = -2;
        match['M']['R'] = -2;
        match['M']['Y'] = -2;
        match['M']['K'] = -4;
        match['M']['M'] = -1;
        match['M']['B'] = -3;
        match['M']['V'] = -1;
        match['M']['H'] = -1;
        match['M']['D'] = -3;
        match['M']['N'] = -1;
        match['B']['A'] = -4;
        match['B']['T'] = -1;
        match['B']['G'] = -1;
        match['B']['C'] = -1;
        match['B']['S'] = -1;
        match['B']['W'] = -3;
        match['B']['R'] = -3;
        match['B']['Y'] = -1;
        match['B']['K'] = -1;
        match['B']['M'] = -3;
        match['B']['B'] = -1;
        match['B']['V'] = -2;
        match['B']['H'] = -2;
        match['B']['D'] = -2;
        match['B']['N'] = -1;
        match['V']['A'] = -1;
        match['V']['T'] = -4;
        match['V']['G'] = -1;
        match['V']['C'] = -1;
        match['V']['S'] = -1;
        match['V']['W'] = -3;
        match['V']['R'] = -1;
        match['V']['Y'] = -3;
        match['V']['K'] = -3;
        match['V']['M'] = -1;
        match['V']['B'] = -2;
        match['V']['V'] = -1;
        match['V']['H'] = -2;
        match['V']['D'] = -2;
        match['V']['N'] = -1;
        match['H']['A'] = -1;
        match['H']['T'] = -1;
        match['H']['G'] = -4;
        match['H']['C'] = -1;
        match['H']['S'] = -3;
        match['H']['W'] = -1;
        match['H']['R'] = -3;
        match['H']['Y'] = -1;
        match['H']['K'] = -3;
        match['H']['M'] = -1;
        match['H']['B'] = -2;
        match['H']['V'] = -2;
        match['H']['H'] = -1;
        match['H']['D'] = -2;
        match['H']['N'] = -1;
        match['D']['A'] = -1;
        match['D']['T'] = -1;
        match['D']['G'] = -1;
        match['D']['C'] = -4;
        match['D']['S'] = -3;
        match['D']['W'] = -1;
        match['D']['R'] = -1;
        match['D']['Y'] = -3;
        match['D']['K'] = -1;
        match['D']['M'] = -3;
        match['D']['B'] = -2;
        match['D']['V'] = -2;
        match['D']['H'] = -2;
        match['D']['D'] = -1;
        match['D']['N'] = -1;
        match['N']['A'] = -2;
        match['N']['T'] = -2;
        match['N']['G'] = -2;
        match['N']['C'] = -2;
        match['N']['S'] = -1;
        match['N']['W'] = -1;
        match['N']['R'] = -1;
        match['N']['Y'] = -1;
        match['N']['K'] = -1;
        match['N']['M'] = -1;
        match['N']['B'] = -1;
        match['N']['V'] = -1;
        match['N']['H'] = -1;
        match['N']['D'] = -1;
        match['N']['N'] = -1;

    }

    /**
     * Initializes a  BLOSUM45, BLOSUM62 or BLOSUM80 scoring matrix.
     *
     * Entries for the BLOSUM80 matrix at a scale of ln(2)/2.0.
     * Source https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM80
     *
     * Entries for the BLOSUM45 matrix at a scale of ln(2)/3.0.
     * Source https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM45
     * @param number 45, 62 or 80
     */
    public static void initialize_BLOSUM_matrix(int number) {
        if (number != 45 && number != 62 && number != 80){
            System.out.println("BLOSUM value must be 45, 62 or 80.");
            System.exit(1);
        }
        match = new int[256][256];
        if (number == 62) {
            match['A']['A'] = 4;
            match['A']['R'] = -1;
            match['A']['N'] = -2;
            match['A']['D'] = -2;
            match['A']['C'] = 0;
            match['A']['Q'] = -1;
            match['A']['E'] = -1;
            match['A']['G'] = 0;
            match['A']['H'] = -2;
            match['A']['I'] = -1;
            match['A']['L'] = -1;
            match['A']['K'] = -1;
            match['A']['M'] = -1;
            match['A']['F'] = -2;
            match['A']['P'] = -1;
            match['A']['S'] = 1;
            match['A']['T'] = 0;
            match['A']['W'] = -3;
            match['A']['Y'] = -2;
            match['A']['V'] = 0;
            match['A']['B'] = -2;
            match['A']['Z'] = -1;
            match['A']['X'] = 0;
            match['A']['*'] = -4;

            match['R']['A'] = -1;
            match['R']['R'] = 5;
            match['R']['N'] = 0;
            match['R']['D'] = -2;
            match['R']['C'] = -3;
            match['R']['Q'] = 1;
            match['R']['E'] = 0;
            match['R']['G'] = -2;
            match['R']['H'] = 0;
            match['R']['I'] = -3;
            match['R']['L'] = -2;
            match['R']['K'] = 2;
            match['R']['M'] = -1;
            match['R']['F'] = -3;
            match['R']['P'] = -2;
            match['R']['S'] = -1;
            match['R']['T'] = -1;
            match['R']['W'] = -3;
            match['R']['Y'] = -2;
            match['R']['V'] = -3;
            match['R']['B'] = -1;
            match['R']['Z'] = 0;
            match['R']['X'] = -1;
            match['R']['*'] = -4;

            match['N']['A'] = -2;
            match['N']['R'] = 0;
            match['N']['N'] = 6;
            match['N']['D'] = 1;
            match['N']['C'] = -3;
            match['N']['Q'] = 0;
            match['N']['E'] = 0;
            match['N']['G'] = 0;
            match['N']['H'] = 1;
            match['N']['I'] = -3;
            match['N']['L'] = -3;
            match['N']['K'] = 0;
            match['N']['M'] = -2;
            match['N']['F'] = -3;
            match['N']['P'] = -2;
            match['N']['S'] = 1;
            match['N']['T'] = 0;
            match['N']['W'] = -4;
            match['N']['Y'] = -2;
            match['N']['V'] = -3;
            match['N']['B'] = 3;
            match['N']['Z'] = 0;
            match['N']['X'] = -1;
            match['N']['*'] = -4;

            match['D']['A'] = -2;
            match['D']['R'] = -2;
            match['D']['N'] = 1;
            match['D']['D'] = 6;
            match['D']['C'] = -3;
            match['D']['Q'] = 0;
            match['D']['E'] = 2;
            match['D']['G'] = -1;
            match['D']['H'] = -1;
            match['D']['I'] = -3;
            match['D']['L'] = -4;
            match['D']['K'] = -1;
            match['D']['M'] = -3;
            match['D']['F'] = -3;
            match['D']['P'] = -1;
            match['D']['S'] = 0;
            match['D']['T'] = -1;
            match['D']['W'] = -4;
            match['D']['Y'] = -3;
            match['D']['V'] = -3;
            match['D']['B'] = 4;
            match['D']['Z'] = 1;
            match['D']['X'] = -1;
            match['D']['*'] = -4;

            match['C']['A'] = 0;
            match['C']['R'] = -3;
            match['C']['N'] = -3;
            match['C']['D'] = -3;
            match['C']['C'] = 9;
            match['C']['Q'] = -3;
            match['C']['E'] = -4;
            match['C']['G'] = -3;
            match['C']['H'] = -3;
            match['C']['I'] = -1;
            match['C']['L'] = -1;
            match['C']['K'] = -3;
            match['C']['M'] = -1;
            match['C']['F'] = -2;
            match['C']['P'] = -3;
            match['C']['S'] = -1;
            match['C']['T'] = -1;
            match['C']['W'] = -2;
            match['C']['Y'] = -2;
            match['C']['V'] = -1;
            match['C']['B'] = -3;
            match['C']['Z'] = -3;
            match['C']['X'] = -2;
            match['C']['*'] = -4;

            match['Q']['A'] = -1;
            match['Q']['R'] = 1;
            match['Q']['N'] = 0;
            match['Q']['D'] = 0;
            match['Q']['C'] = -3;
            match['Q']['Q'] = 5;
            match['Q']['E'] = 2;
            match['Q']['G'] = -2;
            match['Q']['H'] = 0;
            match['Q']['I'] = -3;
            match['Q']['L'] = -2;
            match['Q']['K'] = 1;
            match['Q']['M'] = 0;
            match['Q']['F'] = -3;
            match['Q']['P'] = -1;
            match['Q']['S'] = 0;
            match['Q']['T'] = -1;
            match['Q']['W'] = -2;
            match['Q']['Y'] = -1;
            match['Q']['V'] = -2;
            match['Q']['B'] = 0;
            match['Q']['Z'] = 3;
            match['Q']['X'] = -1;
            match['Q']['*'] = -4;

            match['E']['A'] = -1;
            match['E']['R'] = 0;
            match['E']['N'] = 0;
            match['E']['D'] = 2;
            match['E']['C'] = -4;
            match['E']['Q'] = 2;
            match['E']['E'] = 5;
            match['E']['G'] = -2;
            match['E']['H'] = 0;
            match['E']['I'] = -3;
            match['E']['L'] = -3;
            match['E']['K'] = 1;
            match['E']['M'] = -2;
            match['E']['F'] = -3;
            match['E']['P'] = -1;
            match['E']['S'] = 0;
            match['E']['T'] = -1;
            match['E']['W'] = -3;
            match['E']['Y'] = -2;
            match['E']['V'] = -2;
            match['E']['B'] = 1;
            match['E']['Z'] = 4;
            match['E']['X'] = -1;
            match['E']['*'] = -4;

            match['G']['A'] = 0;
            match['G']['R'] = -2;
            match['G']['N'] = 0;
            match['G']['D'] = -1;
            match['G']['C'] = -3;
            match['G']['Q'] = -2;
            match['G']['E'] = -2;
            match['G']['G'] = 6;
            match['G']['H'] = -2;
            match['G']['I'] = -4;
            match['G']['L'] = -4;
            match['G']['K'] = -2;
            match['G']['M'] = -3;
            match['G']['F'] = -3;
            match['G']['P'] = -2;
            match['G']['S'] = 0;
            match['G']['T'] = -2;
            match['G']['W'] = -2;
            match['G']['Y'] = -3;
            match['G']['V'] = -3;
            match['G']['B'] = -1;
            match['G']['Z'] = -2;
            match['G']['X'] = -1;
            match['G']['*'] = -4;

            match['H']['A'] = -2;
            match['H']['R'] = 0;
            match['H']['N'] = 1;
            match['H']['D'] = -1;
            match['H']['C'] = -3;
            match['H']['Q'] = 0;
            match['H']['E'] = 0;
            match['H']['G'] = -2;
            match['H']['H'] = 8;
            match['H']['I'] = -3;
            match['H']['L'] = -3;
            match['H']['K'] = -1;
            match['H']['M'] = -2;
            match['H']['F'] = -1;
            match['H']['P'] = -2;
            match['H']['S'] = -1;
            match['H']['T'] = -2;
            match['H']['W'] = -2;
            match['H']['Y'] = 2;
            match['H']['V'] = -3;
            match['H']['B'] = 0;
            match['H']['Z'] = 0;
            match['H']['X'] = -1;
            match['H']['*'] = -4;

            match['I']['A'] = -1;
            match['I']['R'] = -3;
            match['I']['N'] = -3;
            match['I']['D'] = -3;
            match['I']['C'] = -1;
            match['I']['Q'] = -3;
            match['I']['E'] = -3;
            match['I']['G'] = -4;
            match['I']['H'] = -3;
            match['I']['I'] = 4;
            match['I']['L'] = 2;
            match['I']['K'] = -3;
            match['I']['M'] = 1;
            match['I']['F'] = 0;
            match['I']['P'] = -3;
            match['I']['S'] = -2;
            match['I']['T'] = -1;
            match['I']['W'] = -3;
            match['I']['Y'] = -1;
            match['I']['V'] = 3;
            match['I']['B'] = -3;
            match['I']['Z'] = -3;
            match['I']['X'] = -1;
            match['I']['*'] = -4;

            match['L']['A'] = -1;
            match['L']['R'] = -2;
            match['L']['N'] = -3;
            match['L']['D'] = -4;
            match['L']['C'] = -1;
            match['L']['Q'] = -2;
            match['L']['E'] = -3;
            match['L']['G'] = -4;
            match['L']['H'] = -3;
            match['L']['I'] = 2;
            match['L']['L'] = 4;
            match['L']['K'] = -2;
            match['L']['M'] = 2;
            match['L']['F'] = 0;
            match['L']['P'] = -3;
            match['L']['S'] = -2;
            match['L']['T'] = -1;
            match['L']['W'] = -2;
            match['L']['Y'] = -1;
            match['L']['V'] = 1;
            match['L']['B'] = -4;
            match['L']['Z'] = -3;
            match['L']['X'] = -1;
            match['L']['*'] = -4;

            match['K']['A'] = -1;
            match['K']['R'] = 2;
            match['K']['N'] = 0;
            match['K']['D'] = -1;
            match['K']['C'] = -3;
            match['K']['Q'] = 1;
            match['K']['E'] = 1;
            match['K']['G'] = -2;
            match['K']['H'] = -1;
            match['K']['I'] = -3;
            match['K']['L'] = -2;
            match['K']['K'] = 5;
            match['K']['M'] = -1;
            match['K']['F'] = -3;
            match['K']['P'] = -1;
            match['K']['S'] = 0;
            match['K']['T'] = -1;
            match['K']['W'] = -3;
            match['K']['Y'] = -2;
            match['K']['V'] = -2;
            match['K']['B'] = 0;
            match['K']['Z'] = 1;
            match['K']['X'] = -1;
            match['K']['*'] = -4;

            match['M']['A'] = -1;
            match['M']['R'] = -1;
            match['M']['N'] = -2;
            match['M']['D'] = -3;
            match['M']['C'] = -1;
            match['M']['Q'] = 0;
            match['M']['E'] = -2;
            match['M']['G'] = -3;
            match['M']['H'] = -2;
            match['M']['I'] = 1;
            match['M']['L'] = 2;
            match['M']['K'] = -1;
            match['M']['M'] = 5;
            match['M']['F'] = 0;
            match['M']['P'] = -2;
            match['M']['S'] = -1;
            match['M']['T'] = -1;
            match['M']['W'] = -1;
            match['M']['Y'] = -1;
            match['M']['V'] = 1;
            match['M']['B'] = -3;
            match['M']['Z'] = -1;
            match['M']['X'] = -1;
            match['M']['*'] = -4;

            match['F']['A'] = -2;
            match['F']['R'] = -3;
            match['F']['N'] = -3;
            match['F']['D'] = -3;
            match['F']['C'] = -2;
            match['F']['Q'] = -3;
            match['F']['E'] = -3;
            match['F']['G'] = -3;
            match['F']['H'] = -1;
            match['F']['I'] = 0;
            match['F']['L'] = 0;
            match['F']['K'] = -3;
            match['F']['M'] = 0;
            match['F']['F'] = 6;
            match['F']['P'] = -4;
            match['F']['S'] = -2;
            match['F']['T'] = -2;
            match['F']['W'] = 1;
            match['F']['Y'] = 3;
            match['F']['V'] = -1;
            match['F']['B'] = -3;
            match['F']['Z'] = -3;
            match['F']['X'] = -1;
            match['F']['*'] = -4;

            match['P']['A'] = -1;
            match['P']['R'] = -2;
            match['P']['N'] = -2;
            match['P']['D'] = -1;
            match['P']['C'] = -3;
            match['P']['Q'] = -1;
            match['P']['E'] = -1;
            match['P']['G'] = -2;
            match['P']['H'] = -2;
            match['P']['I'] = -3;
            match['P']['L'] = -3;
            match['P']['K'] = -1;
            match['P']['M'] = -2;
            match['P']['F'] = -4;
            match['P']['P'] = 7;
            match['P']['S'] = -1;
            match['P']['T'] = -1;
            match['P']['W'] = -4;
            match['P']['Y'] = -3;
            match['P']['V'] = -2;
            match['P']['B'] = -2;
            match['P']['Z'] = -1;
            match['P']['X'] = -2;
            match['P']['*'] = -4;

            match['S']['A'] = 1;
            match['S']['R'] = -1;
            match['S']['N'] = 1;
            match['S']['D'] = 0;
            match['S']['C'] = -1;
            match['S']['Q'] = 0;
            match['S']['E'] = 0;
            match['S']['G'] = 0;
            match['S']['H'] = -1;
            match['S']['I'] = -2;
            match['S']['L'] = -2;
            match['S']['K'] = 0;
            match['S']['M'] = -1;
            match['S']['F'] = -2;
            match['S']['P'] = -1;
            match['S']['S'] = 4;
            match['S']['T'] = 1;
            match['S']['W'] = -3;
            match['S']['Y'] = -2;
            match['S']['V'] = -2;
            match['S']['B'] = 0;
            match['S']['Z'] = 0;
            match['S']['X'] = 0;
            match['S']['*'] = -4;

            match['T']['A'] = 0;
            match['T']['R'] = -1;
            match['T']['N'] = 0;
            match['T']['D'] = -1;
            match['T']['C'] = -1;
            match['T']['Q'] = -1;
            match['T']['E'] = -1;
            match['T']['G'] = -2;
            match['T']['H'] = -2;
            match['T']['I'] = -1;
            match['T']['L'] = -1;
            match['T']['K'] = -1;
            match['T']['M'] = -1;
            match['T']['F'] = -2;
            match['T']['P'] = -1;
            match['T']['S'] = 1;
            match['T']['T'] = 5;
            match['T']['W'] = -2;
            match['T']['Y'] = -2;
            match['T']['V'] = 0;
            match['T']['B'] = -1;
            match['T']['Z'] = -1;
            match['T']['X'] = 0;
            match['T']['*'] = -4;

            match['W']['A'] = -3;
            match['W']['R'] = -3;
            match['W']['N'] = -4;
            match['W']['D'] = -4;
            match['W']['C'] = -2;
            match['W']['Q'] = -2;
            match['W']['E'] = -3;
            match['W']['G'] = -2;
            match['W']['H'] = -2;
            match['W']['I'] = -3;
            match['W']['L'] = -2;
            match['W']['K'] = -3;
            match['W']['M'] = -1;
            match['W']['F'] = 1;
            match['W']['P'] = -4;
            match['W']['S'] = -3;
            match['W']['T'] = -2;
            match['W']['W'] = 11;
            match['W']['Y'] = 2;
            match['W']['V'] = -3;
            match['W']['B'] = -4;
            match['W']['Z'] = -3;
            match['W']['X'] = -2;
            match['W']['*'] = -4;

            match['Y']['A'] = -2;
            match['Y']['R'] = -2;
            match['Y']['N'] = -2;
            match['Y']['D'] = -3;
            match['Y']['C'] = -2;
            match['Y']['Q'] = -1;
            match['Y']['E'] = -2;
            match['Y']['G'] = -3;
            match['Y']['H'] = 2;
            match['Y']['I'] = -1;
            match['Y']['L'] = -1;
            match['Y']['K'] = -2;
            match['Y']['M'] = -1;
            match['Y']['F'] = 3;
            match['Y']['P'] = -3;
            match['Y']['S'] = -2;
            match['Y']['T'] = -2;
            match['Y']['W'] = 2;
            match['Y']['Y'] = 7;
            match['Y']['V'] = -1;
            match['Y']['B'] = -3;
            match['Y']['Z'] = -2;
            match['Y']['X'] = -1;
            match['Y']['*'] = -4;

            match['V']['A'] = 0;
            match['V']['R'] = -3;
            match['V']['N'] = -3;
            match['V']['D'] = -3;
            match['V']['C'] = -1;
            match['V']['Q'] = -2;
            match['V']['E'] = -2;
            match['V']['G'] = -3;
            match['V']['H'] = -3;
            match['V']['I'] = 3;
            match['V']['L'] = 1;
            match['V']['K'] = -2;
            match['V']['M'] = 1;
            match['V']['F'] = -1;
            match['V']['P'] = -2;
            match['V']['S'] = -2;
            match['V']['T'] = 0;
            match['V']['W'] = -3;
            match['V']['Y'] = -1;
            match['V']['V'] = 4;
            match['V']['B'] = -3;
            match['V']['Z'] = -2;
            match['V']['X'] = -1;
            match['V']['*'] = -4;

            match['B']['A'] = -2;
            match['B']['R'] = -1;
            match['B']['N'] = 3;
            match['B']['D'] = 4;
            match['B']['C'] = -3;
            match['B']['Q'] = 0;
            match['B']['E'] = 1;
            match['B']['G'] = -1;
            match['B']['H'] = 0;
            match['B']['I'] = -3;
            match['B']['L'] = -4;
            match['B']['K'] = 0;
            match['B']['M'] = -3;
            match['B']['F'] = -3;
            match['B']['P'] = -2;
            match['B']['S'] = 0;
            match['B']['T'] = -1;
            match['B']['W'] = -4;
            match['B']['Y'] = -3;
            match['B']['V'] = -3;
            match['B']['B'] = 4;
            match['B']['Z'] = 1;
            match['B']['X'] = -1;
            match['B']['*'] = -4;

            match['Z']['A'] = -1;
            match['Z']['R'] = 0;
            match['Z']['N'] = 0;
            match['Z']['D'] = 1;
            match['Z']['C'] = -3;
            match['Z']['Q'] = 3;
            match['Z']['E'] = 4;
            match['Z']['G'] = -2;
            match['Z']['H'] = 0;
            match['Z']['I'] = -3;
            match['Z']['L'] = -3;
            match['Z']['K'] = 1;
            match['Z']['M'] = -1;
            match['Z']['F'] = -3;
            match['Z']['P'] = -1;
            match['Z']['S'] = 0;
            match['Z']['T'] = -1;
            match['Z']['W'] = -3;
            match['Z']['Y'] = -2;
            match['Z']['V'] = -2;
            match['Z']['B'] = 1;
            match['Z']['Z'] = 4;
            match['Z']['X'] = -1;
            match['Z']['*'] = -4;

            match['X']['A'] = 0;
            match['X']['R'] = -1;
            match['X']['N'] = -1;
            match['X']['D'] = -1;
            match['X']['C'] = -2;
            match['X']['Q'] = -1;
            match['X']['E'] = -1;
            match['X']['G'] = -1;
            match['X']['H'] = -1;
            match['X']['I'] = -1;
            match['X']['L'] = -1;
            match['X']['K'] = -1;
            match['X']['M'] = -1;
            match['X']['F'] = -1;
            match['X']['P'] = -2;
            match['X']['S'] = 0;
            match['X']['T'] = 0;
            match['X']['W'] = -2;
            match['X']['Y'] = -1;
            match['X']['V'] = -1;
            match['X']['B'] = -1;
            match['X']['Z'] = -1;
            match['X']['X'] = -1;
            match['X']['*'] = -4;

            match['*']['A'] = -4;
            match['*']['R'] = -4;
            match['*']['N'] = -4;
            match['*']['D'] = -4;
            match['*']['C'] = -4;
            match['*']['Q'] = -4;
            match['*']['E'] = -4;
            match['*']['G'] = -4;
            match['*']['H'] = -4;
            match['*']['I'] = -4;
            match['*']['L'] = -4;
            match['*']['K'] = -4;
            match['*']['M'] = -4;
            match['*']['F'] = -4;
            match['*']['P'] = -4;
            match['*']['S'] = -4;
            match['*']['T'] = -4;
            match['*']['W'] = -4;
            match['*']['Y'] = -4;
            match['*']['V'] = -4;
            match['*']['B'] = -4;
            match['*']['Z'] = -4;
            match['*']['X'] = -4;
            match['*']['*'] = 1;
        }

        if (number == 80) {
            match['A']['A'] = 5;
            match['A']['R'] = -2;
            match['A']['N'] = -2;
            match['A']['D'] = -2;
            match['A']['C'] = -1;
            match['A']['Q'] = -1;
            match['A']['E'] =  -1;
            match['A']['G'] =  0;
            match['A']['H'] = -2;
            match['A']['I'] = -2;
            match['A']['L'] = -2;
            match['A']['K'] = -1;
            match['A']['M'] =  -1;
            match['A']['F'] = -3;
            match['A']['P'] = -1;
            match['A']['S'] = 1;
            match['A']['T'] = 0;
            match['A']['W'] = -3;
            match['A']['Y'] = -2;
            match['A']['V'] = 0;
            match['A']['B'] = -2;
            match['A']['J'] = -2;
            match['A']['Z'] = -1;
            match['A']['X'] = -1;
            match['A']['*'] = -6;

            match['R']['A'] = -2;
            match['R']['R'] = 6;
            match['R']['N'] = -1;
            match['R']['D'] = -2;
            match['R']['C'] = -4;
            match['R']['Q'] = 1;
            match['R']['E'] = -1;
            match['R']['G'] = -3;
            match['R']['H'] = 0;
            match['R']['I'] = -3;
            match['R']['L'] = -3;
            match['R']['K'] = 2;
            match['R']['M'] = -2;
            match['R']['F'] = -4;
            match['R']['P'] = -2;
            match['R']['S'] = -1;
            match['R']['T'] = -1;
            match['R']['W'] = -4;
            match['R']['Y'] = -3;
            match['R']['V'] = -3;
            match['R']['B'] = -1;
            match['R']['J'] = -3;
            match['R']['Z'] = 0;
            match['R']['X'] = -1;
            match['R']['*'] = -6;

            match['N']['A'] = -2;
            match['N']['R'] = -1;
            match['N']['N'] = 6;
            match['N']['D'] = 1;
            match['N']['C'] = -3;
            match['N']['Q'] = 0;
            match['N']['E'] = -1;
            match['N']['G'] = -1;
            match['N']['H'] = 0;
            match['N']['I'] = -4;
            match['N']['L'] = -4;
            match['N']['K'] = 0;
            match['N']['M'] = -3;
            match['N']['F'] = -4;
            match['N']['P'] = -3;
            match['N']['S'] = 0;
            match['N']['T'] = 0;
            match['N']['W'] = -4;
            match['N']['Y'] = -3;
            match['N']['V'] = -4;
            match['N']['B'] = 5;
            match['N']['J'] = -4;
            match['N']['Z'] = 0;
            match['N']['X'] = -1;
            match['N']['*'] = -6;

            match['D']['A'] = -2;
            match['D']['R'] = -2;
            match['D']['N'] = 1;
            match['D']['D'] = 6;
            match['D']['C'] = -4;
            match['D']['Q'] = -1;
            match['D']['E'] = 1;
            match['D']['G'] = -2;
            match['D']['H'] = -2;
            match['D']['I'] = -4;
            match['D']['L'] = -5;
            match['D']['K'] = -1;
            match['D']['M'] = -4;
            match['D']['F'] = -4;
            match['D']['P'] = -2;
            match['D']['S'] = -1;
            match['D']['T'] = -1;
            match['D']['W'] = -6;
            match['D']['Y'] = -4;
            match['D']['V'] = -4;
            match['D']['B'] = 5;
            match['D']['J'] = -5;
            match['D']['Z'] = 1;
            match['D']['X'] = -1;
            match['D']['*'] = -6;

            match['C']['A'] = -1;
            match['C']['R'] = -4;
            match['C']['N'] = -3;
            match['C']['D'] = -4;
            match['C']['C'] = 9;
            match['C']['Q'] = -4;
            match['C']['E'] = -5;
            match['C']['G'] = -4;
            match['C']['H'] = -4;
            match['C']['I'] = -2;
            match['C']['L'] = -2;
            match['C']['K'] = -4;
            match['C']['M'] = -2;
            match['C']['F'] = -3;
            match['C']['P'] =  -4;
            match['C']['S'] = -2;
            match['C']['T'] = -1;
            match['C']['W'] = -3;
            match['C']['Y'] = -3 ;
            match['C']['V'] = -1;
            match['C']['B'] = -4;
            match['C']['J'] = -2 ;
            match['C']['Z'] = -4;
            match['C']['X'] = -1;
            match['C']['*'] = -6;

            match['Q']['A'] = -1;
            match['Q']['R'] = 1;
            match['Q']['N'] = 0;
            match['Q']['D'] = -1;
            match['Q']['C'] = -4;
            match['Q']['Q'] = 6;
            match['Q']['E'] = 2;
            match['Q']['G'] = -2;
            match['Q']['H'] = 1;
            match['Q']['I'] = -3;
            match['Q']['L'] = -3;
            match['Q']['K'] = 1;
            match['Q']['M'] = 0;
            match['Q']['F'] = -4;
            match['Q']['P'] = -2;
            match['Q']['S'] = 0;
            match['Q']['T'] = -1;
            match['Q']['W'] = -3;
            match['Q']['Y'] = -2;
            match['Q']['V'] = -3;
            match['Q']['B'] = 0;
            match['Q']['J'] = -3;
            match['Q']['Z'] = 4;
            match['Q']['X'] = -1;
            match['Q']['*'] = -6;

            match['E']['A'] = -1;
            match['E']['R'] = -1;
            match['E']['N'] = -1;
            match['E']['D'] = 1;
            match['E']['C'] = -5;
            match['E']['Q'] = 2;
            match['E']['E'] = 6;
            match['E']['G'] = -3;
            match['E']['H'] = 0;
            match['E']['I'] = -4;
            match['E']['L'] = -4;
            match['E']['K'] = 1;
            match['E']['M'] = -2;
            match['E']['F'] = -4;
            match['E']['P'] = -2;
            match['E']['S'] = 0;
            match['E']['T'] = -1;
            match['E']['W'] = -4 ;
            match['E']['Y'] = -3;
            match['E']['V'] = -3;
            match['E']['B'] = 1;
            match['E']['J'] = -4;
            match['E']['Z'] = 5;
            match['E']['X'] = -1;
            match['E']['*'] = -6;

            match['G']['A'] = 0;
            match['G']['R'] =  -3;
            match['G']['N'] = -1;
            match['G']['D'] = -2;
            match['G']['C'] = -4;
            match['G']['Q'] =  -2;
            match['G']['E'] = -3;
            match['G']['G'] = 6;
            match['G']['H'] = -3;
            match['G']['I'] =  -5;
            match['G']['L'] = -4;
            match['G']['K'] = -2;
            match['G']['M'] = -4;
            match['G']['F'] = -4;
            match['G']['P'] = -3;
            match['G']['S'] = -1;
            match['G']['T'] = -2;
            match['G']['W'] = -4;
            match['G']['Y'] = -4;
            match['G']['V'] = -4;
            match['G']['B'] = -1;
            match['G']['J'] = -5;
            match['G']['Z'] = -3;
            match['G']['X'] = -1;
            match['G']['*'] = -6;

            match['H']['A'] = -2;
            match['H']['R'] = 0;
            match['H']['N'] = 0;
            match['H']['D'] = -2;
            match['H']['C'] = -4;
            match['H']['Q'] = 1;
            match['H']['E'] = 0;
            match['H']['G'] = -3;
            match['H']['H'] = 8;
            match['H']['I'] = -4;
            match['H']['L'] = -3;
            match['H']['K'] = -1;
            match['H']['M'] = -2;
            match['H']['F'] = -2;
            match['H']['P'] = -3;
            match['H']['S'] = -1;
            match['H']['T'] = -2;
            match['H']['W'] = -3;
            match['H']['Y'] = 2;
            match['H']['V'] = -4;
            match['H']['B'] = -1;
            match['H']['J'] = -4;
            match['H']['Z'] = 0;
            match['H']['X'] = -1;
            match['H']['*'] = -6;

            match['I']['A'] = -2;
            match['I']['R'] = -3;
            match['I']['N'] = -4;
            match['I']['D'] = -4;
            match['I']['C'] =  -2;
            match['I']['Q'] = -3;
            match['I']['E'] = -4;
            match['I']['G'] = -5;
            match['I']['H'] = -4;
            match['I']['I'] = 5;
            match['I']['L'] = 1;
            match['I']['K'] = -3;
            match['I']['M'] = 1;
            match['I']['F'] = -1;
            match['I']['P'] = -4;
            match['I']['S'] = -3;
            match['I']['T'] = -1;
            match['I']['W'] = -3;
            match['I']['Y'] = -2;
            match['I']['V'] = 3;
            match['I']['B'] = -4;
            match['I']['J'] = 3;
            match['I']['Z'] = -4;
            match['I']['X'] = -1;
            match['I']['*'] = -6;

            match['L']['A'] = -2;
            match['L']['R'] = -3;
            match['L']['N'] = -4;
            match['L']['D'] = -5;
            match['L']['C'] = -2;
            match['L']['Q'] = -3;
            match['L']['E'] = -4;
            match['L']['G'] = -4;
            match['L']['H'] = -3;
            match['L']['I'] =  1;
            match['L']['L'] = 4;
            match['L']['K'] = -3;
            match['L']['M'] = 2;
            match['L']['F'] = 0;
            match['L']['P'] = -3;
            match['L']['S'] = -3;
            match['L']['T'] = -2;
            match['L']['W'] = -2;
            match['L']['Y'] = -2;
            match['L']['V'] = 1;
            match['L']['B'] = -4;
            match['L']['J'] = 3;
            match['L']['Z'] = -3;
            match['L']['X'] = -1;
            match['L']['*'] = -6;

            match['K']['A'] = -1;
            match['K']['R'] = 2;
            match['K']['N'] = 0;
            match['K']['D'] = -1;
            match['K']['C'] = -4;
            match['K']['Q'] = 1;
            match['K']['E'] = 1;
            match['K']['G'] = -2;
            match['K']['H'] = -1;
            match['K']['I'] = -3;
            match['K']['L'] = -3;
            match['K']['K'] = 5;
            match['K']['M'] = -2;
            match['K']['F'] = -4;
            match['K']['P'] = -1;
            match['K']['S'] = -1;
            match['K']['T'] = -1;
            match['K']['W'] = -4;
            match['K']['Y'] = -3;
            match['K']['V'] = -3;
            match['K']['B'] = -1;
            match['K']['J'] = -3;
            match['K']['Z'] = 1;
            match['K']['X'] = -1;
            match['K']['*'] = -6;

            match['M']['A'] = -1;
            match['M']['R'] = -2;
            match['M']['N'] = -3;
            match['M']['D'] = -4;
            match['M']['C'] = -2;
            match['M']['Q'] = 0;
            match['M']['E'] = -2;
            match['M']['G'] = -4;
            match['M']['H'] = -2;
            match['M']['I'] = 1;
            match['M']['L'] = 2;
            match['M']['K'] = -2;
            match['M']['M'] = 6;
            match['M']['F'] = 0;
            match['M']['P'] = -3;
            match['M']['S'] = -2;
            match['M']['T'] = -1;
            match['M']['W'] = -2;
            match['M']['Y'] = -2;
            match['M']['V'] = 1;
            match['M']['B'] = -3;
            match['M']['J'] = 2;
            match['M']['Z'] = -1;
            match['M']['X'] = -1;
            match['M']['*'] = -6;

            match['F']['A'] = -3;
            match['F']['R'] = -4;
            match['F']['N'] = -4;
            match['F']['D'] = -4;
            match['F']['C'] = -3;
            match['F']['Q'] = -4;
            match['F']['E'] = -4;
            match['F']['G'] = -4;
            match['F']['H'] = -2;
            match['F']['I'] = -1;
            match['F']['L'] = 0;
            match['F']['K'] = -4;
            match['F']['M'] = 0;
            match['F']['F'] =  6;
            match['F']['P'] =  -4;
            match['F']['S'] = -3;
            match['F']['T'] = -2;
            match['F']['W'] = 0;
            match['F']['Y'] = 3;
            match['F']['V'] = -1;
            match['F']['B'] = -4;
            match['F']['J'] = 0;
            match['F']['Z'] = -4;
            match['F']['X'] = -1;
            match['F']['*'] = -6;

            match['P']['A'] = -1;
            match['P']['R'] = -2;
            match['P']['N'] = -3;
            match['P']['D'] = -2;
            match['P']['C'] = -4 ;
            match['P']['Q'] = -2;
            match['P']['E'] = -2 ;
            match['P']['G'] =  -3;
            match['P']['H'] = -3;
            match['P']['I'] = -4;
            match['P']['L'] = -3 ;
            match['P']['K'] = -1;
            match['P']['M'] = -3;
            match['P']['F'] = -4 ;
            match['P']['P'] = 8;
            match['P']['S'] = -1;
            match['P']['T'] = -2;
            match['P']['W'] = -5;
            match['P']['Y'] = -4;
            match['P']['V'] = -3;
            match['P']['B'] = -2;
            match['P']['J'] = -4;
            match['P']['Z'] = -2;
            match['P']['X'] = -1;
            match['P']['*'] = -6;

            match['S']['A'] = 1;
            match['S']['R'] = -1;
            match['S']['N'] = 0;
            match['S']['D'] = -1;
            match['S']['C'] = -2;
            match['S']['Q'] = 0;
            match['S']['E'] = 0;
            match['S']['G'] = -1;
            match['S']['H'] = -1;
            match['S']['I'] = -3;
            match['S']['L'] = -3;
            match['S']['K'] = -1;
            match['S']['M'] = -2;
            match['S']['F'] = -3;
            match['S']['P'] = -1;
            match['S']['S'] = 5;
            match['S']['T'] = 1;
            match['S']['W'] = -4;
            match['S']['Y'] = -2;
            match['S']['V'] = -2;
            match['S']['B'] = 0;
            match['S']['J'] = -3;
            match['S']['Z'] = 0;
            match['S']['X'] = -1;
            match['S']['*'] = -6;

            match['T']['A'] = 0;
            match['T']['R'] = -1;
            match['T']['N'] = 0;
            match['T']['D'] = -1;
            match['T']['C'] = -1;
            match['T']['Q'] = -1;
            match['T']['E'] = -1;
            match['T']['G'] = -2;
            match['T']['H'] = -2;
            match['T']['I'] = -1;
            match['T']['L'] = -2;
            match['T']['K'] = -1 ;
            match['T']['M'] = -1;
            match['T']['F'] = -2;
            match['T']['P'] = -2;
            match['T']['S'] = 1;
            match['T']['T'] = 5;
            match['T']['W'] = -4;
            match['T']['Y'] = -2;
            match['T']['V'] = 0;
            match['T']['B'] = -1;
            match['T']['J'] = -1;
            match['T']['Z'] = -1 ;
            match['T']['X'] = -1;
            match['T']['*'] = -6;

            match['W']['A'] = -3;
            match['W']['R'] = -4;
            match['W']['N'] = -4;
            match['W']['D'] = -6;
            match['W']['C'] = -3;
            match['W']['Q'] = -3;
            match['W']['E'] = -4;
            match['W']['G'] = -4;
            match['W']['H'] = -3;
            match['W']['I'] =  -3;
            match['W']['L'] = -2;
            match['W']['K'] = -4;
            match['W']['M'] = -2 ;
            match['W']['F'] =  0;
            match['W']['P'] = -5;
            match['W']['S'] = -4;
            match['W']['T'] = -4;
            match['W']['W'] = 11;
            match['W']['Y'] = 2;
            match['W']['V'] = -3;
            match['W']['B'] = -5;
            match['W']['J'] = -3;
            match['W']['Z'] = -3;
            match['W']['X'] = -1;
            match['W']['*'] = -6;

            match['Y']['A'] = -2 ;
            match['Y']['R'] =  -3;
            match['Y']['N'] = -3;
            match['Y']['D'] = -4;
            match['Y']['C'] = -3;
            match['Y']['Q'] = -2;
            match['Y']['E'] = -3;
            match['Y']['G'] = -4;
            match['Y']['H'] = 2;
            match['Y']['I'] = -2;
            match['Y']['L'] = -2;
            match['Y']['K'] = -3;
            match['Y']['M'] = -2;
            match['Y']['F'] = 3;
            match['Y']['P'] = -4;
            match['Y']['S'] =  -2;
            match['Y']['T'] = -2;
            match['Y']['W'] = 2;
            match['Y']['Y'] = 7;
            match['Y']['V'] = -2;
            match['Y']['B'] = -3;
            match['Y']['J'] = -2;
            match['Y']['Z'] = -3;
            match['Y']['X'] = -1;
            match['Y']['*'] = -6;

            match['V']['A'] = 0;
            match['V']['R'] =  -3;
            match['V']['N'] = -4;
            match['V']['D'] = -4;
            match['V']['C'] = -1;
            match['V']['Q'] = -3;
            match['V']['E'] = -3;
            match['V']['G'] = -4;
            match['V']['H'] = -4;
            match['V']['I'] = 3;
            match['V']['L'] = 1;
            match['V']['K'] = -3;
            match['V']['M'] = 1;
            match['V']['F'] = -1;
            match['V']['P'] = -3;
            match['V']['S'] = -2;
            match['V']['T'] = 0;
            match['V']['W'] = -3;
            match['V']['Y'] = -2;
            match['V']['V'] = 4;
            match['V']['B'] = -4;
            match['V']['J'] = 2;
            match['V']['Z'] = -3;
            match['V']['X'] = -1;
            match['V']['*'] = -6;

            match['B']['A'] = -2;
            match['B']['R'] = -1;
            match['B']['N'] = 5;
            match['B']['D'] = 5;
            match['B']['C'] = -4;
            match['B']['Q'] = 0;
            match['B']['E'] = 1;
            match['B']['G'] = -1;
            match['B']['H'] = -1;
            match['B']['I'] = -4;
            match['B']['L'] = -4;
            match['B']['K'] = -1 ;
            match['B']['M'] = -3;
            match['B']['F'] = -4;
            match['B']['P'] = -2;
            match['B']['S'] = 0;
            match['B']['T'] = -1;
            match['B']['W'] = -5 ;
            match['B']['Y'] = -3 ;
            match['B']['V'] = -4;
            match['B']['B'] = 5;
            match['B']['J'] = -4;
            match['B']['Z'] = 0;
            match['B']['X'] = -1;
            match['B']['*'] = -6;

            match['J']['A'] = -2;
            match['J']['R'] = -3;
            match['J']['N'] = -4 ;
            match['J']['D'] = -5;
            match['J']['C'] = -2;
            match['J']['Q'] =  -3;
            match['J']['E'] = -4;
            match['J']['G'] =  -5;
            match['J']['H'] = -4 ;
            match['J']['I'] = 3;
            match['J']['L'] = 3;
            match['J']['K'] = -3;
            match['J']['M'] = 2;
            match['J']['F'] = 0;
            match['J']['P'] = -4;
            match['J']['S'] = -3 ;
            match['J']['T'] = -1;
            match['J']['W'] = -3;
            match['J']['Y'] = -2;
            match['J']['V'] = 2;
            match['J']['B'] = -4;
            match['J']['J'] = 3;
            match['J']['Z'] = -3;
            match['J']['X'] = -1;
            match['J']['*'] = -6;

            match['Z']['A'] = -1;
            match['Z']['R'] = 0;
            match['Z']['N'] = 0;
            match['Z']['D'] = 1;
            match['Z']['C'] = -4;
            match['Z']['Q'] = 4;
            match['Z']['E'] = 5;
            match['Z']['G'] = -3;
            match['Z']['H'] = 0;
            match['Z']['I'] = -4;
            match['Z']['L'] = -3;
            match['Z']['K'] = 1;
            match['Z']['M'] = -1;
            match['Z']['F'] = -4;
            match['Z']['P'] = -2;
            match['Z']['S'] = 0;
            match['Z']['T'] = -1;
            match['Z']['W'] =  -3;
            match['Z']['Y'] = -3;
            match['Z']['V'] = -3;
            match['Z']['B'] = 0;
            match['Z']['J'] = -3;
            match['Z']['Z'] = 5;
            match['Z']['X'] = -1;
            match['Z']['*'] = -6;

            match['X']['A'] = -1;
            match['X']['R'] = -1;
            match['X']['N'] = -1;
            match['X']['D'] = -1;
            match['X']['C'] = -1;
            match['X']['Q'] = -1;
            match['X']['E'] = -1;
            match['X']['G'] = -1;
            match['X']['H'] = -1;
            match['X']['I'] = -1;
            match['X']['L'] = -1;
            match['X']['K'] = -1;
            match['X']['M'] = -1;
            match['X']['F'] = -1;
            match['X']['P'] = -1;
            match['X']['S'] = -1;
            match['X']['T'] = -1;
            match['X']['W'] = -1;
            match['X']['Y'] = -1;
            match['X']['V'] = -1;
            match['X']['B'] = -1;
            match['X']['J'] = -1;
            match['X']['Z'] = -1;
            match['X']['X'] = -1;
            match['X']['*'] = -6;

            match['*']['A'] = -6;
            match['*']['R'] = -6;
            match['*']['N'] = -6;
            match['*']['D'] = -6;
            match['*']['C'] = -6;
            match['*']['Q'] = -6;
            match['*']['E'] = -6;
            match['*']['G'] = -6;
            match['*']['H'] = -6;
            match['*']['I'] = -6;
            match['*']['L'] = -6;
            match['*']['K'] = -6;
            match['*']['M'] = -6;
            match['*']['F'] = -6;
            match['*']['P'] = -6;
            match['*']['S'] = -6;
            match['*']['T'] = -6;
            match['*']['W'] = -6;
            match['*']['Y'] = -6;
            match['*']['V'] = -6;
            match['*']['B'] = -6;
            match['*']['J'] = -6;
            match['*']['Z'] = -6;
            match['*']['X'] = -6;
            match['*']['*'] = 1;
        }

        if (number == 45) {
            match['A']['A'] = 5;
            match['A']['R'] = -2;
            match['A']['N'] = -1;
            match['A']['D'] = -2;
            match['A']['C'] = -1;
            match['A']['Q'] = -1 ;
            match['A']['E'] = -1;
            match['A']['G'] = 0;
            match['A']['H'] = -2;
            match['A']['I'] = -1;
            match['A']['L'] = -1;
            match['A']['K'] = -1;
            match['A']['M'] = -1;
            match['A']['F'] = -2;
            match['A']['P'] = -1;
            match['A']['S'] = 1;
            match['A']['T'] = 0;
            match['A']['W'] = -2;
            match['A']['Y'] = -2;
            match['A']['V'] = 0;
            match['A']['B'] = -1;
            match['A']['J'] = -1;
            match['A']['Z'] = -1;
            match['A']['X'] = -1;
            match['A']['*'] = -5;

            match['R']['A'] = -2;
            match['R']['R'] = 7 ;
            match['R']['N'] = 0;
            match['R']['D'] = -1;
            match['R']['C'] = -3;
            match['R']['Q'] = 1;
            match['R']['E'] = 0;
            match['R']['G'] = -2;
            match['R']['H'] = 0;
            match['R']['I'] = -3;
            match['R']['L'] = -2;
            match['R']['K'] = 3;
            match['R']['M'] = -1;
            match['R']['F'] = -2;
            match['R']['P'] = -2;
            match['R']['S'] = -1;
            match['R']['T'] = -1;
            match['R']['W'] = -2;
            match['R']['Y'] = -1;
            match['R']['V'] = -2;
            match['R']['B'] = -1;
            match['R']['J'] = -3 ;
            match['R']['Z'] = 1;
            match['R']['X'] = -1;
            match['R']['*'] = -5;

            match['N']['A'] =  -1;
            match['N']['R'] = 0;
            match['N']['N'] = 6;
            match['N']['D'] = 2;
            match['N']['C'] = -2;
            match['N']['Q'] = 0;
            match['N']['E'] =  0;
            match['N']['G'] = 0;
            match['N']['H'] = 1;
            match['N']['I'] = -2;
            match['N']['L'] = -3;
            match['N']['K'] =  0;
            match['N']['M'] = -2;
            match['N']['F'] = -2;
            match['N']['P'] =  -2;
            match['N']['S'] = 1;
            match['N']['T'] = 0;
            match['N']['W'] = -4;
            match['N']['Y'] = -2;
            match['N']['V'] = -3;
            match['N']['B'] = 5;
            match['N']['J'] = -3;
            match['N']['Z'] = 0;
            match['N']['X'] = -1;
            match['N']['*'] = -5;

            match['D']['A'] = -2;
            match['D']['R'] = -1;
            match['D']['N'] = 2;
            match['D']['D'] = 7;
            match['D']['C'] = -3;
            match['D']['Q'] = 0;
            match['D']['E'] = 2;
            match['D']['G'] = -1;
            match['D']['H'] = 0;
            match['D']['I'] = -4;
            match['D']['L'] = -3;
            match['D']['K'] = 0;
            match['D']['M'] = -3;
            match['D']['F'] = -4;
            match['D']['P'] = -1;
            match['D']['S'] = 0;
            match['D']['T'] = -1;
            match['D']['W'] = -4;
            match['D']['Y'] = -2;
            match['D']['V'] = -3;
            match['D']['B'] = 6;
            match['D']['J'] = -3;
            match['D']['Z'] = 1;
            match['D']['X'] = -1;
            match['D']['*'] = -5;

            match['C']['A'] =  -1 ;
            match['C']['R'] = -3;
            match['C']['N'] = -2;
            match['C']['D'] =  -3;
            match['C']['C'] = 12;
            match['C']['Q'] = -3 ;
            match['C']['E'] = -3 ;
            match['C']['G'] =  -3;
            match['C']['H'] = -3;
            match['C']['I'] = -3;
            match['C']['L'] = -2 ;
            match['C']['K'] = -3;
            match['C']['M'] = -2;
            match['C']['F'] = -2;
            match['C']['P'] =  -4 ;
            match['C']['S'] = -1;
            match['C']['T'] = -1;
            match['C']['W'] = -5;
            match['C']['Y'] = -3;
            match['C']['V'] = -1;
            match['C']['B'] = -2;
            match['C']['J'] = -2;
            match['C']['Z'] = -3;
            match['C']['X'] = -1;
            match['C']['*'] = -5;

            match['Q']['A'] = -1 ;
            match['Q']['R'] =  1;
            match['Q']['N'] = 0;
            match['Q']['D'] = 0;
            match['Q']['C'] = -3;
            match['Q']['Q'] = 6;
            match['Q']['E'] = 2;
            match['Q']['G'] = -2;
            match['Q']['H'] = 1;
            match['Q']['I'] =  -2;
            match['Q']['L'] = -2 ;
            match['Q']['K'] =  1 ;
            match['Q']['M'] = 0;
            match['Q']['F'] = -4;
            match['Q']['P'] = -1;
            match['Q']['S'] = 0;
            match['Q']['T'] = -1;
            match['Q']['W'] = -2;
            match['Q']['Y'] = -1;
            match['Q']['V'] = -3;
            match['Q']['B'] = 0;
            match['Q']['J'] = -2;
            match['Q']['Z'] = 4;
            match['Q']['X'] = -1;
            match['Q']['*'] = -5;

            match['E']['A'] = 1;
            match['E']['R'] = 0;
            match['E']['N'] = 0;
            match['E']['D'] = 2;
            match['E']['C'] = -3;
            match['E']['Q'] = 2;
            match['E']['E'] = 6;
            match['E']['G'] = -2;
            match['E']['H'] = 0 ;
            match['E']['I'] = -3;
            match['E']['L'] = -2;
            match['E']['K'] = 1;
            match['E']['M'] = -2;
            match['E']['F'] = -3;
            match['E']['P'] = 0;
            match['E']['S'] = 0;
            match['E']['T'] = -1;
            match['E']['W'] = -3;
            match['E']['Y'] = -2;
            match['E']['V'] =  -3;
            match['E']['B'] = 1;
            match['E']['J'] = -3;
            match['E']['Z'] = 5;
            match['E']['X'] = -1;
            match['E']['*'] = -5;

            match['G']['A'] = 0;
            match['G']['R'] = -2;
            match['G']['N'] = 0;
            match['G']['D'] = -1;
            match['G']['C'] = -3;
            match['G']['Q'] = -2;
            match['G']['E'] = -2;
            match['G']['G'] = 7;
            match['G']['H'] = -2;
            match['G']['I'] = -4;
            match['G']['L'] = -3;
            match['G']['K'] = -2;
            match['G']['M'] = -2;
            match['G']['F'] = -3;
            match['G']['P'] = -2;
            match['G']['S'] = 0;
            match['G']['T'] = -2;
            match['G']['W'] = -2;
            match['G']['Y'] = -3;
            match['G']['V'] = -3;
            match['G']['B'] = -1;
            match['G']['J'] = -4;
            match['G']['Z'] = -2;
            match['G']['X'] = -1;
            match['G']['*'] = -5;

            match['H']['A'] = -2 ;
            match['H']['R'] = 0;
            match['H']['N'] = 1;
            match['H']['D'] = 0;
            match['H']['C'] = -3;
            match['H']['Q'] = 1;
            match['H']['E'] = 0;
            match['H']['G'] = -2;
            match['H']['H'] = 10;
            match['H']['I'] = -3;
            match['H']['L'] = -2;
            match['H']['K'] = -1;
            match['H']['M'] = 0;
            match['H']['F'] = -2;
            match['H']['P'] = -2;
            match['H']['S'] = -1;
            match['H']['T'] = -2;
            match['H']['W'] = -3;
            match['H']['Y'] = 2;
            match['H']['V'] = -3;
            match['H']['B'] = 0;
            match['H']['J'] = -2;
            match['H']['Z'] = 0;
            match['H']['X'] = -1;
            match['H']['*'] = -5;

            match['I']['A'] = -1;
            match['I']['R'] = -3;
            match['I']['N'] = -2;
            match['I']['D'] = -4 ;
            match['I']['C'] = -3;
            match['I']['Q'] = -2;
            match['I']['E'] = -3;
            match['I']['G'] = -4;
            match['I']['H'] = -3;
            match['I']['I'] = 5 ;
            match['I']['L'] = 2;
            match['I']['K'] = -3 ;
            match['I']['M'] = 2;
            match['I']['F'] = 0;
            match['I']['P'] = -2 ;
            match['I']['S'] = -2;
            match['I']['T'] = -1;
            match['I']['W'] = -2;
            match['I']['Y'] = 0;
            match['I']['V'] = 3;
            match['I']['B'] = -3;
            match['I']['J'] = 4;
            match['I']['Z'] = -3;
            match['I']['X'] = -1;
            match['I']['*'] = -5;

            match['L']['A'] = -1;
            match['L']['R'] = -2;
            match['L']['N'] = -3;
            match['L']['D'] = -3;
            match['L']['C'] = -2 ;
            match['L']['Q'] = -2;
            match['L']['E'] =  -2;
            match['L']['G'] =  -3 ;
            match['L']['H'] = -2;
            match['L']['I'] = 2;
            match['L']['L'] = 5 ;
            match['L']['K'] = -3;
            match['L']['M'] = 2;
            match['L']['F'] = 1;
            match['L']['P'] =  -3;
            match['L']['S'] = -3;
            match['L']['T'] = -1;
            match['L']['W'] = -2 ;
            match['L']['Y'] = 0;
            match['L']['V'] = 1;
            match['L']['B'] = -3;
            match['L']['J'] = 4;
            match['L']['Z'] = -2;
            match['L']['X'] = -1;
            match['L']['*'] = -5;

            match['K']['A'] = -1;
            match['K']['R'] = 3;
            match['K']['N'] = 0 ;
            match['K']['D'] = 0;
            match['K']['C'] = -3 ;
            match['K']['Q'] = 1;
            match['K']['E'] = 1;
            match['K']['G'] = -2;
            match['K']['H'] = -1 ;
            match['K']['I'] = -3;
            match['K']['L'] = -3;
            match['K']['K'] =  5;
            match['K']['M'] =-1;
            match['K']['F'] = -3 ;
            match['K']['P'] = -1 ;
            match['K']['S'] = -1;
            match['K']['T'] = -1;
            match['K']['W'] = -2 ;
            match['K']['Y'] = -1;
            match['K']['V'] = -2;
            match['K']['B'] = 0;
            match['K']['J'] = -3;
            match['K']['Z'] = 1;
            match['K']['X'] = -1;
            match['K']['*'] = -5;

            match['M']['A'] = -1;
            match['M']['R'] = -1;
            match['M']['N'] = -2 ;
            match['M']['D'] = -3;
            match['M']['C'] = -2;
            match['M']['Q'] = 0;
            match['M']['E'] = -2;
            match['M']['G'] = -2;
            match['M']['H'] = 0;
            match['M']['I'] = 2;
            match['M']['L'] = 2;
            match['M']['K'] = -1 ;
            match['M']['M'] = 6;
            match['M']['F'] = 0;
            match['M']['P'] = -2;
            match['M']['S'] = -2;
            match['M']['T'] = -1;
            match['M']['W'] = -2;
            match['M']['Y'] = 0;
            match['M']['V'] = 1;
            match['M']['B'] = -2;
            match['M']['J'] = 2;
            match['M']['Z'] = -1;
            match['M']['X'] = -1;
            match['M']['*'] = -5;

            match['F']['A'] = -2;
            match['F']['R'] = -2;
            match['F']['N'] = -2;
            match['F']['D'] = -4;
            match['F']['C'] = -2;
            match['F']['Q'] = -4;
            match['F']['E'] = -3;
            match['F']['G'] = -3;
            match['F']['H'] = -2;
            match['F']['I'] = 0;
            match['F']['L'] = 1;
            match['F']['K'] = -3;
            match['F']['M'] = 0;
            match['F']['F'] = 8;
            match['F']['P'] = -3;
            match['F']['S'] = -2;
            match['F']['T'] = -1;
            match['F']['W'] = 1;
            match['F']['Y'] = 3;
            match['F']['V'] = 0;
            match['F']['B'] = -3;
            match['F']['J'] = 1;
            match['F']['Z'] = -3;
            match['F']['X'] = -1;
            match['F']['*'] = -5;

            match['P']['A'] = -1;
            match['P']['R'] = -2;
            match['P']['N'] = -2 ;
            match['P']['D'] = -1;
            match['P']['C'] = -4 ;
            match['P']['Q'] = -1;
            match['P']['E'] = 0;
            match['P']['G'] = -2;
            match['P']['H'] = -2;
            match['P']['I'] = -2;
            match['P']['L'] = -3;
            match['P']['K'] = -1;
            match['P']['M'] = -2;
            match['P']['F'] = -3;
            match['P']['P'] = 9;
            match['P']['S'] = -1;
            match['P']['T'] = -1;
            match['P']['W'] = -3;
            match['P']['Y'] = -3;
            match['P']['V'] = -3;
            match['P']['B'] = -2;
            match['P']['J'] = -3;
            match['P']['Z'] = -1;
            match['P']['X'] = -1;
            match['P']['*'] = -5;

            match['S']['A'] = 1;
            match['S']['R'] = -1;
            match['S']['N'] = 1;
            match['S']['D'] = 0;
            match['S']['C'] = -1;
            match['S']['Q'] = 0;
            match['S']['E'] = 0;
            match['S']['G'] = 0;
            match['S']['H'] = -1;
            match['S']['I'] = -2;
            match['S']['L'] = -3;
            match['S']['K'] = -1;
            match['S']['M'] = -2;
            match['S']['F'] = -2 ;
            match['S']['P'] = -1;
            match['S']['S'] = 4;
            match['S']['T'] = 2;
            match['S']['W'] = -4;
            match['S']['Y'] = -2;
            match['S']['V'] = -1;
            match['S']['B'] = 0;
            match['S']['J'] = -2;
            match['S']['Z'] = 0;
            match['S']['X'] = -1;
            match['S']['*'] = -5;

            match['T']['A'] = 0;
            match['T']['R'] = -1;
            match['T']['N'] = 0;
            match['T']['D'] = -1;
            match['T']['C'] = -1;
            match['T']['Q'] = -1;
            match['T']['E'] = -1;
            match['T']['G'] = -2;
            match['T']['H'] = -2 ;
            match['T']['I'] = -1;
            match['T']['L'] = -1;
            match['T']['K'] = -1;
            match['T']['M'] = -1;
            match['T']['F'] = -1;
            match['T']['P'] = -1;
            match['T']['S'] = 2;
            match['T']['T'] = 5 ;
            match['T']['W'] = -3;
            match['T']['Y'] = -1;
            match['T']['V'] = 0;
            match['T']['B'] = 0;
            match['T']['J'] = -1;
            match['T']['Z'] = -1;
            match['T']['X'] = -1;
            match['T']['*'] = -5;

            match['W']['A'] = -2;
            match['W']['R'] = -2;
            match['W']['N'] = -4;
            match['W']['D'] = -4;
            match['W']['C'] = -5;
            match['W']['Q'] = -2;
            match['W']['E'] = -3;
            match['W']['G'] = -2;
            match['W']['H'] = -3;
            match['W']['I'] = -2;
            match['W']['L'] = -2;
            match['W']['K'] = -2;
            match['W']['M'] = -2;
            match['W']['F'] = 1;
            match['W']['P'] =-3 ;
            match['W']['S'] = -4;
            match['W']['T'] = -3;
            match['W']['W'] = 15;
            match['W']['Y'] = 3;
            match['W']['V'] = -3;
            match['W']['B'] = -4;
            match['W']['J'] = -2;
            match['W']['Z'] = -2;
            match['W']['X'] = -1;
            match['W']['*'] = -5;

            match['Y']['A'] = -2;
            match['Y']['R'] = -1;
            match['Y']['N'] = -2;
            match['Y']['D'] = -2;
            match['Y']['C'] = -3;
            match['Y']['Q'] = -1;
            match['Y']['E'] = -2;
            match['Y']['G'] = -3;
            match['Y']['H'] = 2;
            match['Y']['I'] = 0;
            match['Y']['L'] = 0;
            match['Y']['K'] = -1;
            match['Y']['M'] = 0;
            match['Y']['F'] = 3;
            match['Y']['P'] = -3;
            match['Y']['S'] = -2;
            match['Y']['T'] = -1;
            match['Y']['W'] = 3;
            match['Y']['Y'] = 8;
            match['Y']['V'] = -1;
            match['Y']['B'] = -2;
            match['Y']['J'] = 0;
            match['Y']['Z'] = -2;
            match['Y']['X'] = -1;
            match['Y']['*'] = -5;

            match['V']['A'] = 0;
            match['V']['R'] = -2;
            match['V']['N'] = -3;
            match['V']['D'] = -3;
            match['V']['C'] = -1;
            match['V']['Q'] = -3;
            match['V']['E'] = -3;
            match['V']['G'] = -3;
            match['V']['H'] = -3;
            match['V']['I'] = 3;
            match['V']['L'] = 1;
            match['V']['K'] = -2;
            match['V']['M'] = 1;
            match['V']['F'] = 0;
            match['V']['P'] = -3;
            match['V']['S'] = -1;
            match['V']['T'] = 0;
            match['V']['W'] = -3;
            match['V']['Y'] = -1;
            match['V']['V'] = 5;
            match['V']['B'] = -3;
            match['V']['J'] = 2;
            match['V']['Z'] = -3;
            match['V']['X'] = -1;
            match['V']['*'] = -5;

            match['B']['A'] = -1;
            match['B']['R'] = -1;
            match['B']['N'] = 5;
            match['B']['D'] = 6;
            match['B']['C'] = -2;
            match['B']['Q'] = 0;
            match['B']['E'] = 1;
            match['B']['G'] = -1;
            match['B']['H'] = 0;
            match['B']['I'] = -3;
            match['B']['L'] = -3;
            match['B']['K'] = 0;
            match['B']['M'] = -2;
            match['B']['F'] = -3;
            match['B']['P'] = -2;
            match['B']['S'] = 0;
            match['B']['T'] = 0;
            match['B']['W'] = -4;
            match['B']['Y'] = -2;
            match['B']['V'] = -3;
            match['B']['B'] = 5;
            match['B']['J'] = -3;
            match['B']['Z'] = 1;
            match['B']['X'] = -1;
            match['B']['*'] = -5;

            match['J']['A'] = -1;
            match['J']['R'] = -3;
            match['J']['N'] = -3;
            match['J']['D'] = -3;
            match['J']['C'] = -2;
            match['J']['Q'] =  -2;
            match['J']['E'] = -3;
            match['J']['G'] =  -4;
            match['J']['H'] = -2;
            match['J']['I'] = 4;
            match['J']['L'] = 4;
            match['J']['K'] = -3;
            match['J']['M'] =  2;
            match['J']['F'] =  1 ;
            match['J']['P'] = -3;
            match['J']['S'] = -2;
            match['J']['T'] = -1;
            match['J']['W'] = -2;
            match['J']['Y'] = 0;
            match['J']['V'] = 2 ;
            match['J']['B'] = -3 ;
            match['J']['J'] = 4;
            match['J']['Z'] = -2;
            match['J']['X'] = -1;
            match['J']['*'] = -5;

            match['Z']['A'] = -1;
            match['Z']['R'] = 1;
            match['Z']['N'] = 0;
            match['Z']['D'] = 1;
            match['Z']['C'] = -3;
            match['Z']['Q'] = 4;
            match['Z']['E'] = 5;
            match['Z']['G'] = -2;
            match['Z']['H'] = 0;
            match['Z']['I'] = -3;
            match['Z']['L'] = -2;
            match['Z']['K'] = 1;
            match['Z']['M'] = -1;
            match['Z']['F'] = -3;
            match['Z']['P'] = -1;
            match['Z']['S'] = 0;
            match['Z']['T'] = -1;
            match['Z']['W'] = -2 ;
            match['Z']['Y'] = -2;
            match['Z']['V'] = -3;
            match['Z']['B'] = 1;
            match['Z']['J'] = -2;
            match['Z']['Z'] = 5;
            match['Z']['X'] = -1;
            match['Z']['*'] = -5;

            match['X']['A'] = -1;
            match['X']['R'] = -1;
            match['X']['N'] = -1;
            match['X']['D'] = -1;
            match['X']['C'] = -1;
            match['X']['Q'] = -1;
            match['X']['E'] = -1;
            match['X']['G'] = -1;
            match['X']['H'] = -1;
            match['X']['I'] = -1;
            match['X']['L'] = -1;
            match['X']['K'] = -1;
            match['X']['M'] = -1;
            match['X']['F'] = -1;
            match['X']['P'] = -1;
            match['X']['S'] = -1;
            match['X']['T'] = -1 ;
            match['X']['W'] = -1;
            match['X']['Y'] = -1;
            match['X']['V'] = -1;
            match['X']['B'] = -1;
            match['X']['J'] = -1;
            match['X']['Z'] = -1;
            match['X']['X'] = -1;
            match['X']['*'] = -5;

            match['*']['A'] = -5;
            match['*']['R'] = -5;
            match['*']['N'] = -5;
            match['*']['D'] = -5;
            match['*']['C'] = -5;
            match['*']['Q'] = -5;
            match['*']['E'] = -5;
            match['*']['G'] = -5;
            match['*']['H'] = -5;
            match['*']['I'] = -5;
            match['*']['L'] = -5;
            match['*']['K'] = -5;
            match['*']['M'] = -5;
            match['*']['F'] = -5;
            match['*']['P'] = -5;
            match['*']['S'] = -5;
            match['*']['T'] = -5;
            match['*']['W'] = -5;
            match['*']['Y'] = -5;
            match['*']['V'] = -5;
            match['*']['B'] = -5;
            match['*']['J'] = -5;
            match['*']['Z'] = -5;
            match['*']['X'] = -5;
            match['*']['*'] = 1;
        }
    }

    /**
     * Fills the the similarity and direction matrixes of the two input sequences.
     * First sequence should not be longer than the second sequence.
     *
     * @param s1 The String containing the first sequence
     * @param s2 The String containing the second sequence
     */
    public Alignment align(String s1, String s2) {
        int i, j, stop;
        int m, n, d;
        seq1.setLength(0);
        seq1.append(s1);
        seq2.setLength(0);
        seq2.append(s2);
        m = seq1.length();
        n = seq2.length();
        if (m < MAX_LENGTH) {
            similarity = Integer.MIN_VALUE;
            for (i = 1; i <= m; i++) {
                stop = 2 * BOUND + 1;
                for (j = 1; j <= stop; j++) {
                    up[i][j] = Math.max( up[i-1][j+1] - GAP_EXT , Math.max(matrix[i-1][j+1], left[i-1][j+1]) - GAP_OPEN - GAP_EXT);
                    left[i][j] = Math.max( left[i][j-1] - GAP_EXT , Math.max(matrix[i][j-1], up[i][j-1]) - GAP_OPEN - GAP_EXT);
                    d = match[seq1.charAt(i-1)][seq2.charAt(j+i-2)] + matrix[i-1][j];
                    if (d >= Math.max( up[i][j] , left[i][j])) {
                        matrix[i][j] = d;
                        direction[i][j] = 'M';
                    } else if (left[i][j] > up[i][j]) {
                        matrix[i][j] = left[i][j];
                        direction[i][j] = 'D';
                    } else {
                        matrix[i][j] = up[i][j];
                        direction[i][j] = 'I';
                    }
                    if (matrix[i][j] > similarity) {
                        similarity = matrix[i][j];
                        max_i = i;
                        max_j = j;
                    }
                    //System.out.print(String.format("%4d", matrix[i][j] ));
                    //System.out.print(String.format("%4c",direction[i][j]));
                    //System.out.print(String.format("%4d",left[i][j]));
                    //System.out.print(String.format("%4d",up[i][j]));
                }
                //System.out.println();
            }
        } else {
            System.err.println("Sequences are too large for the aligner.");
            System.exit(0);
        }

        return (new Alignment((short)similarity, -1, max_j, -1, max_i));
        /*System.out.println("m: " + m + " n: " + n);
        System.out.println("Coordinates = "+ max_i + " " + max_j);
        System.out.println(this.get_alignment());
        System.out.println(this.get_cigar());
        System.out.println("offset: "+ this.get_offset());
        System.out.println(similarity_score+ "\n");*/
    }

    /**
     * Calculates the alignment from the similarity matrix.
     * Call align() before this function.
     *
     * @return the alignment string which may contains some gaps.
     */
    public String get_alignment() {
        int i, j;
        int range[] = new int[]{1,max_i, max_j};
        StringBuilder subject = new StringBuilder();
        StringBuilder query = new StringBuilder();
        subject.setLength(0);
        query.setLength(0);
        i = max_i;
        j = max_j;
        if (CLIPPING_STRIGENCY > 0)
            range = calculate_clip_range();
        while (i > 0 && j > 0) {
            if (CLIPPING_STRIGENCY > 0 && i < range[0])
                break;
            if (direction[i][j] == 'I') {
                query.append( seq1.charAt(i-1) );
                subject.append( '-' );
                i = i - 1;
                j = j + 1;
            } else if (direction[i][j] == 'D') {
                query.append( '-' );
                subject.append( seq2.charAt(j+i-2) );
                j = j - 1;
            } else {
                query.append( seq1.charAt(i-1) );
                subject.append( seq2.charAt(j+i-2) );
                i = i - 1;
            }
        }
        if (CLIPPING_STRIGENCY > 0) {
            for (;i > 0 && j > 1; --i, --j) {
                query.append( seq1.charAt(i-1) );
                subject.append( seq2.charAt(j+i-2) );
            }
        }
        for (;i > 0; --i) {
            query.append( seq1.charAt(i-1) );
            subject.append( '-' );
        }
        for (;j > 1; --j) {
            query.append( '-' );
            subject.append( seq2.charAt(i+j-2) );
        }
        return subject.reverse() + "\n" + query.reverse();
    }

    /**
     * Calculates the similarity score of the shorter protein with itself.
     *
     * @return The similarity score of the shorter protein with itself.
     */
    public long perfect_score() {
        char ch;
        int i;
        long score = 0;
        for (i = 0; i < seq1.length(); ++i) {
            ch = seq1.charAt(i);
            score += match[ch][ch];
        }
        return score;
    }

    /**
     * Calculates the score of un-gapped alignment of two sequences.
     *
     * @param s1 The first sequence.
     * @param s2 The second sequence.
     * @return The score of un-gapped alignment of two sequences.
     */
    public long get_match_score(String s1, String s2) {
        int i;
        long score = 0;
        for (i = 0; i < s1.length(); ++i) {
            score += match[s1.charAt(i)][s2.charAt(i)];
        }
        return score;
    }

    /**
     * Calculates the score of un-gapped alignment of two sequences.
     *
     * @param s1 The first sequence.
     * @param s2 The second sequence.
     * @return The score of un-gapped alignment of two sequences.
     */
    public long get_match_score(StringBuilder s1, StringBuilder s2) {
        int i;
        long score = 0;
        for (i = 0; i < s1.length(); ++i) {
            score += match[s1.charAt(i)][s2.charAt(i)];
        }
        return score;
    }

    /**
     * Calculates the score of un-gapped alignment of two sequences in percentage.
     *
     * @param s1 The first sequence.
     * @param s2 The second sequence.
     * @return The score of un-gapped alignment of two sequences in percentage.
     */
    public double get_match_percentage(String s1, String s2) {
        int i;
        long score = 0, p_score = 0;
        for (i = 0; i < s1.length(); ++i) {
            score += match[s1.charAt(i)][s2.charAt(i)];
            p_score += match[s1.charAt(i)][s1.charAt(i)];
        }
        return score * 100.0 / p_score;
    }

    /**
     * @return The similarity score of two sequences, as the largest entry in
     *         similarity matrix.
     */
    public int get_similarity() {
        return similarity;
    }

    /**
     * @return The identity of two sequences after being aligned.
     */
    public double get_identity() {
        return identity;
    }


    /**
     * Calculates the boundaries for soft-clipping of the alignment.
     *
     * @return array [from1, to1, from2, to2] containing the clipping boundaries
     *         in the first and second sequences aligned.
     */
    public int[] calculate_clip_range() {
        int i, j, x, max_ending_here, max_so_far, tmp_start, tmp_stop;
        int range[] = new int[4];
        x = i = max_i;
        j = max_j + max_i - 1;
        while (i > 0 && j > 0) {
            if (direction[i][j - i + 1] == 'I') {
                score_array[x--] = insertion_penalty;
                i = i - 1;
            } else if (direction[i][j - i + 1] == 'D') {
                j = j - 1;
            } else {
                score_array[x--] = seq1.charAt(i-1) == seq2.charAt(j-1) ? 1 : mismatch_penalty;
                i = i - 1;
                j = j - 1;
            }
        }
        for (;i > 0; --i)
            score_array[x--] = 0;
        max_ending_here = max_so_far = score_array[1];
        tmp_start = tmp_stop = 1;
        range[0] = range[1] = 1;
        for (i = 2; i <= max_i; ++i) {
            if (score_array[i] > max_ending_here + score_array[i]) {
                max_ending_here = score_array[i];
                tmp_start = tmp_stop = i;
            } else {
                max_ending_here = max_ending_here + score_array[i];
                tmp_stop = i;
            }
            if (max_so_far < max_ending_here) {
                range[0] = tmp_start;
                range[1] = tmp_stop;
                max_so_far = max_ending_here;
            }
        }
        i = max_i;
        j = max_j + max_i - 1;
        while (i > 0 && j > 0) {
            if (i == range[1]) {
                range[3] = j;
            }
            if (i == range[0]) {
                range[2] = j;
            }
            if (direction[i][j - i + 1] == 'I') {
                i = i - 1;
            } else if (direction[i][j - i + 1] == 'D') {
                j = j - 1;
            } else {
                i = i - 1;
                j = j - 1;
            }
        }
        return range;
    }

    /**
     * Calculates the SAM cigar string of the alignment.
     *
     * @return  The SAM cigar string of the alignment.
     */
    public String get_cigar() {
        int i, j, move_counts, count, identicals = 0;
        int range[];
        char curr_move, prev_move, operation;
        insertions = deletions = 0;
        operation_stack.clear();
        count_stack.clear();
        cigar.setLength(0);
        if (CLIPPING_STRIGENCY > 0) {
            range = calculate_clip_range();
            if (seq1.length() - range[1] > 0) {
                operation_stack.push('S');
                count_stack.push(seq1.length() - range[1]);
            }
            prev_move = 'M';
            move_counts = 1;
        } else {
            range = new int[]{1, max_i, 1, max_j + max_i - 1};
            prev_move = 'M';
            move_counts = range[1] < seq1.length() ? seq1.length() - range[1] + 1 : 1;
        }
        range_len = range[1] - range[0] + 1;
        i = seq1.length();
        j = range[3] + seq1.length() - range[1];
        if (j > seq2.length()) {
            i -= j - seq2.length();
            j -= j - seq2.length();
        }
        for (; i >= range[1]; --i, --j)
            if (seq1.charAt(i-1) == seq2.charAt(j-1))
                identicals++;
        i = range[1] - 1;
        j = range[3] - 1;
        while (i >= range[0]) {
            curr_move = direction[i][j - i + 1];
            if (curr_move == 'I') {
                i = i - 1;
                ++insertions;
            } else if (curr_move == 'D') {
                j = j - 1;
                ++deletions;
            } else {
                if (seq1.charAt(i) == seq2.charAt(j))
                    ++identicals;
                i = i - 1;
                j = j - 1;
            }
            if (prev_move == curr_move)
                ++move_counts;
            else{
                operation_stack.push(prev_move);
                count_stack.push(move_counts);
                move_counts = 1;
            }
            prev_move = curr_move;
            //System.out.println(i+" "+j+ " " +direction[i][j]);
        }
        offset = j - i;
        if (CLIPPING_STRIGENCY > 0) {
            operation_stack.push(prev_move);
            count_stack.push(move_counts);
            if (i > 0) {
                operation_stack.push('S');
                count_stack.push(i);
                offset += i;
            }
        } else {
            if (prev_move == 'M')
                move_counts += i;
            operation_stack.push(prev_move);
            count_stack.push(move_counts);
        }
        while (!operation_stack.isEmpty()) {
            operation = operation_stack.pop();
            count = count_stack.pop();
            cigar.append(count).append(operation);
        }
        for (; i > 0 && j > 0; --i, --j)
            if (seq1.charAt(i-1) == seq2.charAt(j-1))
                identicals++;
        identity = ((double)identicals) / (seq1.length() + deletions);
        return cigar.toString();
    }

    /**
     * Calculates offset as the number of gaps at start of the first sequence
     * after alignment.
     * Can be negative if there are gaps at the start of the second sequence.
     * Should be called only after calling get_cigar().
     *
     * @return The offset as the number of gaps at start of the first sequence
     * after alignment.
     */
    public int get_offset() {
        return offset;
    }

    /**
     * @return The length of the alignment after soft-clipping
     */
    public int get_range_length() {
        return range_len;
    }

    /**
     * @return The total number of insertions in the alignment
     */
    public int get_insertions() {
        return insertions;
    }

    /**
     * @return The total number of deletions in the alignment
     */
    public int get_deletions() {
        return deletions;
    }
}
