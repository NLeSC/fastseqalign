package nl.escience.alignment;

import javax.swing.plaf.synth.SynthTextAreaUI;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Stack;


/**
 * Implements required functionalities for pseudo-global alignment of
 * two nucleotide/peptide sequences 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, the Netherlands
 */
public class LocalSequenceAlignment {

    private static int match[][];
    private int matrix[][];
    private byte direction[][];
    //private char direction[][];
    private int up[][];
    private int left[][];
    private StringBuilder seq1;
    private StringBuilder seq2;
    private StringBuilder cigar;
    private Stack<Character> operation_stack;
    private Stack<Integer> count_stack;
    private int MAX_LENGTH;
    private static int GAP_OPEN;
    private static int GAP_EXT;
    private int max_i;
    private int max_j;
    private int deletions;
    private int insertions;
    private int similarity;
    private double identity;
    private char TYPE;
    private int offset;
    private int range_len;
    private int[] score_array;
    private int mismatch_penalty, insertion_penalty;
    private int CLIPPING_STRIGENCY;
    private boolean DEBUG = false;

    private int align_count;

    private File seq1_file;
    private File seq2_file;
    private File align_stats_file;

    //private FileWriter seq1_file_writer;

    //private FileWriter seq2_file_writer;

    //private FileWriter align_stats_file_writer;
    
    /**
     * Initializes the alignment object. 
     * 
     * @param gap_open The gap opening penalty
     * @param gap_ext The gap extension penalty 
     * @param max_length The maximum possible length of the alignment
     * @param clip The stringency of soft-clipping in the range [0..3]
     * @param type Type of the input sequences(N for nucleotide P for peptide). 
     */
    public LocalSequenceAlignment(int gap_open, int gap_ext, int max_length, int clip, char type){
        int i, j;
        seq1 = new StringBuilder();
        seq2 = new StringBuilder();
        MAX_LENGTH = max_length;
        GAP_OPEN = gap_open;
        GAP_EXT = gap_ext;
        CLIPPING_STRIGENCY = clip;
        TYPE = type;
        cigar = new StringBuilder();
    // initialize matrixes
        matrix = new int[MAX_LENGTH+1][MAX_LENGTH+1];
        direction = new byte[MAX_LENGTH + 1][MAX_LENGTH + 1];
        //direction = new char[MAX_LENGTH + 1][MAX_LENGTH + 1];
        up = new int[MAX_LENGTH+1][MAX_LENGTH+1];
        left = new int[MAX_LENGTH+1][MAX_LENGTH+1];
        score_array = new int[MAX_LENGTH];
        operation_stack = new Stack();
        count_stack = new Stack();
        //direction[0][0] = 'M';
        direction[0][0] = 0;
        matrix[0][0] = 0;
        up[0][0] = left[0][0] = -1000;

        try {
            seq1_file = new File("seq1_file.fasta");
        }
        catch (NullPointerException e){
            System.err.println(e.getMessage());
        }
        try {
            seq2_file = new File("seq2_file.fasta");
        }
        catch (NullPointerException e){
            System.err.println(e.getMessage());
        }

        try {
            align_stats_file = new File("align_stats_file.txt");
        }
        catch (NullPointerException e){
            System.err.println(e.getMessage());
        }

        /*try {
            seq1_file_writer = new FileWriter(seq1_file, true);
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }

        try {
            seq2_file_writer = new FileWriter(seq2_file, true);
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }
        try {
            align_stats_file_writer = new FileWriter(align_stats_file, true);
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }*/

        align_count = 0;




        for (i = 1; i <= MAX_LENGTH; i++) {
                up[i][0] = -1000;
                left[i][0] = -1000;
                matrix[i][0] = 0;
                //direction[i][0] = 'I';
                direction[i][0] = 0b00100010;

            }
        for (j = 1; j <= MAX_LENGTH; j++) {
                up[0][j] = -1000;
                left[0][j] = -1000;
                matrix[0][j] = 0;
                //direction[0][j] = 'D';
                direction[0][j] = 0b000000101;
            } 
        if (TYPE == 'N')
            initialize_NUCC_matrix();
        else if (TYPE == 'P')
            initialize_BLOSUM_matrix();
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
     * Initializes the BLOSUM62 scoring matrix.
     */
    public final void initialize_BLOSUM_matrix() {
        match = new int[256][256];
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

    /**
     * Fills the the similarity and direction matrixes of the two input sequences.
     * First sequence should not be longer than the second sequence.
     * 
     * @param s1 The StringBuilder containing the first sequence
     * @param s2 The StringBuilder containing the second sequence
     */
    public Alignment align(String s1, String s2) {
        int i, j, d;
        int m = s1.length(), n = s2.length();
        seq1.setLength(0); 
        seq1.append(s1);
        seq2.setLength(0); 
        seq2.append(s2);
        if (m < MAX_LENGTH) {
            similarity = Integer.MIN_VALUE;
            if (DEBUG) {
                for (j = 1; j <= n; j++) 
                    System.out.print(String.format("%4c", seq2.charAt(j-1) ));
                System.out.println();
            }
            /*for (i = 1; i <= m; i++) {
                if (DEBUG)
                    System.out.print(seq1.charAt(i-1));
                for (j = 1; j <= n; j++) {
                    up[i][j] = Math.max( up[i-1][j] - GAP_EXT , Math.max(matrix[i-1][j], left[i-1][j]) - GAP_OPEN - GAP_EXT);
                    left[i][j] = Math.max( left[i][j-1] - GAP_EXT , Math.max(matrix[i][j-1], up[i][j-1]) - GAP_OPEN - GAP_EXT);
                    d = match[seq1.charAt(i-1)][seq2.charAt(j-1)] + matrix[i-1][j-1];

                    if (d >= Math.max(up[i][j] , left[i][j])) {
                        matrix[i][j] = d;
                        direction[i][j] = seq1.charAt(i-1) == seq2.charAt(j-1) ? 'M' : 'X';
                    } else if (left[i][j] >= up[i][j]) {
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

                    if (DEBUG) {
                        System.out.print(String.format("%4d", matrix[i][j] ));
                        System.out.print(String.format("%4d", left[i][j] ));
                        System.out.print(String.format("%4d", up[i][j] ));
                        System.out.print(String.format("%4c", direction[i][j] ));
                    }
                }
                if (DEBUG)
                    System.out.println();
            }*/
            for (i = 1; i <= m; i++) {
                if (DEBUG)
                    System.out.print(seq1.charAt(i-1));
                for (j = 1; j <= n; j++) {
                    d = match[seq1.charAt(i-1)][seq2.charAt(j-1)] + matrix[i-1][j-1];

                    if (d >= Math.max(up[i-1][j] , left[i][j-1])) {
                        matrix[i][j] = d;
                        direction[i][j] = 0;
                    } else if (left[i][j-1] >= up[i-1][j]) {
                        matrix[i][j] = left[i][j-1];
                        direction[i][j] = 1;
                    } else {
                        matrix[i][j] = up[i-1][j];
                        direction[i][j] = 2;
                    }
                    if (matrix[i][j] > similarity) {
                        similarity = matrix[i][j];
                        max_i = i;
                        max_j = j;
                    }

                    up[i][j] = Math.max( up[i-1][j] - GAP_EXT , d - GAP_OPEN - GAP_EXT);
                    left[i][j] = Math.max( left[i][j-1] - GAP_EXT , d- GAP_OPEN - GAP_EXT);

                    direction[i][j] |=  ((left[i][j] - GAP_EXT) > (d - GAP_OPEN - GAP_EXT) ? (1<<2) : 0);
                    direction[i][j] |=  ((up[i][j] - GAP_EXT) > (d - GAP_OPEN - GAP_EXT) ? (2<<4) : 0);

                    if (DEBUG) {
                        System.out.print(String.format("%4d", matrix[i][j] ));
                        System.out.print(String.format("%4d", left[i][j] ));
                        System.out.print(String.format("%4d", up[i][j] ));
                        System.out.print(String.format("%4c", direction[i][j] ));
                    }
                }
                if (DEBUG)
                    System.out.println();
            }
        } else {
            System.err.println("Sequences are too large for the aligner. " + m + " " + n + "\n" + seq1.toString() + "\n " + seq2.toString());
            System.exit(0);
        }

        return (new Alignment((short)similarity, -1, max_j, -1, max_i));

        /*align_count++;
        BufferedWriter seq1_writer = null;
        BufferedWriter seq2_writer = null;
        BufferedWriter align_stats_writer = null;

        try {
             seq1_writer = new BufferedWriter(new FileWriter(seq1_file, true));
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }
        try {
            seq2_writer = new BufferedWriter(new FileWriter(seq2_file, true));
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }
        try {
            align_stats_writer = new BufferedWriter(new FileWriter(align_stats_file, true));
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }
        try{
            seq1_writer.append(">" + align_count + "\n");
        }
        catch (IOException e){
            System.err.println("Unable to write in seq 1 file");
            System.exit(1);
        }
        try {
            seq1_writer.append(s1 + "\n");
        }
        catch (IOException e){
            System.err.println("Unable to write in seq 1 file");
            System.exit(1);
        }
        try {
            seq2_writer.append(">" + align_count + "\n");
        }
        catch (IOException e){
            System.err.println("Unable to write in seq 2 file");
            System.exit(1);
        }
        try {
            seq2_writer.append(s2+ "\n");
        }
        catch (IOException e){
            System.err.println("Unable to write in seq 2 file");
            System.exit(1);
        }
       try {
           align_stats_writer.append(">" + align_count + "\n");
       }
       catch (IOException e){
           System.err.println("Unable to write in align stats file");
           System.exit(1);
       }
       try {
           align_stats_writer.append("m: " + m + " n: " + n + "\n");
       }
       catch (IOException e){
           System.err.println("Unable to write in align stats file");
           System.exit(1);
       }
       try {
           align_stats_writer.append("Coordinates = " + max_i + " " + max_j + "\n");
       }
       catch (IOException e){
           System.err.println(e.getMessage());
       }
       try {
           align_stats_writer.append("Alignment: " + this.get_alignment() + "\n");
       }
       catch (IOException e){
           System.err.println("Unable to write in align stats file");
           System.exit(1);
       }
       try {
           align_stats_writer.append("CIGAR: " + this.get_cigar() + "\n");
       }
       catch (IOException e){
           System.err.println("Unable to write in align stats file");
           System.exit(1);
       }
       try {
           align_stats_writer.append("offset: " + this.get_offset() + "\n");
       }
       catch (IOException e){
           System.err.println("Unable to write in align stats file");
           System.exit(1);
       }
       try {
           align_stats_writer.append("Similarity: " + similarity + "\n");
       }
       catch (IOException e){
           System.err.println("Unable to write in align stats file");
           System.exit(1);
       }

       try {
           seq1_writer.close();
       }
       catch (IOException e){
           System.err.println(e.getMessage());
       }
        try {
            seq2_writer.close();
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }
        try {
            align_stats_writer.close();
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }




        if (DEBUG) {
            System.out.println(s2);
            System.out.println(s1);
            System.out.println("m: " + m + " n: " + n);
            System.out.println("Coordinates = "+ max_i + " " + max_j);
            System.out.println(this.get_alignment());
            System.out.println(this.get_cigar());
            System.out.println("offset: "+ this.get_offset());
            System.out.println(similarity + "\n");
        }*/

    }    
    
    /**
     * Calculates the alignment from the similarity matrix. 
     * Call align() before this function.
     * 
     * @return the alignment string which may contains some gaps. 
     */
    public String get_alignment() {
        int i, j;
        int range[];
        StringBuilder subject = new StringBuilder();
        StringBuilder query = new StringBuilder();
        subject.setLength(0);
        query.setLength(0);
        i = max_i;
        j = max_j;
        if (CLIPPING_STRIGENCY > 0)
            range = calculate_clip_range();
        else
            range = new int[]{1, max_i, 1, max_j};
        while (i > 0 && j > 0) {
            if (CLIPPING_STRIGENCY > 0 && i < range[0])
                break;
            if (direction[i][j] == 'I') {
                query.append( seq1.charAt(i-1) );
                subject.append( '-' );
                i = i - 1;
            } else if (direction[i][j] == 'D') {
                query.append( '-' );
                subject.append( seq2.charAt(j-1) );
                j = j - 1;
            } else {
                query.append( seq1.charAt(i-1) );
                subject.append( seq2.charAt(j-1) );
                i = i - 1;
                j = j - 1;
            }
        } 
        if (CLIPPING_STRIGENCY > 0) {
            for (;i > 0 && j > 0; --i, --j) {
                query.append( seq1.charAt(i-1) );
                subject.append( seq2.charAt(j-1) );
            }
        }
        for (;i > 0; --i) {
            query.append( seq1.charAt(i-1) );
            subject.append( '-' );
        }
        for (;j > 0; --j) {
            query.append( '-' );
            subject.append( seq2.charAt(j-1) );
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
    
   // System.out.println("original_score2 " + score + " " + score + " " + (score * 100.0 / p_score));System.out.println("original_score2 " + score + " " + score + " " + (score * 100.0 / p_score));

    /**
     * @return The similarity score of two sequences, as the largest entry in 
     *         similarity matrix.
     */
    public int get_similarity(){
        return similarity;
    }
    
    /**
     * @return The identity of two sequences after being aligned.
     */
    public double get_identity(){
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
        int[] range = new int[4];
        x = i = max_i;
        j = max_j;
        while (i > 0 && j > 0) {
            if (direction[i][j] == 'I') {
                score_array[x--] = insertion_penalty;
                i = i - 1;
            } else if (direction[i][j] == 'D') {
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
        for (i = 2; i <= max_i; ++i){
            if (score_array[i] > max_ending_here + score_array[i]){
                max_ending_here = score_array[i];
                tmp_start = tmp_stop = i;
            } else {
                max_ending_here = max_ending_here + score_array[i];
                tmp_stop = i;
            }
            if (max_so_far < max_ending_here){
                range[0] = tmp_start;
                range[1] = tmp_stop;
                max_so_far = max_ending_here;
            }
        }
        i = max_i;
        j = max_j;
        while (i > 0 && j > 0) {
            if (i == range[1]){
                range[3] = j;
            }
            if (i == range[0]){
                range[2] = j;
            }
            if (direction[i][j] == 'I') {
                i = i - 1;
            } else if (direction[i][j] == 'D') {
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
    public String get_cigar(boolean cigar_with_x) {
        int i, j, move_counts, count, identicals = 0;
        int range[];
        char curr_move, prev_move, operation;
        insertions = deletions = 0;
        operation_stack.clear();
        count_stack.clear();
        cigar.setLength(0);
        if (CLIPPING_STRIGENCY > 0){
            range = calculate_clip_range();
            if (seq1.length() - range[1] > 0){
                operation_stack.push('S');
                count_stack.push(seq1.length() - range[1]);
            }
            prev_move = 'M';
            move_counts = 1;
        } else {
            range = new int[]{1, max_i, 1, max_j};
            prev_move = 'M';
            move_counts = 1;//range[1] < seq1.length() ? seq1.length() - range[1] + 1 : 1;
        }
        range_len = range[1] - range[0] + 1;
        i = seq1.length();
        j = range[3] + seq1.length() - range[1];
        if (j > seq2.length()){
            i -= j - seq2.length();
            j -= j - seq2.length();
        }
        for (; i >= range[1]; --i, --j)
             if (seq1.charAt(i-1) == seq2.charAt(j-1))
                identicals++;
        /*i = range[1] - 1;
        j = range[3] - 1;
        while (i >= range[0]){
            curr_move = direction[i][j];
            if (curr_move == 'I'){
                i = i - 1;
                ++insertions;
            } else if (curr_move == 'D'){
                j = j - 1;
                ++deletions;
            } else {
                if (seq1.charAt(i) == seq2.charAt(j))
                    ++identicals;
                if(!cigar_with_x){
                    curr_move = 'M';
                }
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
        }*/
        i = range[1];
        j = range[3];
        int which = 0;
        int first_time = 1;
        while (i >= range[0]){
            //curr_move = direction[i][j];

           which = (direction[i][j] >> (which<<1)) & 3;
            //System.out.println(String.format("direction:%d, which:%d, i:%d, j:%d", direction[i][j], which, i, j));
            //System.out.flush();
            if (which == 2){
                i = i - 1;
                ++insertions;
                curr_move = 'I';
            } else if (which == 1){
                j = j - 1;
                ++deletions;
                curr_move = 'D';
            } else {
                if (cigar_with_x) {
                    if (seq1.charAt(i - 1) == seq2.charAt(j - 1)) {
                        curr_move = '=';
                        ++identicals;
                    } else {
                        curr_move = 'X';
                    }
                } else {
                    if (seq1.charAt(i - 1) == seq2.charAt(j - 1))
                        ++identicals;
                    curr_move = 'M';
                }
                i = i - 1;
                j = j - 1;
            }

            if (first_time == 1){
                prev_move = curr_move;
                first_time = 0;
            } else {
                if (prev_move == curr_move)
                    ++move_counts;

                else {
                    operation_stack.push(prev_move);
                    count_stack.push(move_counts);
                    move_counts = 1;
                }
                prev_move = curr_move;
            }

        }
        offset = j - i;
        if (CLIPPING_STRIGENCY > 0){
            operation_stack.push(prev_move);
            count_stack.push(move_counts);
            if (i > 0){
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
          while (!operation_stack.isEmpty()){
            operation = operation_stack.pop();
            count = count_stack.pop();
            cigar.append(count).append(operation);
        }
        for (; i > 0 && j > 0; --i, --j)
             if (seq1.charAt(i) == seq2.charAt(j))
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
    public int get_offset(){
        return offset;
    }

    /**
     * @return The length of the alignment after soft-clipping
     */
    public int get_range_length(){
        return range_len;
    }
    
    /**
     * @return The total number of insertions in the alignment
     */
    public int get_insertions(){
        return insertions;
    }
    
    /**
     * @return The total number of deletions in the alignment
     */
    public int get_deletions(){
        return deletions;
    }

    /**
     * Do Striped Smith-Waterman alignment. Warning: No parameter checking is performed. Incorrect arguments are likely to crash the JVM.
     *
     * @param read pointer to the query sequence; the query sequence needs to be numbers
     * @param flattenedMatrix pointer to the substitution matrix; mat needs to be corresponding to the read sequence
     * @param n the square root of the number of elements in mat (mat has n*n elements)
     * @param score_size estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
    your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
     * @param ref pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of function initprofile
     * @param gapOpen the absolute value of gap open penalty.
     * @param gapExtend the absolute value of gap extension penalty.
     * @param flag bitwise FLAG; (from high to low) bit 5: when setted as 1, the function will return the best alignment
    beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
    < filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
    cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
    will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
    setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
    the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
    Note: Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
    and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
    while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
    0-based coordinate.
     * @param filterscore score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
    decription of the flag parameter for detailed usage.)
     * @param filterdistance distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
    the decription of the flag parameter for detailed usage.)
     * @param maskLen The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
    readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
    return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
    alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
    largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
    picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
    reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
    score from the unmasked elements.
     * @return Smith-Waterman alignment
     */
    public static native Alignment alignNative(byte[] read, byte[] flattenedMatch, int n, byte[] ref, int gapOpen, int gapExtend);
    /**
     * Performs striped Smith-Waterman alignment
     *
     * @param read read sequence
     * @param ref reference sequence
     * @return Smith-Waterman alignment
     */
    public static Alignment ssw_align(byte[] read, byte[] ref) {
        //if (flag != flag & 0xFF) throw new IllegalArgumentException("Only lowest 8 bits of flag are meaningful");
        int[] lookup = new int[257]; // lookup[256] is used as sentinal for number of unique bases/matrix size
        java.util.Arrays.fill(lookup, -1);
        lookup[256] = 0;
        byte[] readNum = convertToNumeric(lookup, read);
        byte[] refNum = convertToNumeric(lookup, ref);
        byte[] flattenedMatch = flatten(lookup, match);
        int uniqueBases = lookup[256];
        assert(flattenedMatch.length == uniqueBases * uniqueBases);
        assert(maxValue(readNum) < uniqueBases);
        assert(maxValue(refNum) < uniqueBases);
        Alignment alignment = alignNative(
                readNum, flattenedMatch, uniqueBases,
                refNum, GAP_OPEN, GAP_EXT);
        return alignment;
    }
    private static int maxValue(byte[] array) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < array.length; i++) {
            if (array[i] > max) {
                max = array[i];
            }
        }
        return max;
    }
    /**
     * Converts an ASCII sequence into numeric successive 0-based values
     * @param lookup ASCII to numeric conversion lookup array
     * @param sequence ASCII sequence
     * @return numeric sequence
     */
    private static byte[] convertToNumeric(int[] lookup, byte[] sequence) {
        byte[] numericSeq = new byte[sequence.length];
        for (int i = 0; i < sequence.length; i++) {
            int b = sequence[i];
            if (lookup[b] == -1) {
                lookup[b] = lookup[256]++;
            }
            numericSeq[i] = (byte)lookup[b];
        }
        return numericSeq;
    }
    /**
     * Generates a flattened ssw scoring matrix
     * @param lookup ASCII to numeric conversion lookup array
     * @param matrix scoring matrix
     * @return flattened ssw numeric scoring matrix
     */
    private static byte[] flatten(int[] lookup, int[][] matrix) {
        int size = lookup[256];
        byte[] flattened = new byte[size * size];
        for (int i = 0; i < matrix.length; i++) {
            int newi = lookup[i];
            if (newi == -1) continue;
            for (int j = 0; j < matrix[i].length; j++) {
                int newj = lookup[j];
                if (newj == -1) continue;
                int score = matrix[i][j];
                if (score < Byte.MIN_VALUE || score > Byte.MAX_VALUE) {
                    throw new IllegalArgumentException("Scoring matrix values must fit into signed 8-bit integer");
                }
                flattened[newi * size + newj] = (byte)score;
            }
        }
        return flattened;
    }

}
