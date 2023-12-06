package nl.escience.alignment;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Stack;

import java.io.File;

import java.io.FileReader;

public class Example {


    public static int MAX_ALIGNMENT_LENGTH = 2000;
    public static int GAP_OPEN = 20;
    public static int GAP_EXT = 3;
    public static int ALIGNMENT_BOUND = 5; // This value is between fast and very fast
    public static int CLIPPING_STRINGENCY = 0; // 0: no-clipping // This matches the value of fast
                                                // 1: low
                                                // 2: medium
                                                // 3: high



    public static void main(String[] args) {
        try {
            System.loadLibrary("sswjni");
        } catch (java.lang.UnsatisfiedLinkError e) {
            System.out.println(String.format("Cannot find libsswjni.so. Has the library been built and LD_LIBRARY_PATH or -Djava.library.path set appropriately?\n%s", e));
            throw e;
        }

        String seq1, seq2, seq1_no=null, seq2_no=null;
        boolean cigar_with_x = true;
        BufferedReader seq1_reader = null;
        BufferedReader seq2_reader = null;
        BoundedLocalSequenceAlignment bounded_aligner = new BoundedLocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, ALIGNMENT_BOUND, CLIPPING_STRINGENCY, 'N');
        LocalSequenceAlignment aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIPPING_STRINGENCY, 'N');
        try {
            seq1_reader = new BufferedReader(new FileReader(args[0]));
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }
        try{
            seq2_reader = new BufferedReader(new FileReader(args[1]));
         }
        catch (IOException e){
        System.err.println(e.getMessage());
        }

        long totaltime = 0;
        try {
            while ((seq1 = seq1_reader.readLine()) != null && (seq2 = seq2_reader.readLine()) != null) {
                if (seq1.charAt(0) == '>' && seq2.charAt(0) != '>') {
                    System.out.println("Input files not formatted correctly");
                    System.exit(1);
                }
                if (seq1.charAt(0) == '>' && seq2.charAt(0) == '>') {
                    seq1_no = seq1.substring(1);
                    seq2_no = seq2.substring(1);
                    continue;
                }
                seq1.toUpperCase();
                seq2.toUpperCase();

                long startime = System.currentTimeMillis();
                Alignment aln = aligner.ssw_align(seq1.getBytes(), seq2.getBytes());
                Alignment aln = aligner.align(seq1, seq2);
                //Alignment aln = bounded_aligner.align(seq1, seq2+seq2);
                totaltime += System.currentTimeMillis() - startime;

                //
                if (aln == null) {
                    throw new RuntimeException();
                }
                System.out.println(">"+seq1_no+ "|"+seq2_no);
                System.out.println(aln.toString() + String.format(",CIGAR=%s",aligner.get_cigar(cigar_with_x)));

            }

            System.err.println(totaltime);
            seq1_reader.close();
            seq2_reader.close();
        }
        catch (IOException e){
            System.err.println(e.getMessage());
        }



        }




}