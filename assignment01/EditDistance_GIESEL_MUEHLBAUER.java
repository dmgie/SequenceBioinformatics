package assignment01;

import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * EditDistance_GIESEL_MUEHLBAUER Author(s): GIESEL_MUEHLBAUER Sequence Bioinformatics, WS 22/23
 */
public class EditDistance_GIESEL_MUEHLBAUER {

    public static void main(String[] args) throws IOException {
        if (args.length < 1 || args.length > 2)
            throw new IOException("Usage: EditDistance_GIESEL_MUEHLBAUER infile [outFile]");

        // todo: implement input of FastA records

        // todo: check that all input sequences have the same length, otherwise throw a
        // new IOException("Different lengths");
        var seqs = FastA_GIESEL_MUEHLBAUER.read(args[0]);

        var length = 0;
        for (var pair : seqs) {
            // check if they all have the same length
            if (length == 0)
                length = pair.sequence().length();
            else if (length != pair.sequence().length())
                throw new IOException("Different lengths");
        }

        // var only_seq = get_sequences_as_list(seqs);
        var pairwise = generate_combinations(seqs);


        try (Writer w = (args.length == 2 ? new FileWriter(args[1]) : new OutputStreamWriter(System.out))) {
            // todo: compute distance between any two sequences, using method
            // computeEditDistance(x,y) defined below
            // todo: write distance matrix
            String matrix [] [] = new String [seqs.size()+1] [seqs.size()+1] ;
            for (int i = 1; i <= seqs.size(); i++){
                matrix [i][i] = " 0 ";
                matrix [0][i] = "se"+ String.valueOf(i);
                matrix [i][0] = "se"+ String.valueOf(i)+" ";
            }
            int c = 1;
            int r = 1;
            int max = seqs.size();
            for (var pair : pairwise) {
                var dist = computeEditDistance(pair.seq1(), pair.seq2());
                c++;
                if (c > max){
                    r++;
                    c = r+1;
                }
                matrix[r][c] = String.valueOf(dist);
                matrix[c][r] = String.valueOf(dist);

            }

            for(String[] s : matrix) {
                w.write(Arrays.toString(s)+"\n");
            }
        }
        // example of format:
        // a 0 1 2 3
        // b 1 0 4 5
        // c 2 4 0 6
        // d 3 5 6 0
    }

    // Function to get all pariwise combinations from a list
    // Maybe just take int a List<FastA_GIESEL_MUEHLBAUER.Pair> and return a
    // ArrayList<Tuple>, instead of first needing a list
    public static ArrayList<Tuple> generate_combinations(List<FastA_GIESEL_MUEHLBAUER.Pair> sequences) {
        var combinations = new ArrayList<Tuple>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i + 1; j < sequences.size(); j++) {
                combinations.add(new Tuple(sequences.get(i).sequence(), sequences.get(j).sequence()));
            }
        }
        return combinations;

    }

    private static int computeEditDistance(String x, String y) {
        // todo: implement computation of edit distance
        // simple implementation
        int distance = 0;
        for (int i = 0; i < x.length(); i++) {
            if (x.charAt(i) != y.charAt(i))
                distance++;
        }
        return distance;
    }

    public static record Tuple(String seq1, String seq2) {
    }

}
