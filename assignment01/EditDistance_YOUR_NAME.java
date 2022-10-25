package assignment01;

import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

/**
 * EditDistance_YOUR_NAME Author(s): YOUR_NAME Sequence Bioinformatics, WS 22/23
 */
public class EditDistance_YOUR_NAME {

    public static void main(String[] args) throws IOException {
        if (args.length < 1 || args.length > 2)
            throw new IOException("Usage: EditDistance_YOUR_NAME infile [outFile]");

        // todo: implement input of FastA records

        // todo: check that all input sequences have the same length, otherwise throw a
        // new IOException("Different lengths");
        var seqs = FastA_YOUR_NAME.read(args[0]);

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
        System.out.println(pairwise);

        try (Writer w = (args.length == 2 ? new FileWriter(args[1]) : new OutputStreamWriter(System.out))) {
            // todo: compute distance between any two sequences, using method
            // computeEditDistance(x,y) defined below
            // todo: write distance matrix

            for (var pair : pairwise) {
                var dist = computeEditDistance(pair.seq1(), pair.seq2());
                w.write(dist + " ");
            }
        }
        // example of format:
        // a 0 1 2 3
        // b 1 0 4 5
        // c 2 4 0 6
        // d 3 5 6 0
    }

    // Function to get all pariwise combinations from a list
    // Maybe just take int a List<FastA_YOUR_NAME.Pair> and return a
    // ArrayList<Tuple>, instead of first needing a list
    public static ArrayList<Tuple> generate_combinations(List<FastA_YOUR_NAME.Pair> sequences) {
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
