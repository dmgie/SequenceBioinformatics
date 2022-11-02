package assignment03;

import assignment01.FastA_YOUR_NAME;

import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;

/**
 * setup ILP to solve multiple sequence alignment for three sequences
 * Sequence Bioinformatics, WS22/23
 */
public class AlignmentILP_YOUR_NAME {
	public static void main(String[] args) throws IOException {
		if (args.length != 1 && args.length!=2) {
			throw new IOException("Usage: AlignmentILP_YOUR_NAME input [output]");
		}

		var list = FastA_YOUR_NAME.read(args[0]);

		if (list.size() != 3) {
			throw new IOException("Input file must contain 3 sequences, found: " + list.size());
		}

		// todo: setup and write out ILP based on extended alignment graph, with match score =  4 and mismatch score = 1

		try(var w=(args.length ==2?new FileWriter(args[1]):new OutputStreamWriter(System.out))) {
			w.write("max: ");
			// 1. write the objective function: loop over all pairs of sequences and all pairs of letters

			// 2. write out all the simple mixed cycle constraints between any two sequences

			// 3. write out all the simple mixed cycle constraints between any three sequences

			// 4. write out the binary variable constraints

			// 5. specify all variables as integers
		}
	}

}
