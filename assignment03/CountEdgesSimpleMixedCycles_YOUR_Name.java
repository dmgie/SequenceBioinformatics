package assignment03;

import java.io.IOException;

/**
 * count number of edges between nucleotides in different sequences and the
 * number of mixed cycles Sequence Bioinformatics, WS22/23
 */
public class CountEdgesSimpleMixedCycles_YOUR_Name {
	public static void main(String[] args) throws IOException {
		if (args.length != 3)
			throw new IOException("Usage: CountEdgesSimpleMixedCycles_YOUR_Name aLength bLength cLength");

		var length = new int[] { Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]) };

		System.out.println(CountEdgesSimpleMixedCycles_YOUR_Name.class.getSimpleName());

		System.out.printf("Sequence lengths: %d %d %d%n", length[0], length[1], length[2]);

		// todo: report the number of edges between nucleotides in different sequences

		var numEdgesBetweenDifferenteSequences = count_edges(length);

		System.out.printf("Edges between different sequences: %d%n", numEdgesBetweenDifferenteSequences);

		var numSimpleMixedCycles = 0;

		// todo: implement counting of number of simple mixed cycles

		// first compute the number of simple mixed cycles that use two cycles

		// then compute and add the number of simple mixed cycles that use three cycles

		System.out.printf("Total simple mixed cycles: %d%n", numSimpleMixedCycles);
	}

	public static int count_edges(int[] numlist) {
		// Multiply each length with each length and sum them up
		int sum = 0;
		for (int i = 0; i < numlist.length; i++) {
			for (int j = 0; j < numlist.length; j++) {
				sum += numlist[i] * numlist[j];
			}
		}
		return sum;
	}

	public static void get_mixed_cycles(int[] numlist) {
		// calculate the number of simple mixed cycles given sequence lengths

	}
}
