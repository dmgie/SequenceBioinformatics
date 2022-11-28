package assignment06;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;

/**
 * FastA_GIESEL_MUEHLBAUER Author(s): GIESEL_MUEHLBAUER Sequence Bioinformatics, WS 22/23
 */
public class FastA_GIESEL_MUEHLBAUER {

	// Main function
	public static void main(String[] args) throws IOException {
		// Read the input file
		ArrayList<Pair> input = read("temp.fasta");

		// Print the input
		System.out.println("Input:");
		for (Pair p : input) {
			System.out.println(p);
		}

		write(input, "tempfile.fasta");

		// Split the input file into the header and the sequence
	}

	public static void write(Collection<Pair> list, String fileName) throws IOException {
		try (var w = (fileName != null ? new FileWriter(fileName) : new OutputStreamWriter(System.out))) {
			// write each pair to a file, with a maximum width of 80
			for (Pair p : list) {
				w.write(">" + p.header + "\n" + p.sequence.replaceAll(".{80}(?=.)", "$0\n") + "\n");
			}
		}
	}

	public static ArrayList<Pair> read(String fileName) throws IOException {
		var list = new ArrayList<Pair>();
		try (var r = new BufferedReader(new FileReader(fileName))) {
			// NOTE: Don't know a better way to do it, I'm reading it into the file cause
			// all other ways cause the second sequence not to be included, if I read the
			// file line by line

			// read the entire file until the end into a string
			String all;
			{
				var sb = new StringBuilder();
				String line;
				while ((line = r.readLine()) != null) {
					sb.append(line);
					sb.append('\n');
				}
				all = sb.toString();
			}

			// split at '>'
			String[] parts = all.split(">");

			// for each split part (header+sequence combo), split at first occurence of '\n'
			for (String part : parts) {
				if (part.isEmpty()) {
					continue;
				}
				String[] headerAndSequence = part.split("\n", 2);
				list.add(new Pair(headerAndSequence[0], headerAndSequence[1].replace("\n", "")));
			}

		}
		return list;
	}

	/**
	 * a FastA record consisting of a pair of header and sequence
	 */
	public static record Pair(String header, String sequence) {
	}
}
