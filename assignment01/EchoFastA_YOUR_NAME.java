package assignment01;

import java.io.IOException;

/**
 * EchoFastA_YOUR_NAME Author(s): YOUR_NAME Sequence Bioinformatics, WS 22/23
 */
public class EchoFastA_YOUR_NAME {
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2)
			throw new IOException("Usage: EchoFastA_YOUR_NAME infile [outFile]");

		// todo: read FastA records from infile and echo to outfile or stdout (console)
		var filename = args[0];
		// Read file using buffered reader
		var reader = new java.io.BufferedReader(new java.io.FileReader(filename));
		while (reader.ready()) {
			var line = reader.readLine();
			System.out.println(line);
		}
	}
}
