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
		// Read file using FastA read
		var contents = FastA_YOUR_NAME.read(args[0]);
		if(args.length == 2){
			FastA_YOUR_NAME.write(contents, args[1]);
		}
		else{
			FastA_YOUR_NAME.write(contents, null);
		}


	}
}
