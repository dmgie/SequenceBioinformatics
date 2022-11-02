package assignment01;

import java.io.IOException;

/**
 * EchoFastA_GIESEL_MUEHLBAUER Author(s): GIESEL_MUEHLBAUER Sequence Bioinformatics, WS 22/23
 */
public class EchoFastA_GIESEL_MUEHLBAUER {
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2)
			throw new IOException("Usage: EchoFastA_GIESEL_MUEHLBAUER infile [outFile]");

		// todo: read FastA records from infile and echo to outfile or stdout (console)
		// Read file using FastA read
		var contents = FastA_GIESEL_MUEHLBAUER.read(args[0]);
		if(args.length == 2){
			FastA_GIESEL_MUEHLBAUER.write(contents, args[1]);
		}
		else{
			FastA_GIESEL_MUEHLBAUER.write(contents, null);
		}


	}
}
