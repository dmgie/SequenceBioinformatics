package assignment01;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import assignment01.*;

/**
 * Translate_YOUR_NAME Author(s): YOUR_NAME Sequence Bioinformatics, WS 22/23
 */
public class Translate_YOUR_NAME {
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2)
			throw new IOException("Usage: Translate_YOUR_NAME infile [outFile]");

		// todo: read in FastA pairs
		var list = FastA_YOUR_NAME.read(args[0]);

		// todo: compute translated sequences using translate(sequence) method defined
		// below
		var translated = new ArrayList<FastA_YOUR_NAME.Pair>();

		for (var pair : list) {
			translated.add(new FastA_YOUR_NAME.Pair(pair.header(), translate(pair.sequence())));
		}
		// todo: write translated sequences
		if(args.length == 2){
			FastA_YOUR_NAME.write(translated, args[1]);
		}
		else{
			FastA_YOUR_NAME.write(translated, null);
		}
	}

	public static String translate(String sequence) {
		var buf = new StringBuilder();

		// todo: implement translation of sequence
		var mapping = new HashMap<String, String>();

		// Hashmap of codons to amino acid
		mapping.put("TTT", "F");
		mapping.put("TTC", "F");
		mapping.put("TTA", "L");
		mapping.put("TTG", "L");
		mapping.put("TCT", "S");
		mapping.put("TCC", "S");
		mapping.put("TCA", "S");
		mapping.put("TCG", "S");
		mapping.put("TAT", "Y");
		mapping.put("TAC", "Y");
		mapping.put("TAA", "X");
		mapping.put("TAG", "X");
		mapping.put("TGT", "C");
		mapping.put("TGC", "C");
		mapping.put("TGA", "X");
		mapping.put("TGG", "W");
		mapping.put("CTT", "L");
		mapping.put("CTC", "L");
		mapping.put("CTA", "L");
		mapping.put("CTG", "L");
		mapping.put("CCT", "P");
		mapping.put("CCC", "P");
		mapping.put("CCA", "P");
		mapping.put("CCG", "P");
		mapping.put("CAT", "H");
		mapping.put("CAC", "H");
		mapping.put("CAA", "Q");
		mapping.put("CAG", "Q");
		mapping.put("CGT", "R");
		mapping.put("CGC", "R");
		mapping.put("CGA", "R");
		mapping.put("CGG", "R");
		mapping.put("ATT", "I");
		mapping.put("ATC", "I");
		mapping.put("ATA", "I");
		mapping.put("ATG", "M");
		mapping.put("ACT", "T");
		mapping.put("ACC", "T");
		mapping.put("ACA", "T");
		mapping.put("ACG", "T");
		mapping.put("AAT", "N");
		mapping.put("AAC", "N");
		mapping.put("AAA", "K");
		mapping.put("AAG", "K");
		mapping.put("AGT", "S");
		mapping.put("AGC", "S");
		mapping.put("AGA", "R");
		mapping.put("AGG", "R");
		mapping.put("GTT", "V");
		mapping.put("GTC", "V");
		mapping.put("GTA", "V");
		mapping.put("GTG", "V");
		mapping.put("GCT", "A");
		mapping.put("GCC", "A");
		mapping.put("GCA", "A");
		mapping.put("GCG", "A");
		mapping.put("GAT", "D");
		mapping.put("GAC", "D");
		mapping.put("GAA", "E");
		mapping.put("GAG", "E");
		mapping.put("GGT", "G");
		mapping.put("GGC", "G");
		mapping.put("GGA", "G");
		mapping.put("GGG", "G");

		// go through sequence in steps of 3
		for (int i = 0; i < sequence.length(); i += 3) {
			if (i + 3 > sequence.length()) {
				break;
			}
			// get codon
			var codon = sequence.substring(i, i + 3);

			// get amino acid
			var aminoAcid = mapping.get(codon);
			// append amino acid to buffer
			buf.append(aminoAcid);
		}

		return buf.toString();
	}
}
