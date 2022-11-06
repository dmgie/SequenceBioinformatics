package assignment03;

import assignment03.FastA_YOUR_NAME;

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
			System.out.println(objectiveFunction(list));

			// 2. write out all the simple mixed cycle constraints between any two sequences

			// 3. write out all the simple mixed cycle constraints between any three sequences

			// 4. write out the binary variable constraints

			// 5. specify all variables as integers
		}
	}

	public String objectiveFunction(FastA_YOUR_NAME.Pair seq1, FastA_YOUR_NAME.Pair seq2, FastA_YOUR_NAME.Pair seq3){
		String fun = "";
		String s0 = seq1.sequence();
		String s1 = seq2.sequence();
		String s2 = seq3.sequence();
		for(int i = 0; i < s0.length(); i++){
			for(int j = 0; j < s1.length();j++){
				fun += "+" +String.valueOf(score(s0.charAt(i), s1.charAt(j)))+"*x0"+ String.valueOf(i)+ "_1"+ String.valueOf(j);
			}
		}
		for(int i = 0; i < s0.length(); i++){
			for(int j = 0; j < s2.length();j++){
				fun += "+" +String.valueOf(score(s0.charAt(i), s2.charAt(j)))+"*x0"+ String.valueOf(i)+ "_2"+ String.valueOf(j);
			}
		}
		for(int i = 0; i < s1.length(); i++){
			for(int j = 0; j < s2.length();j++){
				fun += "+" +String.valueOf(score(s1.charAt(i), s2.charAt(j)))+"*x1"+ String.valueOf(i)+ "_2"+ String.valueOf(j);
			}
		}
		return fun;
	}
	 public int score(char x, char y){
		if(x == y){
			return 4;
		}
		else{
			return 1;
		}
	 }

}
