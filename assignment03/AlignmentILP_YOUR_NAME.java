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
			// 1. write the objective function: loop over all pairs of sequences and all pairs of letters
			String obj_fun = objectiveFunction(list.get(0),list.get(1),list.get(2));

			// 2. write out all the simple mixed cycle constraints between any two sequences
			// 3. write out all the simple mixed cycle constraints between any three sequences
			String mix_cyc = mixedCycles(list.get(0),list.get(1),list.get(2));

			// 4. write out the binary variable constraints
			String bin = binary(list.get(0),list.get(1),list.get(2));

			// 5. specify all variables as integers
			String in = makeInt(list.get(0),list.get(1),list.get(2));
			w.write("max: " + obj_fun + mix_cyc + in);
		}
	}

	public static String objectiveFunction(FastA_YOUR_NAME.Pair seq1, FastA_YOUR_NAME.Pair seq2, FastA_YOUR_NAME.Pair seq3){
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
		return fun+";";
	}
	public static int score(char x, char y){
		if(x == y){
			return 4;
		}
		else{
			return 1;
		}
	}

	public static String mixedCycles(FastA_YOUR_NAME.Pair seq1, FastA_YOUR_NAME.Pair seq2, FastA_YOUR_NAME.Pair seq3){
		String mix_cy = "";
		int s0 = seq1.sequence().length();
		int s1 = seq2.sequence().length();
		int s2 = seq3.sequence().length();

		for (int i = 0; i < s0; i++) {
			for (int j = i; j < s0; j++) {
				for (int k = 0; k < s2; k++) {
					for (int l = k; l < s2; l++) {
						if (i != j || k != l) {
							mix_cy += "x0" + String.valueOf(j) + "_2" + String.valueOf(k) + " + x0" + String.valueOf(i) + "_2" + String.valueOf(l) + "< 1;\n";
						}
					}
				}
			}
		}

		for (int i = 0; i < s1; i++) {
			for (int j = i; j < s1; j++) {
				for (int k = 0; k < s2; k++) {
					for (int l = k; l < s2; l++) {
						if (i != j || k != l) {
							mix_cy += "x1" + String.valueOf(j) + "_2" + String.valueOf(k) + " + x1" + String.valueOf(i) + "_2" + String.valueOf(l) + "< 1;\n";
						}
					}
				}
			}
		}

		for (int i = 0; i < s0; i++) {
			for (int j = i; j < s0; j++) {
				for (int k = 0; k < s1; k++) {
					for (int l = k; l < s1; l++) {
						if (i != j || k != l) {
							mix_cy += "x0" + String.valueOf(j) + "_1" + String.valueOf(k) + " + x0" + String.valueOf(i) + "_1" + String.valueOf(l) + "< 1;\n";
						}
						for (int n = 0; n < s2; n++) {
							for (int m = n; m < s2; m++) {
								if (i != j || k != l || m != n) {
									mix_cy += "x0" + String.valueOf(j) + "_1" + String.valueOf(k) + " + x1" + String.valueOf(l) + "_2" + String.valueOf(m) +
											" + x0" + String.valueOf(i) + "_2" + String.valueOf(n) + "< 2;\n";
									mix_cy += "x0" + String.valueOf(j) + "_2" + String.valueOf(m) + " + x1" + String.valueOf(k) + "_2" + String.valueOf(n) +
											" + x0" + String.valueOf(i) + "_1" + String.valueOf(l) + "< 2;\n";
								}
							}
						}
					}
				}
			}
		}
		return mix_cy;
	}
	public static String binary(FastA_YOUR_NAME.Pair seq1, FastA_YOUR_NAME.Pair seq2, FastA_YOUR_NAME.Pair seq3){
		String con = "";
		int s0 = seq1.sequence().length();
		int s1 = seq2.sequence().length();
		int s2 = seq3.sequence().length();
		for(int i = 0; i < s0; i++){
			for(int j = 0; j < s1;j++){
				con += "x0" + String.valueOf(i)+ "_1" + String.valueOf(j) + "<1;\n";
			}
		}
		for(int i = 0; i < s0; i++){
			for(int j = 0; j < s2;j++){
				con += "x0" + String.valueOf(i)+ "_2" + String.valueOf(j) + "<1;\n";
			}
		}
		for(int i = 0; i < s1; i++){
			for(int j = 0; j < s2;j++){
				con += "x1" + String.valueOf(i)+ "_2" + String.valueOf(j) + "<1;\n";
			}
		}
		return con;
	}

	public static String makeInt(FastA_YOUR_NAME.Pair seq1, FastA_YOUR_NAME.Pair seq2, FastA_YOUR_NAME.Pair seq3){
		String int_con = "int";
		int s0 = seq1.sequence().length();
		int s1 = seq2.sequence().length();
		int s2 = seq3.sequence().length();
		for(int i = 0; i < s0; i++){
			for(int j = 0; j < s1;j++){
				int_con += " x0" + String.valueOf(i)+ "_1" + String.valueOf(j) + ",";
			}
		}
		for(int i = 0; i < s0; i++){
			for(int j = 0; j < s2;j++){
				int_con += " x0" + String.valueOf(i)+ "_2" + String.valueOf(j) + ",";
			}
		}
		for(int i = 0; i < s1; i++){
			for(int j = 0; j < s2;j++){
				int_con += " x1" + String.valueOf(i)+ "_2" + String.valueOf(j) + ",";
			}
		}
		return int_con.substring(0, int_con.length()-1)+ ";\n";
	}

}
