package assignment02;

import assignment02.FastA_GIESEL_MUEHLBAUER;

import javax.sql.rowset.spi.SyncResolver;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.RandomAccess;

// TODO: Implement as row-major order vector instead of 2D array

/**
 * GlobalAligner_GIESEL_MUEHLBAUER Sequence Bioinformatics, WS 22/23
 */
public class GlobalAligner_GIESEL_MUEHLBAUER {
	int gap_penalty = 1;
	int match_score = 1;
	int mismatch_score = -1;

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2)
			throw new IOException("Usage: GlobalAligner_GIESEL_MUEHLBAUER infile [quadraticSpace|linearSpace|noDP]");

		var list = FastA_GIESEL_MUEHLBAUER.read(args[0]);

		if (list.size() != 2)
			throw new IOException("Wrong number of input sequences: " + list.size());

		var mode = (args.length == 2 ? args[1] : "quadraticSpace");

		switch (mode) {
		case "quadraticSpace" -> runNeedlemanWunschQuadraticSpace(list.get(0), list.get(1));
		case "linearSpace" -> runNeedlemanWunschLinearSpace(list.get(0), list.get(1));
		case "noDP" -> runNeedlemanWunschRecursively(list.get(0), list.get(1));
		default -> throw new IOException("Unknown mode: " + mode);
		}
	}

	/**
	 * computes the optimal global alignment score and an alignment, using quadratic
	 * space. Prints out the optimal score and a corresponding alignment Also prints
	 * out the number of milliseconds required for the computation
	 *
	 * @param x
	 * @param y
	 */
	public static void runNeedlemanWunschQuadraticSpace(FastA_GIESEL_MUEHLBAUER.Pair x, FastA_GIESEL_MUEHLBAUER.Pair y) {
		// todo: implement, Assignment 2.1
		int gap_penalty = 1;
		int match_score = 1;
		int mismatch_score = -1;

		long start = System.currentTimeMillis();
		// Create matrix with first row and column being 0 row and 0 column
		int[][] matrix = new int[x.sequence().length() + 1][y.sequence().length() + 1];

		//////////////////// Algorithm
		// initialize first row and column with gap penalties
		for (int i = 0; i < matrix.length; i++) {
			matrix[i][0] = i * -1;
		}
		for (int j = 0; j < matrix[0].length; j++) {
			matrix[0][j] = j * -1;
		}

		// Check the three cases:
		// 1. Cell at i-1, j-1 + match/mismatch (look at the sequences at that index)
		// 2. Cell at i-1, j + gap penalty
		// 3. Cell at i, j-1 + gap penalty
		// Take the max of the three cases and put it in the cell at i, j
		for (int i = 1; i < matrix.length; i++) {
			for (int j = 1; j < matrix[0].length; j++) {
				int case1 = matrix[i - 1][j - 1]
						+ (x.sequence().charAt(i - 1) == y.sequence().charAt(j - 1) ? match_score : mismatch_score);
				int case2 = matrix[i - 1][j] - gap_penalty;
				int case3 = matrix[i][j - 1] - gap_penalty;
				matrix[i][j] = Math.max(case1, Math.max(case2, case3)); // NOTE: Math.max does not take multiple
																		// arguments, only two
			}
		}

		// Print out the optimal score and a corresponding alignment
		System.out.println("Optimal score: " + matrix[matrix.length - 1][matrix[0].length - 1]);

		//////////////////// TRACEBACK
		String xAligned = "";
		String yAligned = "";
		// From the last cell, go backwards through the matrix to find the alignment
		// Look at the three cases:
		// 1. Cell at i-1, j-1 - match/mismatch (look at the sequences at that index)
		// 2. Cell at i-1, j + gap penalty
		// 3. Cell at i, j-1 + gap penalty
		// If the cell is equal to case 1, add the characters at that index to each
		// alignment
		// If the cell is equal to case 2, add a gap to the y alignment and the
		// character at that index to the x alignment
		// If the cell is equal to case 3, add a gap to the x alignment and the
		// character at that index to the y alignment

		int i = matrix.length - 1;
		int j = matrix[0].length - 1;
		while (i > 0 && j > 0) {
			int case1 = matrix[i - 1][j - 1] + (x.sequence().charAt(i - 1) == y.sequence().charAt(j - 1) ? 1 : -1);
			int case2 = matrix[i - 1][j] - gap_penalty;
			int case3 = matrix[i][j - 1] - gap_penalty;
			if (matrix[i][j] == case1) {
				xAligned = x.sequence().charAt(i - 1) + xAligned;
				yAligned = y.sequence().charAt(j - 1) + yAligned;
				i--;
				j--;
			} else if (matrix[i][j] == case2) {
				xAligned = x.sequence().charAt(i - 1) + xAligned;
				yAligned = "-" + yAligned;
				i--;
			} else if (matrix[i][j] == case3) {
				xAligned = "-" + xAligned;
				yAligned = y.sequence().charAt(j - 1) + yAligned;
				j--;
			}
		}
		System.out.println("Alignment: ");
		System.out.println(xAligned);
		System.out.println(yAligned);

		///// Timing
		long end = System.currentTimeMillis();
		System.out.println("Time: " + (end - start) + " ms");

	}

	/**
	 * computes the optimal global alignment score and an alignment, using linear
	 * space. Prints out the optimal score and a corresponding alignment. Also
	 * prints out the number of milliseconds required for the computation
	 *
	 * @param x
	 * @param y
	 */
	public static void runNeedlemanWunschLinearSpace(FastA_GIESEL_MUEHLBAUER.Pair x, FastA_GIESEL_MUEHLBAUER.Pair y) {
		long start = System.currentTimeMillis();
		ArrayList<ArrayList<Integer>> indices = linear_space(x.sequence(), y.sequence(),0,0,true);
		System.out.println(indices);
		print_alignment(x,y,indices);
		long end = System.currentTimeMillis();
		System.out.println("Time: " + (end - start) + " ms");

	}

	/**
	 * Given two (sub)sequences and their offset from the start of the entire sequence, run the needleman wunsch algorithm to compute it in
	 * linear space, returns the paired indices
	 *
	 * @param x
	 * @param y
	 * @param x_offset
	 * @param y_offset
	 */
	public static ArrayList<ArrayList<Integer>> linear_space(String x, String y, int x_offset, int y_offset, Boolean first_call) {

		// Initialise the constants used throughout the program
		int gap_penalty = 1;
		int match_score = 1;
		int mismatch_score = -1;

		int xLength = x.length();
		int yLength = y.length();

		int middle_column = Math.ceilDiv(xLength, 2);

		//////////////////////// Initialise arrays to be used
		// make two arrays that are equal to the number for rows, if x is number columns
		int[] prev_column = new int[yLength + 1];
		int[] current_column = new int[yLength + 1];
		// keep two more array, relating to the rows that we came from (linear-space
		// c-matrix)
		int[] c_prev = new int[yLength + 1];
		int[] c_current = new int[yLength + 1];

		// initialise the first prev column with index * gap_penalty
		for (int i = 0; i < prev_column.length; i++) {
			prev_column[i] = i * -gap_penalty;
		}

		//////////////////////// Calculate the actual matrix/columns
		// Loop through the x sequence (i.e number of columns) and compute the current
		// column. i = columns, j = rows
		// For this, the first value of current_column is current index*gap_penalty
		// Then, for each value in the current column, compute the three cases:
		// 1. Cell at i-1, j-1 + match/mismatch (look at the sequences at that index)
		// 2. Cell at i-1, j - gap penalty // left
		// 3. Cell at i, j-1 - gap penalty // up
		// Take the max of the three cases and put it in the cell at i, j
		// Then, set the current column to the previous column and repeat

		for (int i = 1; i < xLength + 1; i++) {
			// middle-column shenanigans
			// check if middle column has been reached (+ 1 since we initialise prev_column
			// to be the current row number)
			if (i == middle_column + 1) {
				// if we have, initialise c_prev with row numbers
				for (int j = 0; j < c_prev.length; j++) {
					c_prev[j] = j;
				}
			}

			// ########Normal Needlman Wunsch#########
			// Initialise first values in each column
			current_column[0] = i * -gap_penalty;
			c_current[0] = 0;

			// loop through rows
			for (int j = 1; j < yLength + 1; j++) {
				int case1 = prev_column[j - 1]
						+ (x.charAt(i - 1) == y.charAt(j - 1) ? match_score : mismatch_score);
				int case2 = prev_column[j] - gap_penalty;
				int case3 = current_column[j - 1] - gap_penalty;
				int max = Math.max(case1, Math.max(case2, case3));
				current_column[j] = max;
				// Fill in the c_current array depending on which case was chosen
				if (current_column[j] == case1) {
					c_current[j] = c_prev[j - 1]; // look diagonal, previous column, previous row
				} else if (current_column[j] == case2) {
					c_current[j] = c_prev[j]; // look directly left i.e previous column
				} else if (current_column[j] == case3) {
					c_current[j] = c_current[j - 1]; // look directly up i.e previous row, same column
				}
			}

			// only swap before last column
			if (i < xLength) {
				prev_column = current_column;
				current_column = new int[yLength + 1];

				c_prev = c_current;
				c_current = new int[yLength + 1];
			}

		}
		if (first_call) {
			// print last value in current_column
			System.out.println("Optimal score: " + current_column[yLength]);
		}

		ArrayList<ArrayList<Integer>> indices = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> cur_ind = new ArrayList<Integer>();
		if (y.length()>=2){
			cur_ind.add(middle_column + x_offset);
			cur_ind.add(c_current[yLength] + y_offset);
			indices.add(cur_ind);
		}



		if (x.length() > 2 && y.length()>=2) {



			String new1_x = x.substring(0, middle_column);
			String new1_y = y.substring(0, c_current[yLength]);
			indices.addAll(linear_space(new1_x, new1_y, x_offset, y_offset, false));
			String new2_y;
			String new2_x;


			if (c_current[yLength] == 0 ){
				new2_x = x.substring(middle_column - 1, x.length());
				new2_y = y;
				x_offset += middle_column - 1;
			}
			else{
				new2_x = x.substring(middle_column - 1, x.length());
				new2_y = y.substring(c_current[yLength] - 1, y.length());
				x_offset += middle_column - 1;
				y_offset += c_current[yLength] - 1;
			}
			indices.addAll(linear_space(new2_x, new2_y, x_offset, y_offset, false));
		}
		return indices;

	}
	/**
	 * Given two sequences and a list of indices prints out the alignment
	 *
	 * @param x
	 * @param y
	 * @param indices
	 */
	public static void print_alignment(FastA_GIESEL_MUEHLBAUER.Pair x, FastA_GIESEL_MUEHLBAUER.Pair y, ArrayList<ArrayList<Integer>> indices) {
		String seqx = x.sequence();
		String seqy = y.sequence();
		String s1 = "";
		String s2 = "";
		int d = 0;
		ArrayList<Integer> y_used = new ArrayList<Integer>();
		for(int i = 1; i<=seqx.length(); i++){

			if (i%80 == 0){

				System.out.println(s1);
				System.out.println(s2);
				System.out.println();
				s1 = "";
				s2 = "";
			}

			for(int q = 0; q<indices.size();q++){
				if(indices.get(q).get(0) == i){
					int indy = indices.get(q).get(1);
					int com = i+d;
					if (indy == com){
						s1 += seqx.charAt(i-1);
						if(y_used.contains(indy)){
							break;
						}
						s2 += seqy.charAt(indy-1);
						y_used.add(indy);
					}
					else if(indy > com){
						int dif = indy - com;
						s1 += seqx.charAt(i-1);
						if(y_used.contains(com)){
							break;
						}
						s2 += seqy.charAt(com-1);
						y_used.add(com);
						for(int w = 0; w<dif; w++){
							s1 += "-";
							if(y_used.contains(com+w+1)){
								break;
							}
							s2 += seqy.charAt(com+w);
							y_used.add(com+w+1);
						}
						d += dif;
					}
					else if (com > indy){
						int dif = com - indy;

						s1 += seqx.charAt(i-1);
						if(!y_used.contains(indy)){
							s2 += seqy.charAt(i+d-1);
							y_used.add(i+d);
						}


						for(int w = 0; w<dif;w++){
							s1 += seqx.charAt(i+w);
							s2 += "-";
						}
						d -= dif;
					}
					break;
				}
				else if(q == indices.size()-1){
					if (i != seqx.length()){
						s1 += seqx.charAt(i-1);
						s2 += "-";
					}


				}
			}
		}
		for (int k = 1; k <= seqy.length(); k++){
			if (!y_used.contains(k)){
				s2 += seqy.charAt(k-1);
			}
		}

		System.out.println(s1);
		System.out.println(s2);

	}

	/**
	 * computes the optimal global alignment score using a recursion and no table.
	 * Prints out the optimal score. Also prints out the number of milliseconds
	 * required for the computation
	 *
	 * @param x
	 * @param y
	 */
	public static void runNeedlemanWunschRecursively(FastA_GIESEL_MUEHLBAUER.Pair x, FastA_GIESEL_MUEHLBAUER.Pair y) {
		long start = System.currentTimeMillis();
		int score = computeF(x.sequence().length()-1,y.sequence().length()-1,x.sequence(),y.sequence());
		long end = System.currentTimeMillis();
		System.out.println(score);
		System.out.println("Time: " + (end - start) + " ms");

	}

	public static int computeF(int i, int j, String x, String y) {
		int gap_penalty = 1;
		if(i == -1 && j == -1){
			return 0;
		}
		else if(i == -1){
			return -j*gap_penalty;
		}
		else if(j == -1){
			return -i*gap_penalty;
		}
		else{
			int m = Math.max(computeF(i-1, j-1, x, y) + s(x.charAt(i), y.charAt(j)), computeF(i-1, j, x, y)-gap_penalty);
			int max = Math.max(m, computeF(i, j-1, x, y)-gap_penalty);
			return max;
		}
	}
	public static int s(char xi, char yj){
		if(xi == yj){
			return 1;
		}
		else{
			return -1;
		}
	}
}
