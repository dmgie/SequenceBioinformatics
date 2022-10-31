package assignment02;

import assignment02.FastA_YOUR_NAME;
import java.io.IOException;
import java.util.List;

// TODO: Implement as row-major order vector instead of 2D array

/**
 * GlobalAligner_YOUR_NAME Sequence Bioinformatics, WS 22/23
 */
public class GlobalAligner_YOUR_NAME {
	int gap_penalty = 1;
	int match_score = 1;
	int mismatch_score = -1;

	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2)
			throw new IOException("Usage: GlobalAligner_YOUR_NAME infile [quadraticSpace|linearSpace|noDP]");

		var list = FastA_YOUR_NAME.read(args[0]);

		if (list.size() != 2)
			throw new IOException("Wrong number of input sequences: " + list.size());

		var mode = (args.length == 2 ? args[1] : "quadraticSpace");

		switch (mode) {
		case "quadraticSpace" -> runNeedlemanWunschQuadraticSpace(list.get(0), list.get(1));
		case "linearSpace" -> runNeedlemanWunschLinearSpace(list.get(0), list.get(1));
		case "noDP" -> runNeedlemanWunschQuadraticSpace(list.get(0), list.get(1));
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
	public static void runNeedlemanWunschQuadraticSpace(FastA_YOUR_NAME.Pair x, FastA_YOUR_NAME.Pair y) {
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

		// print matrix, for debugging
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				System.out.print(matrix[i][j] + " ");
			}
			System.out.println();
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
	public static void runNeedlemanWunschLinearSpace(FastA_YOUR_NAME.Pair x, FastA_YOUR_NAME.Pair y) {
		// todo: implement, Assignment 2.2
		// initialise a list of list of size 2

		// Get the middle column

		linear_space(x, y);
	}

	/**
	 * Given two sequences, run the needleman wunsch algorithm to compute it in
	 * linear space
	 *
	 * @param x
	 * @param y
	 */
	public static void linear_space(FastA_YOUR_NAME.Pair x, FastA_YOUR_NAME.Pair y) {
		// Initialise the constants used throughout the program
		int gap_penalty = 1;
		int match_score = 1;
		int mismatch_score = -1;

		int xLength = x.sequence().length();
		int yLength = y.sequence().length();

		int middle_column = Math.floorDiv(xLength, 2);

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
						+ (x.sequence().charAt(i - 1) == y.sequence().charAt(j - 1) ? match_score : mismatch_score);
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
		// print last value in current_column
		System.out.println("Optimal score: " + current_column[yLength]);
		// get the index of
		System.out.println("Index of traceback: " + c_current[yLength] + " " + middle_column);
	}

	/**
	 * computes the optimal global alignment score using a recursion and no table.
	 * Prints out the optimal score. Also prints out the number of milliseconds
	 * required for the computation
	 *
	 * @param x
	 * @param y
	 */
	public static void runNeedlemanWunschRecursively(FastA_YOUR_NAME.Pair x, FastA_YOUR_NAME.Pair y) {
		// todo: implement using recursive function computeF, , Assignment 2.3

	}

	public static int computeF(int i, int j) {
		// todo: implement
		return 0;
	}
}
