package assignment08;

import assignment08.FastA_GIESEL_MUEHLBAUER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * proof-of-concept implementation for the basic minimap algorithm Sequence
 * Bioinformatics, WS 22/23 Minimap_GIESEL_MUEHLBAUER, 12.22
 */
public class Minimap_GIESEL_MUEHLBAUER {
	/**
	 * run the basic minimap algorithm
	 *
	 * @param args commandline arguments
	 * @throws IOException if arguments are incorrect or sequences not found
	 */
	public static void main(String[] args) throws IOException {
		int w;
		int k;
		ArrayList<FastA_GIESEL_MUEHLBAUER.Pair> targets;
		ArrayList<FastA_GIESEL_MUEHLBAUER.Pair> queries;
		if (args.length == 2) {
			w = 10;
			k = 15;
			targets = FastA_GIESEL_MUEHLBAUER.read(args[0]);
			queries = FastA_GIESEL_MUEHLBAUER.read(args[1]);
		} else if (args.length == 4) {
			w = Integer.parseInt(args[0]);
			k = Integer.parseInt(args[1]);
			targets = FastA_GIESEL_MUEHLBAUER.read(args[2]);
			queries = FastA_GIESEL_MUEHLBAUER.read(args[3]);
		} else
			throw new IOException("Usage: Minimap_GIESEL_MUEHLBAUER [w k] targets-file queries-file");

		if (true) { // this tests whether hashing and reverse-complement methods work as expected
			var testSequence = "CACGGTAGA";
			var hash = h(sk(testSequence, 1, 5, 0));
			if (hash != 107)
				throw new RuntimeException("Hashing broken, expected 107, got: " + hash);
			var hashReverseComplement = h(sk(testSequence, 1, 5, 1));
			if (hashReverseComplement != 91)
				throw new RuntimeException(
						"Hashing reverse-complement broken, expected 91, got: " + hashReverseComplement);
		}

		System.err.printf(Minimap_GIESEL_MUEHLBAUER.class.getSimpleName() + " w=%d k=%d targets=%d queries=%d%n", w, k,
				targets.size(), queries.size());

		var targetIndex = computeTargetIndex(targets, w, k);

		for (var record : queries) {
			var query = record.sequence();
			System.err.println("\nQuery: " + record.header());
			var matches = mapQuerySequence(targetIndex, query, w, k, 500);
			System.err.println("Matches: " + matches.size());
			for (var match : matches) {
				System.err.println(("Target: %d, query: %d - %d, target: %d - %d, reverse: %d".formatted(match.t() + 1,
						match.qMin() + 1, match.qMax() + 1, match.tMin() + 1, match.tMax() + 1, match.r())));
				System.err.println(query.substring(match.qMin(), match.qMax()));
				var target = targets.get(match.t()).sequence();
				System.err.println(sk(target, match.tMin(), match.tMax() - match.tMin(), match.r()));
			}
		}
	}

	/**
	 * extracts the k-mer at given position (as described in script)
	 *
	 * @param sequence the DNA sequence
	 * @param pos      the position
	 * @param k        the k-mer size
	 * @param r        0 for forward and 1 for reverse complement
	 * @return k-mer at given position or its reverse complement
	 */
	public static String sk(String sequence, int pos, int k, int r) {
		if (r == 0)
			return sequence.substring(pos, pos + k);
		else {
			var buf = new StringBuilder();
			// get the reverse compliment of the kmer string
			// start from end of kmer string and go backwards
			for (int i = pos + k - 1; i >= pos; i--) {
				if (sequence.charAt(i) == 'A') {
					buf.append('T');
				} else if (sequence.charAt(i) == 'T') {
					buf.append('A');
				} else if (sequence.charAt(i) == 'C') {
					buf.append('G');
				} else if (sequence.charAt(i) == 'G') {
					buf.append('C');
				}
			}

			return buf.toString();
		}
	}

	/**
	 * computes the h value for a k-mer (as described in script)
	 *
	 * @param s DNA string of length k
	 * @return h value
	 */
	public static int h(String s) {
		var value = 0;
		// Create a hashing algorithm
		// for a kmer (s) of length k, we can use the following formula
		// h(s) = s[0] * 4^(k-1) + s[1] * 4^(k-2) + ... + s[k-1] * 4^0
		// where s[i] is the i-th character of s.

		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == 'A') {
				value += 0 * Math.pow(4, s.length() - i - 1);
			} else if (s.charAt(i) == 'C') {
				value += 1 * Math.pow(4, s.length() - i - 1);
			} else if (s.charAt(i) == 'G') {
				value += 2 * Math.pow(4, s.length() - i - 1);
			} else if (s.charAt(i) == 'T') {
				value += 3 * Math.pow(4, s.length() - i - 1);
			}
		}

		return value;
	}

	/**
	 * computes a minimizer sketch for a given sequence and parameter k (algorithm 1
	 * in the script)
	 *
	 * @param s the DNA sequence
	 * @param w the word size
	 * @param k the k-mer size
	 * @return sorted set of all minimizers
	 */
	public static Set<Minimizer> minimizerSketch(String s, int w, int k) {
		var sketch = new HashSet<Minimizer>();

		for (int i = 0; i < s.length() - w - k + 1; i++) {
			var m = Double.MAX_VALUE;
			for (int j = 0; j <= w - 1; j++) {
				// no tuple destructuring in java :(
				int u = h(sk(s, i + j, k, 0));
				int v = h(sk(s, i + j, k, 1));
				if (u != v) {
					m = Math.min(m, Math.min(u, v));
				}
			}
			for (int j = 0; j <= w - 1; j++) {
				int u = h(sk(s, i + j, k, 0));
				int v = h(sk(s, i + j, k, 1));
				if (u < v && u == m) {
					sketch.add(new Minimizer((int) m, i + j, 0));
				} else if (u > v && v == m) {
					sketch.add(new Minimizer((int) m, i + j, 1));
				}
			}
		}

		return sketch;
	}

	/**
	 * Compute a hash map of h-values to minimizer locations in the target sequences
	 * (algorithm 3 in the script)
	 *
	 * @param targets the target sequences
	 * @param w       the word size
	 * @param k       the k-mer size
	 * @return the
	 */
	public static HashMap<Integer, Set<Location>> computeTargetIndex(ArrayList<FastA_GIESEL_MUEHLBAUER.Pair> targets,
			int w, int k) {
		var targetIndex = new HashMap<Integer, Set<Location>>();

		// todo: implement computation of target index as described in script (algorithm
		for (var t = 1; t < targets.size(); t++) {
			// get sequence & minimizer
			var target = targets.get(t).sequence();
			var minimizers = minimizerSketch(target, w, k);
			// iterate over minimizers and get h value
			for (var minimizer : minimizers) {
				if (!targetIndex.containsKey(minimizer.h())) {
					targetIndex.put(minimizer.h(), new HashSet<Location>());
				}
				targetIndex.get(minimizer.h()).add(new Location(t, minimizer.pos(), minimizer.r()));
			}
		}

		return targetIndex;
	}

	/**
	 * compute all matches of query to any of the target sequences (algorithm 4 in
	 * script)
	 *
	 * @param targetIndex the target index computed using #computeTargetIndex()
	 * @param query       the query DNA sequences
	 * @param w           the word size
	 * @param k           the k-mer size
	 * @param epsilon     the dialog shift allowed between two chained dialogs
	 */
	public static ArrayList<Match> mapQuerySequence(HashMap<Integer, Set<Location>> targetIndex, String query, int w,
			int k, int epsilon) {

		// compute array of k-mer hits:
		var A = new ArrayList<KMerHit>();
		var m = minimizerSketch(query, w, k);

		// todo: compute array of k-mer hits (as described in script, algorithm 4, part
		// 1)
		for (var minimizer : m) {
			// for each element in hashmap, check if h value is equal to minimizer h value
			for (var entry : targetIndex.entrySet()) {
				// if they're the same add the KmerHit(num of target seq,strand,num diagonal,pos
				// in
				// target) to the array
				var a = entry.getKey();
				var b = targetIndex.get(a);
				if (a == minimizer.h()) {
					A.add(new KMerHit(entry.getKey(), minimizer.r(), minimizer.pos(),
							entry.getValue().iterator().next().t()));
				}
			}
		}

		A.sort(KMerHit::compareTo);

		// chain k-mer hits into matches and return the matches
		var result = new ArrayList<Match>();
		var b = 0;
		for (var e = 0; e < A.size(); e++) {
			// todo: compute matches or ``clusters'' (as described in script, algorithm 4,
			// part 2)
			if (e == A.size() || A.get(e + 1).t != A.get(e).t || A.get(e + 1).c - A.get(e).c > epsilon
					|| A.get(e + 1).r != A.get(e).r) {
				// c = maximal colinear subset of A[b..e]
				var c = new ArrayList<KMerHit>();

				// loop over the array entries between index b and e to determine the minimum
				// and
				// maximum coordinates for both the query and the target sequence, and print or
				// return these.

				var minQuery = Integer.MAX_VALUE;
				var maxQuery = Integer.MIN_VALUE;
				var minTarget = Integer.MAX_VALUE;
				var maxTarget = Integer.MIN_VALUE;

				for (var i = b; i <= e; i++) {
					var kmerHit = A.get(i);
					if (kmerHit.pos < minQuery) {
						minQuery = kmerHit.pos;
					}
					if (kmerHit.pos > maxQuery) {
						maxQuery = kmerHit.pos;
					}
					if (kmerHit.pos < minTarget) {
						minTarget = kmerHit.pos;
					}
					if (kmerHit.pos > maxTarget) {
						maxTarget = kmerHit.pos;
					}
				}

				// add match to results
				result.add(new Match(A.get(e).t, A.get(e).r, minQuery, maxQuery, minTarget, maxTarget));
				b = e + 1;

			}

		}
		return result;
	}

	// PLEASE DO NOT CHANGE ANYTHING BELOW HERE:

	/**
	 * a minimizer
	 *
	 * @param h   hash value
	 * @param pos position in sequence
	 * @param r   0 for forward and 1 for reverse strand
	 */
	public static record Minimizer(int h, int pos, int r) {
	}

	/**
	 * a minimizer location
	 *
	 * @param t   index of target sequence 0...T-1, with T the number of target
	 *            sequences
	 * @param pos position in sequence t
	 * @param r   0 for forward and 1 for reverse strand
	 */
	public static record Location(int t, int pos, int r) {
	}

	/**
	 * a k-mer hit
	 *
	 * @param t   is the number of the target sequence
	 * @param r   is the *relative* strand (0 if both query and target minimizers on
	 *            same strand, otherwise 1)
	 * @param c   the dialog number
	 * @param pos the position in the target sequence
	 */
	public static record KMerHit(int t, int r, int c, int pos) implements Comparable<KMerHit> {
		@Override
		public int compareTo(KMerHit other) {
			if (t < other.t)
				return -1;
			else if (t > other.t)
				return 1;
			else if (r < other.r)
				return -1;
			else if (r > other.r)
				return 1;
			else if (c < other.c)
				return -1;
			else if (c > other.c)
				return 1;
			else
				return Integer.compare(pos, other.pos);
		}
	}

	/**
	 * a query-target sequence match
	 *
	 * @param t    the target index
	 * @param r    is the *relative* strand (0 for matches on same strand, 1 for on
	 *             different strands)
	 * @param qMin minimum position of match in query
	 * @param qMax maximum position of match in query
	 * @param tMin minimum position of match in target
	 * @param tMax maximum position of match in target
	 *
	 */
	public static record Match(int t, int r, int qMin, int qMax, int tMin, int tMax) {
	}
}
