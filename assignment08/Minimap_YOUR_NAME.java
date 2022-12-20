package assignment08;

import assignment08.FastA_YOUR_NAME;

import java.io.IOException;
import java.util.*;

/**
 * proof-of-concept implementation for the basic minimap algorithm
 * Sequence Bioinformatics, WS 22/23
 * Minimap_YOUR_NAME, 12.22
 */
public class Minimap_YOUR_NAME {
	/**
	 * run the basic minimap algorithm
	 * @param args commandline arguments
	 * @throws IOException if arguments are incorrect or sequences not found
	 */
	public static void main(String[] args) throws IOException {
		int w;
		int k;
		ArrayList<FastA_YOUR_NAME.Pair> targets;
		ArrayList<FastA_YOUR_NAME.Pair> queries;
		if(args.length==2) {
			w=10;
			k=15;
			targets=FastA_YOUR_NAME.read(args[0]);
			queries=FastA_YOUR_NAME.read(args[1]);
		}
		else if(args.length==4) {
			 w=Integer.parseInt(args[0]);
			 k=Integer.parseInt(args[1]);
			 targets= FastA_YOUR_NAME.read(args[2]);
			 queries=FastA_YOUR_NAME.read(args[3]);
		}
		else
			throw new IOException("Usage: Minimap_YOUR_NAME [w k] targets-file queries-file");


		if(true) { // this tests whether hashing and reverse-complement methods work as expected
			var testSequence= "CACGGTAGA";
			var hash=h(sk(testSequence, 1, 5, 0));
			if (hash != 107)
				throw new RuntimeException("Hashing broken, expected 107, got: "+hash);
			var hashReverseComplement=h(sk(testSequence, 1, 5, 1));
			if (hashReverseComplement != 91)
				throw new RuntimeException("Hashing reverse-complement broken, expected 91, got: "+hashReverseComplement);
		}

		System.err.printf(Minimap_YOUR_NAME.class.getSimpleName() + " w=%d k=%d targets=%d queries=%d%n", w,k,targets.size(),queries.size());

		var targetIndex=computeTargetIndex(targets,w,k);

		for(var record:queries) {
			var query= record.sequence();
			System.err.println("\nQuery: "+record.header());
			var matches=mapQuerySequence(targetIndex,query,w,k,500);
			System.err.println("Matches: "+matches.size());
			for(var match:matches) {
				System.err.println(("Target: %d, query: %d - %d, target: %d - %d, reverse: %d"
						.formatted(match.t()+1,match.qMin()+ 1, match.qMax()+ 1,match.tMin()+ 1, match.tMax()+ 1,match.r())));
				System.err.println(query.substring(match.qMin(), match.qMax()));
				var target=targets.get(match.t()).sequence();
				System.err.println(sk(target,match.tMin(),match.tMax()-match.tMin(), match.r()));
			}
		}
	}

	/**
	 * extracts the k-mer at given position  (as described in script)
	 * @param sequence the DNA sequence
	 * @param pos the position
	 * @param k the k-mer size
	 * @param r 0 for forward and 1 for reverse complement
	 * @return k-mer at given position or its reverse complement
	 */
	public static String sk(String sequence,int pos,int k,int r) {
		String kmer = sequence.substring(pos,pos+k);
		if(r==0)
			return kmer;
		else {
			var buf=new StringBuilder();
			// todo: implement reverse-complement of k-mer here
			for(int i = 0; i < kmer.length(); i++ ){
				if(kmer.charAt(i) == 'A'){
					buf.append('T');
				}
				else if(kmer.charAt(i) == 'T'){
					buf.append('A');
				}
				else if(kmer.charAt(i) == 'G'){
					buf.append('C');
				}
				else if(kmer.charAt(i) == 'C'){
					buf.append('G');
				}
			}
			buf.reverse();
			return buf.toString();
		}
	}

	/**
	 * computes the h value for a k-mer (as described in script)
	 * @param s DNA string of length k
	 * @return h value
	 */
	public static int h(String s) {
		var value = 0;
		// todo: implement hashing as described in script
		int len = s.length();
		for(int i = 0; i < len; i++ ) {
			if(s.charAt(i) == 'C'){
				value += Math.pow(4,len-i-1);
			}
			else if(s.charAt(i) == 'G'){
				value += 2 * Math.pow(4,len-i-1);
			}
			else if(s.charAt(i) == 'T'){
				value += 3 * Math.pow(4,len-i-1);
			}
		}
		return value;
	}

	/**
	 * computes a minimizer sketch for a given sequence and parameter k (algorithm 1 in the script)
	 * @param s the DNA sequence
	 * @param w the word size
	 * @param k the k-mer size
	 * @return sorted set of all minimizers
	 */
	public static Set<Minimizer> minimizerSketch(String s, int w, int k) {
		var sketch=new HashSet<Minimizer>();
		int l = s.length();
		// todo: implement computation of minimizer sketch as described in script (algorithm 1)
		for(int i = 1; i<=l-w-k+1;i++){
			double m = Double.POSITIVE_INFINITY;
			for(int j = 0; j<=w-1; j++){
				int u = h(sk(s,i+j,k,0));
				int v = h(sk(s,i+j,k,1));
				if(u!=v){
					m = Math.min(m,Math.min(u,v));
				}
			}
			for(int j = 0; j <= w-1; j++){
				int u = h(sk(s,i+j,k,0));
				int v = h(sk(s,i+j,k,1));
				if(u<v&&u==m){
					sketch.add(new Minimizer((int) m,i+j,0));
				}
				else if(u>v&&v==m){
					sketch.add(new Minimizer((int) m,i+j,1));
				}
			}
		}

		return sketch;
	}

	/**
	 * Compute a hash map of h-values to minimizer locations in the target sequences (algorithm 3 in the script)
	 * @param targets the target sequences
	 * @param w the word size
	 * @param k the k-mer size
	 * @return the
	 */
	public static HashMap<Integer,Set<Location>> computeTargetIndex(ArrayList<FastA_YOUR_NAME.Pair> targets, int w, int k) {
		var targetIndex= new HashMap<Integer,Set<Location>>();
		int T = targets.size();

		// todo: implement computation of target index as described in script (algorithm 1)
		for(int t = 0; t<T; t++){
			var M = minimizerSketch(targets.get(t).sequence(), w,k);
			for (Minimizer minimizer: M) {
				Set<Location> hSet = targetIndex.get(minimizer.h);
				Location loc = new Location(t,minimizer.pos, minimizer.r);
				if(Objects.isNull(hSet)){
					hSet = new HashSet<Location>();
				}
				hSet.add(loc);
				targetIndex.put(minimizer.h, hSet);
			}
		}

		return targetIndex;
	}

	/**
	 * compute all matches of query to any of the target sequences (algorithm 4 in script)
	 * @param targetIndex the target index computed using #computeTargetIndex()
	 * @param query the query DNA sequences
	 * @param w the word size
	 * @param k the k-mer size
	 * @param epsilon the dialog shift allowed between two chained dialogs
	 */
	public static ArrayList<Match> mapQuerySequence(HashMap<Integer,Set<Location>> targetIndex, String query, int w, int k, int epsilon) {

		// compute array of k-mer hits:
		var A=new ArrayList<KMerHit>();

		// todo: compute array of k-mer hits (as described in script, algorithm 4, part 1)
		Set<Minimizer> M = minimizerSketch(query,w,k);

		for(Minimizer min : M){
			if(!Objects.isNull(targetIndex.get(min.h))){
				for (Location loc: targetIndex.get(min.h)){
					if(min.r==loc.r){
						A.add(new KMerHit(loc.t,0,min.pos-loc.pos,loc.pos));
					}
					else{
						A.add(new KMerHit(loc.t,1,min.pos+loc.pos,loc.pos));
					}
				}
			}
		}

		A.sort(KMerHit::compareTo);

		// chain k-mer hits into matches and return the matches
		var result=new ArrayList<Match>();
		var b=0;
		for(var e=0;e<A.size();e++) {
			// todo: compute matches or ``clusters'' (as described in script, algorithm 4, part;s 2 and 3

			if(e+1 == A.size()||A.get(e+1).t!=A.get(e).t||A.get(e+1).r!=A.get(e).r||A.get(e+1).c-A.get(e).c >= epsilon){
				int bq = 0;
				int eq = 0;
				int bt = 0;
				int et = 0;
				if(A.get(e).r == 0){
					bt = A.get(b).pos;
					et = A.get(e).pos+k;
					bq = A.get(b).c + A.get(b).pos;
					eq = A.get(e).c + A.get(e).pos+k;
				}
				else{
					bt = A.get(b).pos;
					et = A.get(e).pos+k;
					bq = A.get(e).c - A.get(e).pos;
					eq = A.get(b).c - A.get(b).pos+k;
				}
				Match C = new Match(A.get(e).t, A.get(e).r,bq,eq,bt,et);
				b = e+1;
				result.add(C);
			}
		}
		return result;
	}


	// PLEASE DO NOT CHANGE ANYTHING BELOW HERE:

	/**
	 * a minimizer
	 * @param h hash value
	 * @param pos position in sequence
	 * @param r 0 for forward and 1 for reverse strand
	 */
	public static record Minimizer(int h,int pos, int r) {
	}

	/**
	 * a minimizer location
	 * @param t index of target sequence 0...T-1, with T the number of target sequences
	 * @param pos position in sequence t
	 * @param r 0 for forward and 1 for reverse strand
	 */
	public static record Location(int t,int pos,int r) {}

	/**
	 * a k-mer hit
	 * @param t is the number of the target sequence
	 * @param r is the *relative* strand (0 if both query and target minimizers on same strand, otherwise 1)
	 * @param c the dialog number
	 * @param pos the position in the target sequence
	 */
	public static record KMerHit(int t, int r, int c, int pos) implements Comparable<KMerHit>{
		@Override
		public int compareTo(KMerHit other) {
			if(t<other.t)
				return -1;
			else if(t>other.t)
				return 1;
			else if(r<other.r)
				return -1;
			else if(r>other.r)
				return 1;
			else if(c<other.c)
				return -1;
			else if(c>other.c)
				return 1;
			else return Integer.compare(pos, other.pos);
		}
	}

	/**
	 * a query-target sequence match
	 * @param t the target index
	 * @param r is the *relative* strand (0 for matches on same strand, 1 for on different strands)
	 * @param qMin minimum position of match in query
	 * @param qMax maximum position of match in query
	 * @param tMin minimum position of match in target
	 * @param tMax maximum position of match in target

	 */
	public static record Match(int t, int r, int qMin, int qMax, int tMin, int tMax) {
	}
}

