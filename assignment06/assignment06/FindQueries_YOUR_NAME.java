package assignment06;

import assignment01.FastA_YOUR_NAME;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;

/**
 * find occurrences of queries in a text
 * Sequence Bioinformatics, WS 22/23
 * FindQueries_YOUR_NAME 11.22
 */
public class FindQueries_YOUR_NAME {
	public static void main (String[] args) throws IOException {
		System.out.println(FindQueries_YOUR_NAME.class.getSimpleName());

		if(args.length!=2)
			throw new IOException("Usage: FindQueries_YOUR_NAME text queries");

		var textItems= FastA_YOUR_NAME.read(args[0]);
		if(textItems.size()!=1)
			throw new IOException("text must contain 1 sequence, found: "+textItems.size());

		System.out.println("Text: "+textItems.get(0).sequence());

		var suffixTree=new NaiveSuffixTree(textItems.get(0).sequence());

		var queryItems=FastA_YOUR_NAME.read(args[1]);

		for(var item:queryItems) {
			System.out.println("Query "+item.sequence());
			System.out.println("Contained: "+contains(suffixTree,item.sequence()));
			System.out.print("Occurrences:");
			for(var pos:find(suffixTree,item.sequence())) {
				System.out.print(" "+pos);
			}
			System.out.println();
		}
	}

	/**
	 * determines whether text contains query
	 * @param suffixTree the suffix tree representing the text
	 * @param query the query
	 * @return true, if query in text
	 */
	public static boolean contains(NaiveSuffixTree suffixTree, String query) {
		// todo: please implement this
		return true;
	}

	/**
	 * find and return all occurrences of query in text
	 * @param suffixTree the suffix tree representing the text
	 * @param query the query
	 * @return all positions in text at which query occurs
	 */
	public static Collection<Integer> find(NaiveSuffixTree suffixTree, String query) {
		// todo: please implement this
		return Collections.emptyList();
	}
}
