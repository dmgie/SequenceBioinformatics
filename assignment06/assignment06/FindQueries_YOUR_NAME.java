package assignment06;

import assignment06.FastA_GIESEL_MUEHLBAUER;

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

		var textItems= FastA_GIESEL_MUEHLBAUER.read(args[0]);
		if(textItems.size()!=1)
			throw new IOException("text must contain 1 sequence, found: "+textItems.size());

		System.out.println("Text: "+textItems.get(0).sequence());

		var suffixTree=new NaiveSuffixTree(textItems.get(0).sequence());

		var queryItems=FastA_GIESEL_MUEHLBAUER.read(args[1]);

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
		if (query.length()<1){
			throw new RuntimeException("No query given");
		}
		NaiveSuffixTree.Node root = suffixTree.getRoot();
		suffixTree.printTree();
		NaiveSuffixTree.Node child = root.getChild(query.charAt(0));
		return nodeContains(child,query);
	}

	public static boolean nodeContains(NaiveSuffixTree.Node node, String query){
		String label = node.getLetters();
		if(query.length()>label.length()){
			if(label.equals(query.substring(0,label.length()))){
				String newQuery = query.substring(label.length()-1);
				NaiveSuffixTree.Node child = node.getChild(newQuery.charAt(0));
				return nodeContains(child,newQuery);
			}
			else{
				return false;
			}
		}
		else{
			if(query.equals(label.substring(0,query.length()))){
				return true;
			}
			else{
				return false;
			}
		}
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
