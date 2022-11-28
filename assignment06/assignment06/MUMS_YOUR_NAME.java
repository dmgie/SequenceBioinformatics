package assignment06;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * computes all MUMs in a text
 * Sequence Bioinformatics, WS 22/23
 * MUMS_YOUR_NAME, 11.22
 */
public class MUMS_YOUR_NAME {

	public static void main (String[] args) throws IOException {
		System.out.println(MUMS_YOUR_NAME.class.getSimpleName());

		if (args.length != 1)
			throw new IOException("Usage: MUMS_YOUR_NAME fasta-file");

		var textItems = FastA_GIESEL_MUEHLBAUER.read(args[0]);
		if (textItems.size() != 2)
			throw new IOException("fasta-file must contain 2 sequence, found: " + textItems.size());

		// todo: please implement

		// 1. setup suffix tree for appropriately concatenated sequences
		String seq1 = textItems.get(0).sequence();
		String seq2 = textItems.get(1).sequence();
		var suffixTree=new NaiveSuffixTree(seq1+"%"+seq2);

		// 2. implement algorithm to report all MUMs (of any size)
		int end1 = seq1.length();
		NaiveSuffixTree.Node root = suffixTree.getRoot();
		HashMap<String, ArrayList<Integer>> potMUMs = findPotMUMs(root,end1,"");
		System.out.println(potMUMs);



		// output should be:
		// MUM "GC" at 2 and 2 (1-based)
		// for input AGCT and GGCC
	}
	public static HashMap<String, ArrayList<Integer>> findPotMUMs(NaiveSuffixTree.Node node, int end1, String label){
		ArrayList<NaiveSuffixTree.Node> directChild = new ArrayList<NaiveSuffixTree.Node>();
		HashMap<String, ArrayList<Integer>> potMUMs = new HashMap<String, ArrayList<Integer>>();
		for (var child : node.getChildren()) {
			int suffPos = child.getSuffixPos();
			if (suffPos != -1){
					directChild.add(child);
			}
			else{
				String lab = label + child.getLetters();
				potMUMs.putAll(findPotMUMs(child,end1,lab));
			}
		}
		if (directChild.size() == 2){
			int suf1 = directChild.get(0).getSuffixPos();
			int suf2 = directChild.get(1).getSuffixPos();

			boolean difSeq = (suf1<end1&&suf2>end1)||(suf1>end1&&suf2<end1);
			if(difSeq){
				ArrayList<Integer> potMUM = new ArrayList<Integer>();
				potMUM.add(suf1);
				potMUM.add(suf2);

				potMUMs.put(label, potMUM);
			}
		}
		return potMUMs;
	}
}
