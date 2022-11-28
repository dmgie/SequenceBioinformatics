package assignment06;

import java.util.*;

/**
 * Naive implementation of suffix tree
 * PLEASE DON'T CHANGE THE IMPLEMENTATION OF THIS (UNLESS TO FIX BUGS)
 * Sequence Bioinformatics, WS 22/23
 * Daniel_Huson 11.22
 */
public class NaiveSuffixTree {
	private final Node root;

	/**
	 * constructs a suffix tree for the given text
	 * @param text the text
	 */
	public NaiveSuffixTree(String text) {
		if(!text.endsWith("$"))
			text = text + "$";
		root = new Node("", -1);

		for (var suffixPos = 0; suffixPos < text.length(); suffixPos++) {
			var node = root;
			var pos = suffixPos;

			while (true) {
				var firstLetter = text.charAt(pos);
				var next = node.getChild(firstLetter);

				var letters = text.substring(pos);
				if (next == null) {
					node.setChild(firstLetter, new Node(letters, suffixPos));
					break;
				} else {
					var nextLetters = next.getLetters();
					var countSame = 0;
					var top = Math.min(letters.length(), nextLetters.length());
					while (countSame < top) {
						if (letters.charAt(countSame) != nextLetters.charAt(countSame))
							break;
						countSame++;
					}
					if (countSame == letters.length()) // can't happen because of unique terminal $ at end of string
						throw new RuntimeException("Suffix already contained");
					if (countSame < nextLetters.length()) { // need to insert new node and leaf edge
						var middle = new Node(letters.substring(0, countSame));
						var leaf = new Node(letters.substring(countSame), suffixPos);
						node.setChild(firstLetter, middle);
						middle.setChild(letters.charAt(countSame), leaf);
						var remaining=nextLetters.substring(countSame);
						next.setLetters(remaining);
						middle.setChild(remaining.charAt(0),next);
						break;
					} else { // move to next node
						pos += countSame;
						node = next;
					}
				}
			}
		}
	}

	/**
	 * accesses the root of the tree
	 * @return root
	 */
	public Node getRoot() {
		return root;
	}

	/**
	 * gets all children of a node
	 * @param node the node
	 * @return array list containing all children
	 */
	public ArrayList<Node> getChildren(Node node) {
		return new ArrayList<>(node.getChildren());
	}

	/**
	 * prints the suffix tree
	 */
	public void printTree() {
		root.printSubTreeRec( 0);
	}

	/**
	 * a suffix tree node, with children, letters on edge to node and position of
	 * suffix, if leaf
	 */
	public static final class Node {
		private final Map<Character, Node> children = new HashMap<>();
		private String letters;
		private final int suffixPos;

		/**
		 * constructs an internal node
		 *
		 * @param letters edge label
		 */
		private Node(String letters) {
			this(letters, -1);
		}

		/**
		 * constructs a leaf node
		 *
		 * @param letters   edge label
		 * @param suffixPos suffix position
		 */
		private Node(String letters, int suffixPos) {
			this.letters = letters;
			this.suffixPos = suffixPos;
		}

		/**
		 * set a child
		 *
		 * @param ch   first letter of label to child
		 * @param node child
		 */
		public void setChild(Character ch, Node node) {
			children.put(ch, node);
		}

		/**
		 * get a child
		 *
		 * @param ch first letter of label to child
		 * @return child or null
		 */
		public Node getChild(Character ch) {
			return children.get(ch);
		}

		/**
		 * gets letters on label to node
		 *
		 * @return letters
		 */
		public String getLetters() {
			return letters;
		}

		/**
		 * set the letters
		 * @param letters new letters
		 */
		public void setLetters(String letters) {
			this.letters=letters;
		}

		/**
		 * gets all children
		 *
		 * @return children
		 */
		public Collection<Node> getChildren() {
			return children.values();
		}

		/**
		 * gets suffix position
		 *
		 * @return suffix position or -1
		 */
		public int getSuffixPos() {
			return suffixPos;
		}

		@Override
		public String toString() {
			return letters+" "+children+" "+suffixPos;
		}

		/**
		 * recursively print the subtree rooted at the given node
		 * @param level used for indenting
		 */
		private void printSubTreeRec(int level) {
			System.out.print("->".repeat(level) +"\""+letters+"\"");
			if(suffixPos!=-1)
				System.out.print(" pos: "+suffixPos);
			System.out.println();
			for (var child : getChildren()) {
				child.printSubTreeRec(level + 1);
			}
		}
	}
}
