package assignment07;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * implementation of mash sketch and distance calculation
 * Sequence Bioinformatics, WS 22/23
 * Mash_YOUR_NAME
 */
public class Mash_YOUR_NAME {
    public static void main(String[] args) throws IOException {
        if (args.length != 3)
            throw new IOException("Usage: Mash_YOUR_NAME k s directory-of-input-files");
        else
            System.err.printf("Mash k=%s s=%s input-directory=%s%n", args[0], args[1], args[2]);

        // get k and s:
        var k = Integer.parseInt(args[0]);
        var s = Integer.parseInt(args[1]);

        // load all FastA files - use the implementation of FastA provided at the bottom of the file
        var genomes = new ArrayList<FastA>();
        var directory = args[2];
        var files = (new File(directory)).list();
        if (files == null || files.length == 0)
            throw new IOException("Not a directory or no files found: " + args[2]);

        for (var file : files) {
            genomes.add(FastA.merge(FastA.load(directory + File.separator + file)));
        }
        // don't change anything above here (except for adding your name!)

        // represent a sketch as a sorted set of integers
        // the following three lines compute a sketch for each genome, for the given k-mer size and sketch size,
        // using a method called computeSketch(), which you need to implement below

        var sketches = genomes.stream() // replace stream() by parallelStream() to run your algorithm in parallel!
                .map(genome -> computeSketch(k, s, genome.getSequence()))
                .collect(Collectors.toCollection(ArrayList::new));
        // todo: task 1: implement the method computeSketch() below


        // compute the Jaccard index for all pairs of genomes:
        // todo: task 2: please implement the method computeJaccardIndex() below

        var jaccard = new double[genomes.size()][genomes.size()];
        for (var i = 0; i < genomes.size(); i++) {
            for (var j = i + 1; j < genomes.size(); j++) {
                jaccard[i][j] = jaccard[j][i] = computeJaccardIndex(s, sketches.get(i), sketches.get(j));
            }
        }


        // compute the Mash distances from the Jaccard indices, as discussed in the lecture:
        // todo: task 3: please implement the method computeMashDistance()

        var mashDistance = new double[genomes.size()][genomes.size()];
        for (var i = 0; i < genomes.size(); i++) {
            for (var j = i + 1; j < genomes.size(); j++) {
                mashDistance[i][j] = computeMashDistance(k, jaccard[i][j]);
            }
        }

        // output in format that can be read by SplitsTree6:
        System.out.println(mashDistance.length);
        for (var i = 0; i < mashDistance.length; i++) {
            System.out.printf("%s", genomes.get(i).getHeader());
            for (var j = 0; j < mashDistance.length; j++) {
                System.out.printf("\t%.5f", mashDistance[i][j]);
            }
            System.out.println();
        }
    }

    /**
     * computes a bottom sketch of size s for k-mers.
     * Use String method substring() to extract k-mer and hashCode() to compute k-mer hash value
     * Implement the method computeReverseComplement() to compute the reverse complement of a k-mer
     *
     * @param k
     * @param s
     * @param genome
     * @return sketch
     */
    public static SortedSet<Integer> computeSketch(int k, int s, String genome) {
        final SortedSet<Integer> sketch = new TreeSet<>();
        int hash;
        String g = genome.toUpperCase();
        for(int i = 0; i <= g.length()-k; i++){
            int j = i+k;
            String kmer = g.substring(i,j);
            String rev_comp_kmer = computeReverseComplement(kmer);
            int kmer_hash = kmer.hashCode();
            int rc_hash = rev_comp_kmer.hashCode();
            if(kmer_hash<rc_hash){
                hash = kmer_hash;
            }
            else{
                hash = rc_hash;
            }
            if(sketch.size()<s){
                sketch.add(hash);
            }
            else if(sketch.last() > hash){
                sketch.remove(sketch.last());
                sketch.add(hash);
            }
        }
        System.out.println(sketch);
        return sketch;
    }

    /**
     * gets the reverse compliment of a sequence
     *
     * @param dna
     * @return reverse complement
     */
    public static String computeReverseComplement(String dna) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < dna.length(); i++ ){
            if(dna.charAt(i) == 'A'){
                sb.append('T');
            }
            else if(dna.charAt(i) == 'T'){
                sb.append('A');
            }
            else if(dna.charAt(i) == 'G'){
                sb.append('C');
            }
            else if(dna.charAt(i) == 'C'){
                sb.append('G');
            }
        }

        sb.reverse();
        return sb.toString();
    }


    /**
     * compute the Jaccard index for two sketches
     *
     * @param s
     * @param sketchA
     * @param sketchB
     * @return Jaccard index
     */
    public static double computeJaccardIndex(int s, SortedSet<Integer> sketchA, SortedSet<Integer> sketchB) {
        // todo: please implement
        final SortedSet<Integer> sketch_union = new TreeSet<>();
        sketch_union.addAll(sketchB);
        final SortedSet<Integer> sketch_intersect = new TreeSet<>();
        //compute union and intersection
        for(int hashA: sketchA){
            if(!sketch_union.contains(hashA)&&hashA < sketch_union.last()){
                sketch_union.remove(sketch_union.last());
                sketch_union.add(hashA);
            }
        }
        for(int hashA: sketchA){
            if(sketchB.contains(hashA) && sketch_union.contains(hashA)){
                sketch_intersect.add(hashA);
            }
        }
        double z = sketch_intersect.size();
        double n = s;
        double ji = z/n;
        System.out.println("JI:");
        System.out.println(ji);
        return ji;
    }

    /**
     * computes the mash distance from the Jaccard index
     *
     * @param k
     * @param jaccardIJ
     * @return mash distance
     */
    public static double computeMashDistance(int k, double jaccardIJ) {
        // todo: please implement
        double n_k = k;
        double d = -1/n_k * Math.log10(2*jaccardIJ/(1+jaccardIJ));
        return d;
    }

    // Code for reading FastA, please don't change anything below here

    public static class FastA {
        private final String header;
        private final String sequence;

        public FastA(String header, String sequence) {
            this.header = header;
            this.sequence = sequence;
        }

        public String getHeader() {
            return header;
        }

        public String getSequence() {
            return sequence;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof FastA fastA)) return false;
            return Objects.equals(header, fastA.header) &&
                    Objects.equals(sequence, fastA.sequence);
        }

        @Override
        public int hashCode() {
            return Objects.hash(header, sequence);
        }

        /**
         * load fastA records from a file
         *
         * @param fileName
         * @return list of loaded records
         * @throws IOException
         */
        public static ArrayList<FastA> load(String fileName) throws IOException {
            var result = new ArrayList<FastA>();
            var headerBuffer = new StringBuilder();
            var sequenceBuffer = new StringBuilder();

            Files.lines(new File(fileName).toPath()).forEach(line -> {
                if (line.startsWith(">")) {
                    if (headerBuffer.length() > 0) {
                        result.add(new FastA(headerBuffer.toString().replaceAll("^>\\s*", ""), sequenceBuffer.toString()));
                    }
                    headerBuffer.setLength(0);
                    headerBuffer.append(line);
                    sequenceBuffer.setLength(0);
                } else
                    sequenceBuffer.append(line.replaceAll("\\s*", ""));
            });
            if (headerBuffer.length() > 0)
                result.add(new FastA(headerBuffer.toString().replaceAll("^>\\s*", ""), sequenceBuffer.toString()));

            if (result.size() == 0)
                throw new IOException("No sequences found in file: " + fileName);

            return result;
        }

        public static FastA merge(ArrayList<FastA> list) {
            if (list.size() == 0)
                return null;
            var header = list.get(0).getHeader();
            var sequence = list.stream().map(FastA::getSequence).collect(Collectors.joining());
            return new FastA(header, sequence);
        }
    }
}


