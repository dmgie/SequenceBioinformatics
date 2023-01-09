

assembly() {

	input_file=$1
	output_file=$2
	# Minimap - create overlaps
	./minimap2 -x ava-pb -t8 $input_file $input_file | gzip -1 > "$output_file"

	### Explanation
	# -x ava-pb = pacbio reads format
	# -t = number of threads
	# gzip = gzip to make a compressed file

	# Miniasm - layout above overlaps
	# get output file name without extension
	output_file_name=$(basename "$input_file" | cut -d. -f1)
	./miniasm -f $input_file $output_file > "$output_file_name.gfa"

	### Explanation
	# -f = input file for the reads 
	# second file is the output from minimap (overlaps)
}


# Create read subset
make_subset() {
	NUM_READS=$1
	# Assume every read is made up of 4 lines (fastq format)
	NUM_LINES=$(($NUM_READS * 4))
	input_file=$2
	output_file=$3

	# Use zcat if it ends with gz
	if [[ $input_file == *.gz ]]; then
		zcat "$input_file" | head -n $NUM_READS | sed -n '1~4p' > "$output_file"
	else
		head -n $NUM_READS $input_file | sed -n '1~4p' > "$output_file"
	fi

}

get_statistics() {
	# Get statistics from the assembly
	input_file=$1
	# Get the number of contigs (gfa file, miniasm output)
	num_contigs=$(grep -c "S" $input_file)
	# Get the total length of the num_contigs
	total_length=$(grep "S" $input_file | awk '{sum+=$3} END {print sum}')
	# Get the N50
	N50=$(awk -v total=$total_length 'BEGIN {print total/2}')
	echo "Number of contigs: $num_contigs"
	echo "Total length: $total_length"
	echo "N50: $N50"

	# Or also use quast
	# quast.py $input_file -o quast_output

}

# Run normally
# assembly pacbio-reads.fastq.gz pacbio-reads-overlaps.paf.gz
get_statistics pacbio-reads.gfa

# # Make subset
# make_subset 50000 pacbio-reads.fastq.gz pacbio-50000.fastq

# # Run on subset
# assembly pacbio-50000.fastq.gz pacbio-50000-overlaps.paf.gz
# get_statistics pacbio-50000-overlaps.paf.gz 


