# Navigate to the directory containing the raw data for the CRISPRi experiment
cd /home/ryandward/R/acinetobacter_baumannii_CRISPRi_seq/Raw_Data

# Iterate over all gzip-compressed files in the current directory that contain "R1" in the name, corresponding to the forward reads
for i in *R1*gz; do
	# Assign the file name to the forward variable
	export forward=$i
	
	# Assign the file name with "R1" replaced with "R2" to the reverse variable, corresponding to the reverse reads
	export reverse=$(echo $forward | sed 's/R1/R2/g')
	
	# Assign the file name with everything after the first underscore removed and ".seal.tsv" appended to the end to the set variable
	# This will be used as the output file name for the seal.sh script
	export set=$(echo $forward | sed 's/\_.*//g').seal.tsv
	
	# Echo a command to run the seal.sh script, which is a bbtools binary file for counting barcodes in a CRISPRi experiment
	# The seal.sh script is given various options, including the input files for the forward and reverse reads (in1=$forward, in2=$reverse)
	# The output file name is specified using the stats option (stats=Seal_Stats/$set)
	echo seal.sh overwrite=t in1=$forward in2=$reverse k=20 mm=f t=4 ambiguous=toss hdist=0 trd=t rename fbm int=f ref=ABA_barcodes.fasta stats=Seal_Stats/$set
done > map_it_all.sh

# Run the contents of map_it_all.sh in parallel using 6 threads, and run it in the background
# This will allow the seal.sh script to be run simultaneously on multiple input files, potentially reducing the overall runtime
parallel -j6 'echo {} | sh' :::: map_it_all.sh  &

# Iterate over all files in the Seal_Stats directory that have a ".tsv" extension, corresponding to the output of the seal.sh script
for i in Seal_Stats/*tsv; do
	# Perform an awk command on each file
	awk '$0!~"#"{
		# Set the name variable to the filename, and remove the ".seal.tsv" suffix and "Seal_Stats/" prefix
		name=FILENAME; gsub(".seal.tsv", "", name); gsub("Seal_Stats/", "", name);
		# Print the first and second columns divided by 2, and the name
		# The first column is the barcode sequence and the second column is the count of occurrences of the barcode
		# Dividing the count by 2 is necessary because the counts include both forward and reverse reads
		print $1, $2/2, name;
	}' $i
done |
# Sort the output by the third column in reverse numerical order and the second column in reverse numerical order
sort -k3,3V -k2,2nr |
# Write the output to a gzip-compressed file called
