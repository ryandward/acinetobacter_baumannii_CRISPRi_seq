cd /home/ryandward/R/acinetobacter_baumannii_CRISPRi_seq/Raw_Data

for i in *R1*gz; do
	export forward=$i;
	export reverse=$(echo $forward | sed 's/R1/R2/g');
	export set=$(echo $forward | sed 's/\_.*//g').seal.tsv;
	echo seal.sh overwrite=t in1=$forward in2=$reverse k=20 mm=f t=4 ambiguous=toss hdist=0 trd=t rename fbm int=f ref=ABA_barcodes.fasta stats=Seal_Stats/$set
done > map_it_all.sh

parallel -j6 'echo {} | sh' :::: map_it_all.sh  &

for i in Seal_Stats/*tsv; do awk '$0!~"#"{name=FILENAME; gsub(".seal.tsv", "", name); gsub("Seal_Stats/", "", name); print $1, $2/2, name;}' $i; done |
sort -k3,3V -k2,2nr |
pigz --best > ../all_counts_seal.tsv.gz
