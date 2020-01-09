#!/bin/bash -norc
#retrieve RS flanking sequences
#required Entrez EDirect E-utilities to be installed on your computer (https://www.ncbi.nlm.nih.gov/books/NBK179288/) 

flank_length=100       #user specified flank sequence length (ie. 100bp)
efetch -db snp -id 19,268,12516 -format xml | #retrieve one or more dbSNP rs id comma-delimited
xtract -pattern DocumentSummary -element CHRPOS ACC SNP_ID ALLELE -block GENES/GENE_E -element NAME | #extract xml elements
tr -s ':' '\t' |
while read chr pos accn snp_id allele gene
do
    flank_start=$((pos-flank_length))
    flank_end=$((pos+flank_length))

	lft=$(efetch -db nuccore -format fasta -id "$accn" \
	            -seq_start "$flank_start" -seq_stop "$((pos-1))" < /dev/null |
	          grep -v '>' | tr -d '\n')

	rgt=$(efetch -db nuccore -format fasta -id "$accn" \
	            -seq_start "$((pos+1))" -seq_stop "$flank_end" < /dev/null |
	          grep -v '>' | tr -d '\n')

	echo ">gnl|dbSNP|rs=$snp_id|allele=$allele|gene=$gene|chr=$chr|chr_acc=$accn|flank_start=$flank_start|rs_pos=$pos|flank_end=$flank_end"
	echo "$lft"
	echo "$allele"
	echo "$rgt"
	echo "||" #FASTA delimiter

done
