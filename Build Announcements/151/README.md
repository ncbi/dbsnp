
# ANNOUNCEMENT: dbSNP Human Build 151 Double in Size to 660 Million Reference SNP (rs)


## April 24, 2018


 
## REMINDER: Important dbSNP changes and notifications 
https://www.ncbi.nlm.nih.gov/mailman/pipermail/dbsnp-announce/2018q2/000186.html
 

dbSNP human build 151 for both GRCh38.p7 and GRCh37.p13 assemblies is now available.  This build include new submissions from TopMed (https://www.nhlbi.nih.gov/research/resources/nhlbi-precision-medicine-initiative/topmed) and GnomAD (http://gnomad.broadinstitute.org/about) which more than double the number of available dbSNP reference SNP (rs) from 324 million to 660 million.   Allele frequency data are available for more than 500 million rs (see summary below) with most being rare (MAF < =0.001).

### Build Summary:

|dbSNP ID|Build Total|
|-------------------------------------|-------------|
|Total Submitted SNP (ss) - redundant |1,803,358,848|
|Total Reference SNP (rs) - non-redundant|660,773,127|

|Genomic mapping|GRCh37.p13|GRCh38.p7|
|--------|-----------|-----------|
|Assembly|648,992,551|660,440,048|
 
|Refseq Annotation|GRCh37.p13|GRCh38.p7|
|-----------------|----------|---------|
|Gene ID|30,194|38,811|
|mRNA Accession|106,113|163,679|
|Protein Accession|82,936|115,774|

|Function Class|GRCh37.p13|GRCh38.p7|
|--------------|----------|---------|
|CDS-INDEL|273,510|313,439|
|CDS-SYNON|3,366,422|3,650,883|
|FRAMESHIFT|411,981|424,833|
|INTRON|250,243,824|348,938,103|
|MISSENSE|6,938,964|7,506,129|
|NCRNA|5,669,033|13,602,650|
|NEARGENE-3|3,818,163|6,400,653|
|NEARGENE-5|15,808,046|25,853,318|
|SPLICE-3|95,699|117,225|
|SPLICE-5|111,212|141,498|
|STOP-GAIN|244,372|279,741|
|STOP-LOSS|10,913|13,661|
|UTR-3|6,129,919|8,668,133|
|UTR-5|2,419,529|4,640,915|


### RS Allele Frequency Counts:

|Minor Allele Frequency (MAF)|TOPMED|GnomAD|1000 Genomes|ExAC|GO-ESP|
|-----------------------|-----------|-----------|----------|---------|---------|
|<=0.001|470,535,424|198,960,749|54,686,241|8,667,575|1,527,303|
|>0.001 and <=0.01|20,069,273|15,967,086|15,944,038|276,851|197,816|
|>0.01 and <0.1|11,410,633|7,940,223|7,852,416|99,605|99,694|
|>=0.1|8,225,418|5,799,895|6,365,867|85,773|63,191|
|Total RS with Frequency|510,240,748|228,667,953|84,848,562|9,129,804|1,888,004|

### Entrez Search
	https://www.ncbi.nlm.nih.gov/snp

### FTP

			GRCh37.p13
			[ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13)


			GRCh38.p7
			[ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7)

## NOTES
dbSNP have several projects planned to analyze its content to identify attributes to aid variant discovery, annotation, prioritization, and interpretation for later release.   The proposed analysis includes: 
###	Different reference allele between GRCh37 and GRCh38
###	Variant type (SNV, INDEL, etc.) change due to assembly differences
###	Reference allele is the minor allele on GRCh37 or GRCh38
###	RS has non-reciprocal mapping between GRCh37 and GRCh38
###	Identify suspect variants including paralogous SNV
###	Identify singleton variants (variants observed in 1 sample)
###	Add new genomic feature annotations:	
				o	segmental duplication
				o	pseudogene
				o	promoter
				o	enhancer
				o	CpG island
				o	methylation site
				
We will post progress and the results here and on GitHub (https://github.com/ncbi/dbsnp).    

Please email snp-admin@ncbi.nlm.nih.gov for any questions, suggestions, and comments.

Regards,

dbSNP Production Team
