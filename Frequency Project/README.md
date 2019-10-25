# dbSNP (https://www.ncbi.nlm.nih.gov/snp)
## ****This project is subject to change due to work in progress.  Please follow this space for updates.****

dbSNP aggregate allele frequency data from multiple sources including:

* HapMap
* 1000 Genomes
* GO-ESP
* ExAC
* GnomAD
* TOPMED
* and many others

Example: (https://www.ncbi.nlm.nih.gov/snp/rs328#frequency_tab) 

dbSNP is currently designing new services to allow searching and retrieving frequency data.   Please send your comments and suggestions to snp-admin@ncbi.nlm.nih.gov or submit a request on GitHub (https://github.com/ncbi/dbsnp/issues).

Thank you for your interest.

Regards,

dbSNP Team



---------------------------------------------------------------------------------------------------
# ASHG 2019 Presentation

### Open access to dbGaP new aggregated allele frequency for variant interpretation.

The slides are available on the dbSNP homepage (https://www.ncbi.nlm.nih.gov/snp/)

dbSNP Frequency content: https://ftp.ncbi.nlm.nih.gov/pub/factsheets/CoLabs_dbGaP_Frequency.pdf

## Abstract:

NCBI database of Genotypes and Phenotypes (dbGaP) contains the results of over 1,200 studies investigating the interaction of genotype and phenotype. The database has over two million subjects and hundreds of millions of variants along with thousands of phenotypes and molecular assay data. This unprecedented volume and variety of data promise huge opportunities to identify genetic factors that influence health and disease. With this possibility, NIH has recently updated the Genomic Summary Results (GSR) access restriction to allow responsible sharing and use of the dbGaP GSR data (https://grants.nih.gov/grants/guide/notice-files/NOT-OD-19-023.html).

In fulfilling the updated GSR policy and to improve variant interpretation for health and disease, NCBI has undertaken the challenging task to compute allele frequency for variants in dbGaP across approved un-restricted studies and provide the data as ‘open-access’ to the public. The work involved harmonizing and normalizing heterogeneous data and file formats either from GWAS chip array or direct sequencing. Using dbSNP and dbGaP workflows the data were QA/QC and were transformed to standard VCF format as input into an automated pipeline to aggregate, remap and cluster to existing dbSNP rs, and compute allele frequency. Allele frequencies are calculated for 12 major populations including European, Hispanic, African, Asian, and others that were computed using GRAF-pop (Jin et al., 2019).

The initially released data (pending) included MAF for about 500M sites with data in dbSNP and +20M novel sites from +150 thousand subjects across more than 60 studies. dbGaP MAF data are consistent with MAF data previously reported in GnomAD for the same variants. Moreover, dbGaP has frequency data for novel and existing variants in dbSNP and ClinVar but not reported in 1000Genomes, GnomAD, ExAC, or TopMed. The data volume will grow and can potentially reach over a billion variants from millions of subjects combined across all dbGaP studies. New studies will be added to future dbSNP build release for ‘de novo’ allele frequency calculation across all studies. This presentation will describe the available resources (Web, FTP, and API) and how researchers, clinicians, and developers can incorporate these data into their workflows and applications to understanding human variation and disease.

## Acknowledgments:
Work at NCBI is supported by the NIH Intramural Research Program and the National Library of Medicine.
