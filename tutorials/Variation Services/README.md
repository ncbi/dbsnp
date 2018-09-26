# Sample Python script for querying and retrieving dbSNP data from NCBI SPDI service (https://api.ncbi.nlm.nih.gov/variation/v0/)
## # ---Annotate VCF with RS ID and INFO
# python spdi_batch.py -i test_vcf.vcf -t VCF
#
# ---Retrieve RS JSON objects
# python spdi_batch.py -i test_rs.txt -t RS
#
# ---Convert HGVS to SPDI 
# python spdi_batch.py -i test_hgvs.txt -t HGVS
#
# ---Convert HGVS to RS
# python spdi_batch.py -i test_hgvs.txt -t HGVS_RS
##
##
============================
