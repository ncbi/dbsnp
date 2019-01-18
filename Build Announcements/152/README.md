## dbSNP build 152 release notes
https://www.ncbi.nlm.nih.gov/mailman/pipermail/dbsnp-announce/2018q4/000193.html
```text
Organism name: Homo sapiens
Taxonomy ID: 9606

1. RefSnp (RS)

Total RS: 683309324

1.1 RS counts by location:
chr1	51344938
chr2	54963502
chr3	44984330
chr4	43257612
chr5	40588494
chr6	37950204
chr7	36319227
chr8	34468936
chr9	28512286
chr10	30298369
chr11	31059211
chr12	30016334
chr13	22168655
chr14	20193282
chr15	18891124
chr16	20788663
chr17	18415589
chr18	17534420
chr19	14060284
chr20	14403548
chr21	8635154
chr22	8971776
chrX	25569389
chrY	1275185
chrM	3164
Alt Only	29553
Patch	7655
Not On	39879
PAR	798063
Unplaced	93094
Unlocalized	184778

NOTE: Assembly term (ALT, PAR, etc.) definitions (https://www.ncbi.nlm.nih.gov/grc/help/definitions/).

1.2 RS counts by type:
Live RS	655379774
Unsupported RS	103739
Withdrawn RS	6854924
Locationless RS	1286
Merged RS	20969601

RS Type Definitions:
Live = RS has location on reference sequence 
Merged = RS merged to existing RS due to improved clustering algorithm or possibly a change to the reference sequence that would result in identical canonical alleles (e.g. updated repeat regions).
Unsupported = No Submitted SNP (SS) matched any of the RS alleles. (Same causes as merging.)
Locationless = An older RS where the location couldn't be obtained from SS and was not available for the build.
Withdrawn = All SS that belong to the RS cluster were withdrawn.

All above RS including non-Live records have history for traceability.

2. SubSnp (SS)

Total SS: 1828331768
Unmapped SS: 60193
```
