Scripts and tutorials for using dbSNP data
============================

### directory layout

    .
    ├── refsnp-sample.json.gz                   # Sample data containing one RefSNP JSON example for rs268 for testing rsjson_demo.py  
    ├── rsjson_demo.py               			# Sample Python script to parse RefSNP (rs) JSON object.   The script 
    |               								produces a tab-delimited output containing the assembly version, sequence ID, 
    |                                             	position, reference allele, variant allele and ClinVar clinical significance, 
    |                                             	if available. NOTE: this script was tested using Python 2.7.12.
    └── README.md
