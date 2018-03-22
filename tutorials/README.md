Scripts and tutorials for using dbSNP data

JSON Alpha release files based on dbSNP build 150 are available on the FTP site (ftp://ftp.ncbi.nih.gov/snp/.redesign/latest_release/JSON). 
============================

### directory layout

    .
    ├── refsnp-sample.json.gz                # Sample data containing one RefSNP JSON example for rs268 for testing rsjson_demo.py  
    ├── rsjson_demo.py                       # Sample Python script to parse RefSNP (rs) JSON object.   The script
    |                                          produces a tab-delimited output containing the assembly version, sequence ID, 
    |                                          position, reference allele, variant allele and ClinVar clinical significance, 
    |                                          if available. NOTE: this script was tested using Python 2.7.12.
    └── README.md
