Scripts and tutorials for using dbSNP data

dbSNP build release JSON files are available on the FTP site ([ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON](ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON)). 
============================

### directory layout

    .
    ├── Variation Services                   # Tutorial for working with SPDI Variation Service
    ├── eUtils.ipynb                         # Sample dbSNP eUtils query 
    ├── extract_flank.sh                     # Script using eUtils to get rs flanking sequences 
    ├── MafGraph.ipynb                       # eUtils query and MAF parsing and graphing
    ├── hadoop_json_annotation.py            # parse dbSNP RS JSON object and extract the rs annotation using Hadoop
    ├── hadoop_json_clinical.py              # parse dbSNP RS JSON object and extract clinical rs data using Hadoop
    ├── hadoop_json_merge.py                 # parse dbSNP RS JSON object and extract rs merge history using Hadoop
    ├── hadoop_json_placement.py             # parse dbSNP RS JSON object and extract rs mapping information (ie. position)
    ├── refsnp-sample.json.gz                # Sample data containing one RefSNP JSON example for rs268 for testing             rsjson_demo.py  
    ├── rsjson_demo.py                       # Sample Python script to parse RefSNP (rs) JSON object.   The script
    |                                          produces a tab-delimited output containing the assembly version, sequence ID, 
    |                                          position, reference allele, variant allele and ClinVar clinical significance, 
    |                                          if available. NOTE: this script was tested using Python 2.7.12.
    ├── rsjson_allele_info_demo.py           # Extract allele information  position, mrna and protein SPDI reference allele (inserted) and variant (deleted) sequence
    ├── rsjson_getss_info_demo.py            # Extract submission information (ss, local_snp_id, etc.)

    └── README.md
    
## Run and explore notebook interactively on Binder server.  It may take a few minutes for Binder server to start up.

|Notebook|Description|Binder|
|---|---|---|
|eUtils.ipynb|dbSNP eUtils query|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=%2Ftutorials%2FeUtils.ipynb)|
|MafGraph.ipynb|eUtils query and MAF parsing and graphing|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ncbi/dbsnp/master?filepath=%2Ftutorials%2FMafGraph.ipynb)|
