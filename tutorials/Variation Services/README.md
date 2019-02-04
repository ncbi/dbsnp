## Sample Python script for querying and retrieving dbSNP data from NCBI SPDI service (https://api.ncbi.nlm.nih.gov/variation/v0/)
#### See also [Jupyter notebook](https://github.com/ncbi/dbsnp/tree/master/tutorials/Variation%20Services/Jupyter_Notebook) examples. 
## ================================================
### Annotate VCF with RS ID  
## ================================================
```python spdi_batch.py -i test_vcf.vcf -t VCF```
#### Input
``` 
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
NC_000001.11    10019   .       TA      T       .       .       RS=775809821;dbSNPBuildID=144;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.11    10039   .       A       C       .       .       RS=978760828;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10043   .       T       A       .       .       RS=1008829651;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
```
#### Output
``` 
NC_000001.11    10019   rs775809821     TA      T       .       .       RS=775809821;dbSNPBuildID=144;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.11    10039   rs978760828     A       C       .       .       RS=978760828;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10043   rs1008829651    T       A       .       .       RS=1008829651;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
```
## ================================================
### Convert HGVS to SPDI 
## ================================================
```python spdi_batch.py -i test_hgvs.txt -t HGVS```
#### Input
``` 
NC_000008.10:g.19819724C>G
NC_000008.11:g.19962213C>G
NG_008855.1:g.28143C>G
NM_000237.2:c.1421C>G
```
#### Output
``` 
NC_000008.10:g.19819724C>G      NC_000008.10:19819723:C:G
NC_000008.11:g.19962213C>G      NC_000008.11:19962212:C:G
NG_008855.1:g.28143C>G  NG_008855.1:28142:C:G
NM_000237.2:c.1421C>G   NM_000237.2:1790:C:G
```
## ================================================
### Convert HGVS to RS ID
## ================================================
```python spdi_batch.py -i test_hgvs.txt -t HGVS_RS```
 #### Input
``` 
NC_000008.10:g.19819724C>G
NC_000008.11:g.19962213C>G
NG_008855.1:g.28143C>G
NM_000237.2:c.1421C>G
```
#### Output
``` 
NC_000008.10:g.19819724C>G      NC_000008.10:19819723:C:G       328
NC_000008.11:g.19962213C>G      NC_000008.11:19962212:C:G       328
NG_008855.1:g.28143C>G  NG_008855.1:28142:C:G   328
NM_000237.2:c.1421C>G   NM_000237.2:1790:C:G    328
```
## ================================================
### Retrieve RS JSON objects
## ================================================
```python spdi_batch.py -i test_rs.txt -t RS```
#### Input
``` 
328
rs775809821
rs1052373574
```
#### Output
``` 
{
  "refsnp_id": "328",
  "create_date": "2000-09-19T17:02Z",
  "last_update_date": "2018-05-11T06:02Z",
  "last_update_build_id": "151",
  "dbsnp1_merges": [
    {
      "merged_rsid": "3735962",
      "revision": "108",
      "merge_date": "2002-10-9T00:18Z"
    },
    {
      "merged_rsid": "17482566",
      "revision": "123",
      "merge_date": "2004-10-8T05:17Z"
    },
    {
      "merged_rsid": "52834251",
      "revision": "128",
      "merge_date": "2007-09-21T16:13Z"
    }
  ],
  "citations": [
    1731801,
    1907278,
    2216713,
    16642433,
    16700901,
    17157861,
    17291198,
    17903299,
    18193044,
    18275964,
    18280754,
    18513389,
    18596051,
    18660489,
    18678614,
    18922999,
    19018513,
    19041386,
    19060910,
    19131662,
    19148283,
    19185284,
    19197348,
    19200524,
    19263529,
    19299407,
    19336475,
    19379518,
    19408013,
    19435741,
    19474294,
    19489872,
    19501493,
    19557453,
    19602472,
    19729614,
    19736300,
    19773416,
    19802338,
    19878569,
    19884647,
    19951432,
    20018036,
    20018039,
    20031591,
    20150529,
    20400780,
    20410100,
    20421590,
    20429872,
    20565774,
    20650961,
    21282362,
    21288825,
    21316679,
    21407270,
    21424820,
    21466885,
    21738485,
    21767357,
    21777205,
    21840003,
    21860704,
    21995669,
    22024213,
    22171074,
    22239554,
    22318170,
    22328972,
    22425169,
    22879966,
    23101478,
    23236364,
    23497168,
    24319689,
    24886709,
    24959828,
    25176936,
    25205864,
    25430627,
    25500319,
    25579610,
    25587205,
    25626708,
    25671407,
    25788903,
    25864161,
    26101956,
    26370976,
    26446360,
    26483159,
    26786614,
    26791477,
    26820803,
    26971241,
    26975783,
    26999119,
    27277665,
    27386434,
    27535653,
    27716211,
    28115978,
    28143480,
    28167353,
    28293042,
    28577571,
    28623937,
    28639428,
    28685248,
    28705542,
    29303622
  ],
  "lost_obs_movements": [],
  "present_obs_movements": [
    {
      "component_ids": [
        {
          "type": "subsnp",
          "value": "198888197"
        },
        {
          "type": "subsnp",
          "value": "217321427"
        },
        {
          "type": "subsnp",
          "value": "217397633"
        },
        {
          "type": "subsnp",
          "value": "217399174"
        },
        {
          "type": "subsnp",
          "value": "217407233"
        },
        {
          "type": "subsnp",
          "value": "217418296"
        },
        {
          "type": "subsnp",
          "value": "217419027"
        },
        {
          "type": "subsnp",
          "value": "217422429"
        },
        {
          "type": "subsnp",
          "value": "244294501"
        },
        {
          "type": "subsnp",
          "value": "254171515"
        },
        {
          "type": "subsnp",
          "value": "279724045"
        },
        {
          "type": "subsnp",
          "value": "410878568"
        },
        {
          "type": "subsnp",
          "value": "485584695"
        },
        {
          "type": "subsnp",
          "value": "491921994"
        },
        {
          "type": "subsnp",
          "value": "1397520212"
        },
        {
          "type": "subsnp",
          "value": "1594862344"
        },
        {
          "type": "subsnp",
          "value": "2094987020"
        }
      ],
      "observation": {
        "seq_id": "NC_000008.9",
        "position": 19864003,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "frequency",
          "value": "1000Genomes.1:41010104"
        },
        {
          "type": "frequency",
          "value": "ALSPAC.1:22797219"
        },
        {
          "type": "frequency",
          "value": "ExAC.1:9205345"
        },
        {
          "type": "frequency",
          "value": "GnomAD.1:204721173"
        },
        {
          "type": "frequency",
          "value": "GnomAD_exomes.1:6037335"
        },
        {
          "type": "frequency",
          "value": "TWINSUK.1:22797219"
        },
        {
          "type": "subsnp",
          "value": "223585670"
        },
        {
          "type": "subsnp",
          "value": "234352504"
        },
        {
          "type": "subsnp",
          "value": "241227319"
        },
        {
          "type": "subsnp",
          "value": "342253806"
        },
        {
          "type": "subsnp",
          "value": "484193264"
        },
        {
          "type": "subsnp",
          "value": "490960926"
        },
        {
          "type": "subsnp",
          "value": "491410902"
        },
        {
          "type": "subsnp",
          "value": "536381272"
        },
        {
          "type": "subsnp",
          "value": "655035593"
        },
        {
          "type": "subsnp",
          "value": "779528090"
        },
        {
          "type": "subsnp",
          "value": "780867941"
        },
        {
          "type": "subsnp",
          "value": "782542564"
        },
        {
          "type": "subsnp",
          "value": "783552875"
        },
        {
          "type": "subsnp",
          "value": "834998628"
        },
        {
          "type": "subsnp",
          "value": "985272683"
        },
        {
          "type": "subsnp",
          "value": "1067495952"
        },
        {
          "type": "subsnp",
          "value": "1075340057"
        },
        {
          "type": "subsnp",
          "value": "1328915353"
        },
        {
          "type": "subsnp",
          "value": "1582593788"
        },
        {
          "type": "subsnp",
          "value": "1584057350"
        },
        {
          "type": "subsnp",
          "value": "1620133815"
        },
        {
          "type": "subsnp",
          "value": "1663127848"
        },
        {
          "type": "subsnp",
          "value": "1689111743"
        },
        {
          "type": "subsnp",
          "value": "1711194717"
        },
        {
          "type": "subsnp",
          "value": "1752723251"
        },
        {
          "type": "subsnp",
          "value": "1917826331"
        },
        {
          "type": "subsnp",
          "value": "1928562440"
        },
        {
          "type": "subsnp",
          "value": "1946231552"
        },
        {
          "type": "subsnp",
          "value": "1959093919"
        },
        {
          "type": "subsnp",
          "value": "2024980591"
        },
        {
          "type": "subsnp",
          "value": "2095209247"
        },
        {
          "type": "subsnp",
          "value": "2153202052"
        },
        {
          "type": "subsnp",
          "value": "2470946277"
        },
        {
          "type": "subsnp",
          "value": "2634720469"
        },
        {
          "type": "subsnp",
          "value": "2634720470"
        },
        {
          "type": "subsnp",
          "value": "2634720471"
        },
        {
          "type": "subsnp",
          "value": "2708962560"
        },
        {
          "type": "subsnp",
          "value": "2737022600"
        },
        {
          "type": "subsnp",
          "value": "2748007794"
        },
        {
          "type": "subsnp",
          "value": "2864093419"
        },
        {
          "type": "subsnp",
          "value": "2985433067"
        },
        {
          "type": "subsnp",
          "value": "2986076219"
        },
        {
          "type": "subsnp",
          "value": "3002804512"
        },
        {
          "type": "subsnp",
          "value": "3022826115"
        },
        {
          "type": "subsnp",
          "value": "3348082059"
        },
        {
          "type": "subsnp",
          "value": "3625947311"
        },
        {
          "type": "subsnp",
          "value": "3630013668"
        },
        {
          "type": "subsnp",
          "value": "3630013669"
        },
        {
          "type": "subsnp",
          "value": "3632621006"
        },
        {
          "type": "subsnp",
          "value": "3635162189"
        },
        {
          "type": "subsnp",
          "value": "3640869479"
        },
        {
          "type": "subsnp",
          "value": "3644964726"
        }
      ],
      "observation": {
        "seq_id": "NC_000008.10",
        "position": 19819723,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "frequency",
          "value": "TOPMED.2:384683375"
        },
        {
          "type": "subsnp",
          "value": "244238714"
        },
        {
          "type": "subsnp",
          "value": "252841585"
        },
        {
          "type": "subsnp",
          "value": "2301288406"
        },
        {
          "type": "subsnp",
          "value": "3026281130"
        },
        {
          "type": "subsnp",
          "value": "3555884012"
        }
      ],
      "observation": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "subsnp",
          "value": "10467174"
        }
      ],
      "observation": {
        "seq_id": "NT_030737.7",
        "position": 3540947,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "subsnp",
          "value": "329"
        },
        {
          "type": "subsnp",
          "value": "3173350"
        },
        {
          "type": "subsnp",
          "value": "4921960"
        },
        {
          "type": "subsnp",
          "value": "16343000"
        },
        {
          "type": "subsnp",
          "value": "24648907"
        },
        {
          "type": "subsnp",
          "value": "48420139"
        },
        {
          "type": "subsnp",
          "value": "69043156"
        },
        {
          "type": "subsnp",
          "value": "71648660"
        },
        {
          "type": "subsnp",
          "value": "74808885"
        },
        {
          "type": "subsnp",
          "value": "181341878"
        },
        {
          "type": "subsnp",
          "value": "181834342"
        },
        {
          "type": "subsnp",
          "value": "181835906"
        },
        {
          "type": "subsnp",
          "value": "182258758"
        }
      ],
      "observation": {
        "seq_id": "NT_167187.1",
        "position": 7677869,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "C"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "subsnp",
          "value": "198888197"
        },
        {
          "type": "subsnp",
          "value": "217321427"
        },
        {
          "type": "subsnp",
          "value": "217397633"
        },
        {
          "type": "subsnp",
          "value": "217399174"
        },
        {
          "type": "subsnp",
          "value": "217407233"
        },
        {
          "type": "subsnp",
          "value": "217418296"
        },
        {
          "type": "subsnp",
          "value": "217419027"
        },
        {
          "type": "subsnp",
          "value": "217422429"
        },
        {
          "type": "subsnp",
          "value": "244294501"
        },
        {
          "type": "subsnp",
          "value": "254171515"
        },
        {
          "type": "subsnp",
          "value": "279724045"
        },
        {
          "type": "subsnp",
          "value": "410878568"
        },
        {
          "type": "subsnp",
          "value": "485584695"
        },
        {
          "type": "subsnp",
          "value": "491921994"
        },
        {
          "type": "subsnp",
          "value": "1397520212"
        },
        {
          "type": "subsnp",
          "value": "1594862344"
        },
        {
          "type": "subsnp",
          "value": "2094987020"
        }
      ],
      "observation": {
        "seq_id": "NC_000008.9",
        "position": 19864003,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "frequency",
          "value": "1000Genomes.1:41010104"
        },
        {
          "type": "frequency",
          "value": "ALSPAC.1:22797219"
        },
        {
          "type": "frequency",
          "value": "ExAC.1:9205345"
        },
        {
          "type": "frequency",
          "value": "GnomAD.1:204721173"
        },
        {
          "type": "frequency",
          "value": "GnomAD_exomes.1:6037335"
        },
        {
          "type": "frequency",
          "value": "TWINSUK.1:22797219"
        },
        {
          "type": "subsnp",
          "value": "223585670"
        },
        {
          "type": "subsnp",
          "value": "234352504"
        },
        {
          "type": "subsnp",
          "value": "241227319"
        },
        {
          "type": "subsnp",
          "value": "342253806"
        },
        {
          "type": "subsnp",
          "value": "484193264"
        },
        {
          "type": "subsnp",
          "value": "490960926"
        },
        {
          "type": "subsnp",
          "value": "491410902"
        },
        {
          "type": "subsnp",
          "value": "536381272"
        },
        {
          "type": "subsnp",
          "value": "655035593"
        },
        {
          "type": "subsnp",
          "value": "779528090"
        },
        {
          "type": "subsnp",
          "value": "780867941"
        },
        {
          "type": "subsnp",
          "value": "782542564"
        },
        {
          "type": "subsnp",
          "value": "783552875"
        },
        {
          "type": "subsnp",
          "value": "834998628"
        },
        {
          "type": "subsnp",
          "value": "985272683"
        },
        {
          "type": "subsnp",
          "value": "1067495952"
        },
        {
          "type": "subsnp",
          "value": "1075340057"
        },
        {
          "type": "subsnp",
          "value": "1328915353"
        },
        {
          "type": "subsnp",
          "value": "1582593788"
        },
        {
          "type": "subsnp",
          "value": "1584057350"
        },
        {
          "type": "subsnp",
          "value": "1620133815"
        },
        {
          "type": "subsnp",
          "value": "1663127848"
        },
        {
          "type": "subsnp",
          "value": "1689111743"
        },
        {
          "type": "subsnp",
          "value": "1711194717"
        },
        {
          "type": "subsnp",
          "value": "1752723251"
        },
        {
          "type": "subsnp",
          "value": "1917826331"
        },
        {
          "type": "subsnp",
          "value": "1928562440"
        },
        {
          "type": "subsnp",
          "value": "1946231552"
        },
        {
          "type": "subsnp",
          "value": "1959093919"
        },
        {
          "type": "subsnp",
          "value": "2024980591"
        },
        {
          "type": "subsnp",
          "value": "2095209247"
        },
        {
          "type": "subsnp",
          "value": "2153202052"
        },
        {
          "type": "subsnp",
          "value": "2470946277"
        },
        {
          "type": "subsnp",
          "value": "2634720469"
        },
        {
          "type": "subsnp",
          "value": "2634720470"
        },
        {
          "type": "subsnp",
          "value": "2634720471"
        },
        {
          "type": "subsnp",
          "value": "2708962560"
        },
        {
          "type": "subsnp",
          "value": "2737022600"
        },
        {
          "type": "subsnp",
          "value": "2748007794"
        },
        {
          "type": "subsnp",
          "value": "2864093419"
        },
        {
          "type": "subsnp",
          "value": "2985433067"
        },
        {
          "type": "subsnp",
          "value": "2986076219"
        },
        {
          "type": "subsnp",
          "value": "3002804512"
        },
        {
          "type": "subsnp",
          "value": "3022826115"
        },
        {
          "type": "subsnp",
          "value": "3348082059"
        },
        {
          "type": "subsnp",
          "value": "3625947311"
        },
        {
          "type": "subsnp",
          "value": "3630013668"
        },
        {
          "type": "subsnp",
          "value": "3630013669"
        },
        {
          "type": "subsnp",
          "value": "3632621006"
        },
        {
          "type": "subsnp",
          "value": "3635162189"
        },
        {
          "type": "subsnp",
          "value": "3640869479"
        },
        {
          "type": "subsnp",
          "value": "3644964726"
        }
      ],
      "observation": {
        "seq_id": "NC_000008.10",
        "position": 19819723,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "clinvar",
          "value": "RCV000001598.2"
        },
        {
          "type": "clinvar",
          "value": "RCV000385586.1"
        },
        {
          "type": "frequency",
          "value": "TOPMED.2:384683375"
        },
        {
          "type": "subsnp",
          "value": "244238714"
        },
        {
          "type": "subsnp",
          "value": "252841585"
        },
        {
          "type": "subsnp",
          "value": "2301288406"
        },
        {
          "type": "subsnp",
          "value": "3026281130"
        },
        {
          "type": "subsnp",
          "value": "3555884012"
        }
      ],
      "observation": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "subsnp",
          "value": "10467174"
        }
      ],
      "observation": {
        "seq_id": "NT_030737.7",
        "position": 3540947,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    },
    {
      "component_ids": [
        {
          "type": "subsnp",
          "value": "329"
        },
        {
          "type": "subsnp",
          "value": "3173350"
        },
        {
          "type": "subsnp",
          "value": "4921960"
        },
        {
          "type": "subsnp",
          "value": "16343000"
        },
        {
          "type": "subsnp",
          "value": "24648907"
        },
        {
          "type": "subsnp",
          "value": "48420139"
        },
        {
          "type": "subsnp",
          "value": "69043156"
        },
        {
          "type": "subsnp",
          "value": "71648660"
        },
        {
          "type": "subsnp",
          "value": "74808885"
        },
        {
          "type": "subsnp",
          "value": "181341878"
        },
        {
          "type": "subsnp",
          "value": "181834342"
        },
        {
          "type": "subsnp",
          "value": "181835906"
        },
        {
          "type": "subsnp",
          "value": "182258758"
        }
      ],
      "observation": {
        "seq_id": "NT_167187.1",
        "position": 7677869,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "allele_in_cur_release": {
        "seq_id": "NC_000008.11",
        "position": 19962212,
        "deleted_sequence": "C",
        "inserted_sequence": "G"
      },
      "rsids_in_cur_release": [],
      "last_added_to_this_rs": "151"
    }
  ],
  "primary_snapshot_data": {
    "placements_with_allele": [
      {
        "seq_id": "NC_000008.11",
        "is_ptlp": true,
        "placement_annot": {
          "seq_type": "refseq_chromosome",
          "mol_type": "genomic",
          "seq_id_traits_by_assembly": [
            {
              "assembly_name": "GRCh38.p7",
              "assembly_accession": "GCF_000001405.33",
              "is_top_level": true,
              "is_alt": false,
              "is_patch": false,
              "is_chromosome": true
            }
          ],
          "is_aln_opposite_orientation": false,
          "is_mismatch": false
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NC_000008.11",
                "position": 19962212,
                "deleted_sequence": "C",
                "inserted_sequence": "C"
              }
            },
            "hgvs": "NC_000008.11:g.19962213C="
          },
          {
            "allele": {
              "spdi": {
                "seq_id": "NC_000008.11",
                "position": 19962212,
                "deleted_sequence": "C",
                "inserted_sequence": "G"
              }
            },
            "hgvs": "NC_000008.11:g.19962213C>G"
          }
        ]
      },
      {
        "seq_id": "NC_000008.10",
        "is_ptlp": false,
        "placement_annot": {
          "seq_type": "refseq_chromosome",
          "mol_type": "genomic",
          "seq_id_traits_by_assembly": [
            {
              "assembly_name": "GRCh37.p13",
              "assembly_accession": "GCF_000001405.25",
              "is_top_level": true,
              "is_alt": false,
              "is_patch": false,
              "is_chromosome": true
            }
          ],
          "is_aln_opposite_orientation": false,
          "is_mismatch": false
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NC_000008.10",
                "position": 19819723,
                "deleted_sequence": "C",
                "inserted_sequence": "C"
              }
            },
            "hgvs": "NC_000008.10:g.19819724C="
          },
          {
            "allele": {
              "spdi": {
                "seq_id": "NC_000008.10",
                "position": 19819723,
                "deleted_sequence": "C",
                "inserted_sequence": "G"
              }
            },
            "hgvs": "NC_000008.10:g.19819724C>G"
          }
        ]
      },
      {
        "seq_id": "NG_008855.1",
        "is_ptlp": false,
        "placement_annot": {
          "seq_type": "refseq_genomic",
          "mol_type": "genomic",
          "seq_id_traits_by_assembly": [],
          "is_aln_opposite_orientation": false,
          "is_mismatch": false
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NG_008855.1",
                "position": 28142,
                "deleted_sequence": "C",
                "inserted_sequence": "C"
              }
            },
            "hgvs": "NG_008855.1:g.28143C="
          },
          {
            "allele": {
              "spdi": {
                "seq_id": "NG_008855.1",
                "position": 28142,
                "deleted_sequence": "C",
                "inserted_sequence": "G"
              }
            },
            "hgvs": "NG_008855.1:g.28143C>G"
          }
        ]
      },
      {
        "seq_id": "NM_000237.2",
        "is_ptlp": false,
        "placement_annot": {
          "seq_type": "refseq_mrna",
          "mol_type": "rna",
          "seq_id_traits_by_assembly": [],
          "is_aln_opposite_orientation": false,
          "is_mismatch": false
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NM_000237.2",
                "position": 1790,
                "deleted_sequence": "C",
                "inserted_sequence": "C"
              }
            },
            "hgvs": "NM_000237.2:c.1421C="
          },
          {
            "allele": {
              "spdi": {
                "seq_id": "NM_000237.2",
                "position": 1790,
                "deleted_sequence": "C",
                "inserted_sequence": "G"
              }
            },
            "hgvs": "NM_000237.2:c.1421C>G"
          }
        ]
      },
      {
        "seq_id": "NP_000228.1",
        "is_ptlp": false,
        "placement_annot": {
          "seq_type": "refseq_prot",
          "mol_type": "protein",
          "seq_id_traits_by_assembly": [],
          "is_aln_opposite_orientation": false,
          "is_mismatch": false
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NP_000228.1",
                "position": 473,
                "deleted_sequence": "S",
                "inserted_sequence": "S"
              }
            },
            "hgvs": "NP_000228.1:p.Ser474="
          },
          {
            "allele": {
              "spdi": {
                "seq_id": "NP_000228.1",
                "position": 473,
                "deleted_sequence": "S",
                "inserted_sequence": "*"
              }
            },
            "hgvs": "NP_000228.1:p.Ser474Ter"
          }
        ]
      }
    ],
    "allele_annotations": [
      {
        "frequency": [
          {
            "study_name": "1000Genomes",
            "study_version": 1,
            "local_row_id": 41010104,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 4545,
            "total_count": 5008,
            "common": false
          },
          {
            "study_name": "ALSPAC",
            "study_version": 1,
            "local_row_id": 22797219,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 3443,
            "total_count": 3854,
            "common": false
          },
          {
            "study_name": "ExAC",
            "study_version": 1,
            "local_row_id": 9205345,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 109942,
            "total_count": 121282,
            "common": false
          },
          {
            "study_name": "GnomAD",
            "study_version": 1,
            "local_row_id": 204721173,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 27364,
            "total_count": 29952,
            "common": false
          },
          {
            "study_name": "GnomAD_exomes",
            "study_version": 1,
            "local_row_id": 6037335,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 195507,
            "total_count": 215148,
            "common": false
          },
          {
            "study_name": "TOPMED",
            "study_version": 2,
            "local_row_id": 384683375,
            "observation": {
              "seq_id": "NC_000008.11",
              "position": 19962212,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 114303,
            "total_count": 125568,
            "common": false
          },
          {
            "study_name": "TWINSUK",
            "study_version": 1,
            "local_row_id": 22797219,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "C"
            },
            "allele_count": 3308,
            "total_count": 3708,
            "common": false
          }
        ],
        "clinical": [],
        "submissions": [
          "329",
          "3173350",
          "4921960",
          "10467174",
          "16343000",
          "24648907",
          "48420139",
          "69043156",
          "71648660",
          "74808885",
          "181341878",
          "181834342",
          "181835906",
          "182258758",
          "198888197",
          "217321427",
          "217397633",
          "217399174",
          "217407233",
          "217418296",
          "217419027",
          "217422429",
          "223585670",
          "234352504",
          "241227319",
          "244238714",
          "244294501",
          "252841585",
          "254171515",
          "279724045",
          "342253806",
          "410878568",
          "484193264",
          "485584695",
          "490960926",
          "491410902",
          "491921994",
          "536381272",
          "655035593",
          "779528090",
          "780867941",
          "782542564",
          "783552875",
          "834998628",
          "985272683",
          "1067495952",
          "1075340057",
          "1328915353",
          "1397520212",
          "1582593788",
          "1584057350",
          "1594862344",
          "1620133815",
          "1663127848",
          "1689111743",
          "1711194717",
          "1752723251",
          "1917826331",
          "1928562440",
          "1946231552",
          "1959093919",
          "2024980591",
          "2094987020",
          "2095209247",
          "2153202052",
          "2301288406",
          "2470946277",
          "2634720469",
          "2634720470",
          "2634720471",
          "2708962560",
          "2737022600",
          "2748007794",
          "2864093419",
          "2985433067",
          "2986076219",
          "3002804512",
          "3022826115",
          "3026281130",
          "3348082059",
          "3555884012",
          "3625947311",
          "3630013668",
          "3630013669",
          "3632621006",
          "3635162189",
          "3640869479",
          "3644964726"
        ],
        "assembly_annotation": [
          {
            "seq_id": "NC_000008.11",
            "annotation_release": "Homo sapiens Annotation Release 108",
            "genes": [
              {
                "name": "lipoprotein lipase",
                "id": 4023,
                "locus": "LPL",
                "is_pseudo": false,
                "orientation": "plus",
                "sequence_ontology": [],
                "rnas": [
                  {
                    "id": "NM_000237.2",
                    "transcript_change": {
                      "seq_id": "NM_000237.2",
                      "position": 1789,
                      "deleted_sequence": "TCA",
                      "inserted_sequence": "TCA"
                    },
                    "sequence_ontology": [
                      {
                        "name": "coding_sequence_variant",
                        "accession": "SO:0001580"
                      }
                    ],
                    "product_id": "NP_000228.1",
                    "protein": {
                      "variant": {
                        "spdi": {
                          "seq_id": "NP_000228.1",
                          "position": 473,
                          "deleted_sequence": "S",
                          "inserted_sequence": "S"
                        }
                      },
                      "sequence_ontology": []
                    }
                  }
                ]
              }
            ]
          }
        ]
      },
      {
        "frequency": [
          {
            "study_name": "1000Genomes",
            "study_version": 1,
            "local_row_id": 41010104,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 463,
            "total_count": 5008,
            "common": false
          },
          {
            "study_name": "ALSPAC",
            "study_version": 1,
            "local_row_id": 22797219,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 411,
            "total_count": 3854,
            "common": false
          },
          {
            "study_name": "ExAC",
            "study_version": 1,
            "local_row_id": 9205345,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 11340,
            "total_count": 121282,
            "common": false
          },
          {
            "study_name": "GnomAD",
            "study_version": 1,
            "local_row_id": 204721173,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 2588,
            "total_count": 29952,
            "common": false
          },
          {
            "study_name": "GnomAD_exomes",
            "study_version": 1,
            "local_row_id": 6037335,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 19641,
            "total_count": 215148,
            "common": false
          },
          {
            "study_name": "TOPMED",
            "study_version": 2,
            "local_row_id": 384683375,
            "observation": {
              "seq_id": "NC_000008.11",
              "position": 19962212,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 11265,
            "total_count": 125568,
            "common": false
          },
          {
            "study_name": "TWINSUK",
            "study_version": 1,
            "local_row_id": 22797219,
            "observation": {
              "seq_id": "NC_000008.10",
              "position": 19819723,
              "deleted_sequence": "C",
              "inserted_sequence": "G"
            },
            "allele_count": 400,
            "total_count": 3708,
            "common": false
          }
        ],
        "clinical": [
          {
            "accession_version": "RCV000001598.2",
            "allele_id": 16573,
            "measure_set_id": 1534,
            "variant_identifiers": [
              {
                "organization": "OMIM",
                "accession": "609708.0014"
              }
            ],
            "refsnp_id": "328",
            "create_date": "2012-08-13T00:00Z",
            "update_date": "2017-12-15T00:00Z",
            "last_evaluated_date": "1992-01-15T00:00Z",
            "review_status": "no_assertion_criteria_provided",
            "disease_names": [
              "LIPOPROTEIN LIPASE POLYMORPHISM"
            ],
            "clinical_significances": [
              "benign"
            ],
            "disease_ids": [],
            "origins": [
              "germline"
            ],
            "collection_method": [
              "literature-only"
            ],
            "citations": [
              1731801,
              1907278,
              21042222
            ],
            "gene_ids": [
              "4023"
            ]
          },
          {
            "accession_version": "RCV000385586.1",
            "allele_id": 16573,
            "measure_set_id": 1534,
            "variant_identifiers": [
              {
                "organization": "OMIM",
                "accession": "609708.0014"
              }
            ],
            "refsnp_id": "328",
            "create_date": "2016-12-5T00:00Z",
            "update_date": "2017-12-15T00:00Z",
            "last_evaluated_date": "2016-06-14T00:00Z",
            "review_status": "criteria_provided_single_submitter",
            "disease_names": [
              "Hyperlipoproteinemia, type I"
            ],
            "clinical_significances": [
              "likely-benign"
            ],
            "disease_ids": [
              {
                "organization": "GeneReviews",
                "accession": "NBK1308"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000030379"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000507679"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000511672"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000511673"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000519367"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000558516"
              },
              {
                "organization": "Genetic Testing Registry (GTR)",
                "accession": "GTR000558518"
              },
              {
                "organization": "MedGen",
                "accession": "C0023817"
              },
              {
                "organization": "Orphanet",
                "accession": "444490"
              },
              {
                "organization": "OMIM",
                "accession": "238600"
              }
            ],
            "origins": [
              "germline"
            ],
            "collection_method": [
              "clinical-testing"
            ],
            "citations": [
              21042222
            ],
            "gene_ids": [
              "4023"
            ]
          }
        ],
        "submissions": [
          "329",
          "3173350",
          "4921960",
          "10467174",
          "16343000",
          "24648907",
          "48420139",
          "69043156",
          "71648660",
          "74808885",
          "181341878",
          "181834342",
          "181835906",
          "182258758",
          "198888197",
          "217321427",
          "217397633",
          "217399174",
          "217407233",
          "217418296",
          "217419027",
          "217422429",
          "223585670",
          "234352504",
          "241227319",
          "244238714",
          "244294501",
          "252841585",
          "254171515",
          "279724045",
          "342253806",
          "410878568",
          "484193264",
          "485584695",
          "490960926",
          "491410902",
          "491921994",
          "536381272",
          "655035593",
          "779528090",
          "780867941",
          "782542564",
          "783552875",
          "834998628",
          "985272683",
          "1067495952",
          "1075340057",
          "1328915353",
          "1397520212",
          "1582593788",
          "1584057350",
          "1594862344",
          "1620133815",
          "1663127848",
          "1689111743",
          "1711194717",
          "1752723251",
          "1917826331",
          "1928562440",
          "1946231552",
          "1959093919",
          "2024980591",
          "2094987020",
          "2095209247",
          "2153202052",
          "2301288406",
          "2470946277",
          "2634720469",
          "2634720470",
          "2634720471",
          "2708962560",
          "2737022600",
          "2748007794",
          "2864093419",
          "2985433067",
          "2986076219",
          "3002804512",
          "3022826115",
          "3026281130",
          "3348082059",
          "3555884012",
          "3625947311",
          "3630013668",
          "3630013669",
          "3632621006",
          "3635162189",
          "3640869479",
          "3644964726"
        ],
        "assembly_annotation": [
          {
            "seq_id": "NC_000008.11",
            "annotation_release": "Homo sapiens Annotation Release 108",
            "genes": [
              {
                "name": "lipoprotein lipase",
                "id": 4023,
                "locus": "LPL",
                "is_pseudo": false,
                "orientation": "plus",
                "sequence_ontology": [],
                "rnas": [
                  {
                    "id": "NM_000237.2",
                    "transcript_change": {
                      "seq_id": "NM_000237.2",
                      "position": 1789,
                      "deleted_sequence": "TCA",
                      "inserted_sequence": "TGA"
                    },
                    "sequence_ontology": [
                      {
                        "name": "coding_sequence_variant",
                        "accession": "SO:0001580"
                      }
                    ],
                    "product_id": "NP_000228.1",
                    "protein": {
                      "variant": {
                        "spdi": {
                          "seq_id": "NP_000228.1",
                          "position": 473,
                          "deleted_sequence": "S",
                          "inserted_sequence": "*"
                        }
                      },
                      "sequence_ontology": [
                        {
                          "name": "stop_gained",
                          "accession": "SO:0001587"
                        }
                      ]
                    }
                  }
                ]
              }
            ]
          }
        ]
      }
    ],
    "support": [
      {
        "id": {
          "type": "subsnp",
          "value": "ss329"
        },
        "revision_added": "36",
        "create_date": "2000-09-19T17:02Z",
        "submitter_handle": "DEBNICK"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3173350"
        },
        "revision_added": "98",
        "create_date": "2001-08-15T12:13Z",
        "submitter_handle": "WIAF-CSNP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss4921960"
        },
        "revision_added": "108",
        "create_date": "2002-08-28T10:40Z",
        "submitter_handle": "YUSUKE"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss10467174"
        },
        "revision_added": "116",
        "create_date": "2003-07-11T22:49Z",
        "submitter_handle": "BCM_SSAHASNP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss16343000"
        },
        "revision_added": "120",
        "create_date": "2004-02-27T09:30Z",
        "submitter_handle": "IMCJ-GDT"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss24648907"
        },
        "revision_added": "123",
        "create_date": "2004-09-20T23:19Z",
        "submitter_handle": "PERLEGEN"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss48420139"
        },
        "revision_added": "126",
        "create_date": "2006-03-13T14:28Z",
        "submitter_handle": "APPLERA_GI"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss69043156"
        },
        "revision_added": "127",
        "create_date": "2007-05-18T14:24Z",
        "submitter_handle": "PERLEGEN"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss71648660"
        },
        "revision_added": "127",
        "create_date": "2007-05-18T14:25Z",
        "submitter_handle": "SI_EXO"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss74808885"
        },
        "revision_added": "128",
        "create_date": "2007-08-16T17:42Z",
        "submitter_handle": "AFFY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss181341878"
        },
        "revision_added": "132",
        "create_date": "2010-07-4T16:01Z",
        "submitter_handle": "PAGE_STUDY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss181834342"
        },
        "revision_added": "132",
        "create_date": "2010-07-4T16:01Z",
        "submitter_handle": "PAGE_STUDY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss181835906"
        },
        "revision_added": "132",
        "create_date": "2010-07-4T16:01Z",
        "submitter_handle": "PAGE_STUDY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss182258758"
        },
        "revision_added": "132",
        "create_date": "2010-07-4T16:01Z",
        "submitter_handle": "PAGE_STUDY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss198888197"
        },
        "revision_added": "132",
        "create_date": "2010-07-4T16:01Z",
        "submitter_handle": "BUSHMAN"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217321427"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217397633"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217399174"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217407233"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217418296"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217419027"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss217422429"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss223585670"
        },
        "revision_added": "132",
        "create_date": "2010-07-14T09:28Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss234352504"
        },
        "revision_added": "132",
        "create_date": "2010-07-15T10:12Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss241227319"
        },
        "revision_added": "132",
        "create_date": "2010-07-15T10:12Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss244238714"
        },
        "revision_added": "132",
        "create_date": "2010-05-27T17:38Z",
        "submitter_handle": "OMICIA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss244294501"
        },
        "revision_added": "132",
        "create_date": "2010-07-4T16:01Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss252841585"
        },
        "revision_added": "132",
        "create_date": "2010-08-10T18:40Z",
        "submitter_handle": "OMIM-CURATED-RECORDS"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss254171515"
        },
        "revision_added": "134",
        "create_date": "2011-05-9T22:13Z",
        "submitter_handle": "BL"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss279724045"
        },
        "revision_added": "137",
        "create_date": "2012-05-4T12:43Z",
        "submitter_handle": "GMI"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss342253806"
        },
        "revision_added": "134",
        "create_date": "2011-05-9T22:13Z",
        "submitter_handle": "NHLBI-ESP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss410878568"
        },
        "revision_added": "135",
        "create_date": "2011-09-17T04:05Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss484193264"
        },
        "revision_added": "137",
        "create_date": "2012-05-4T12:43Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss485584695"
        },
        "revision_added": "137",
        "create_date": "2012-05-4T12:43Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss490960926"
        },
        "revision_added": "137",
        "create_date": "2012-05-4T12:43Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss491410902"
        },
        "revision_added": "137",
        "create_date": "2012-05-4T12:43Z",
        "submitter_handle": "EXOME_CHIP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss491921994"
        },
        "revision_added": "137",
        "create_date": "2012-05-4T12:43Z",
        "submitter_handle": "CLINSEQ_SNP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss536381272"
        },
        "revision_added": "146",
        "create_date": "2015-09-8T16:25Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss655035593"
        },
        "revision_added": "138",
        "create_date": "2013-04-25T23:59Z",
        "submitter_handle": "SSMP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss779528090"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss780867941"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss782542564"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss783552875"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss834998628"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss985272683"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "EVA-GONL"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1067495952"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "JMKIDD_LAB"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1075340057"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "JMKIDD_LAB"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1328915353"
        },
        "revision_added": "142",
        "create_date": "2014-08-21T15:37Z",
        "submitter_handle": "1000GENOMES"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1397520212"
        },
        "revision_added": "146",
        "create_date": "2015-09-8T16:25Z",
        "submitter_handle": "HAMMER_LAB"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1582593788"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_GENOME_DK"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1584057350"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_FINRISK"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1594862344"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_DECODE"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1620133815"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_UK10K_ALSPAC"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1663127848"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_UK10K_TWINSUK"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1689111743"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_EXAC"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1711194717"
        },
        "revision_added": "144",
        "create_date": "2015-04-1T10:43Z",
        "submitter_handle": "EVA_MGP"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1752723251"
        },
        "revision_added": "146",
        "create_date": "2015-09-8T16:25Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1917826331"
        },
        "revision_added": "147",
        "create_date": "2016-02-12T12:05Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1928562440"
        },
        "revision_added": "147",
        "create_date": "2016-02-12T12:05Z",
        "submitter_handle": "WEILL_CORNELL_DGM"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1946231552"
        },
        "revision_added": "147",
        "create_date": "2016-02-12T12:05Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss1959093919"
        },
        "revision_added": "147",
        "create_date": "2016-02-12T12:05Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2024980591"
        },
        "revision_added": "149",
        "create_date": "2016-09-14T10:20Z",
        "submitter_handle": "JJLAB"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2094987020"
        },
        "revision_added": "150",
        "create_date": "2016-12-20T12:49Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2095209247"
        },
        "revision_added": "150",
        "create_date": "2016-12-20T12:49Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2153202052"
        },
        "revision_added": "150",
        "create_date": "2016-12-20T12:49Z",
        "submitter_handle": "USC_VALOUEV"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2301288406"
        },
        "revision_added": "150",
        "create_date": "2016-12-20T12:49Z",
        "submitter_handle": "HUMAN_LONGEVITY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2470946277"
        },
        "revision_added": "150",
        "create_date": "2016-12-20T12:49Z",
        "submitter_handle": "TOPMED"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2634720469"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2634720470"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2634720471"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2708962560"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "GRF"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2737022600"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "GNOMAD"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2748007794"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "GNOMAD"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2864093419"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "GNOMAD"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2985433067"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "AFFY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss2986076219"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "AFFY"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3002804512"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "SWEGEN"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3022826115"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3026281130"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "BIOINF_KMB_FNS_UNIBA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3348082059"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "CSHL"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3555884012"
        },
        "revision_added": "151",
        "create_date": "2017-11-8T11:52Z",
        "submitter_handle": "TOPMED"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3625947311"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3630013668"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3630013669"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3632621006"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3635162189"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3640869479"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "subsnp",
          "value": "ss3644964726"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ILLUMINA"
      },
      {
        "id": {
          "type": "frequency",
          "value": "1000Genomes.1:41010104"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "1000Genomes"
      },
      {
        "id": {
          "type": "frequency",
          "value": "ALSPAC.1:22797219"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ALSPAC"
      },
      {
        "id": {
          "type": "frequency",
          "value": "ExAC.1:9205345"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "ExAC"
      },
      {
        "id": {
          "type": "frequency",
          "value": "GnomAD.1:204721173"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "GnomAD"
      },
      {
        "id": {
          "type": "frequency",
          "value": "GnomAD_exomes.1:6037335"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "GnomAD_exomes"
      },
      {
        "id": {
          "type": "frequency",
          "value": "TOPMED.2:384683375"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "TOPMED"
      },
      {
        "id": {
          "type": "frequency",
          "value": "TWINSUK.1:22797219"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": "TWINSUK"
      },
      {
        "id": {
          "type": "clinvar",
          "value": "RCV000001598.2"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": ""
      },
      {
        "id": {
          "type": "clinvar",
          "value": "RCV000385586.1"
        },
        "revision_added": "151",
        "create_date": "2018-05-11T06:02Z",
        "submitter_handle": ""
      }
    ],
    "anchor": "NC_000008.11:0019962212:1:snv",
    "variant_type": "snv"
  }
}
....and MORE
```
 
