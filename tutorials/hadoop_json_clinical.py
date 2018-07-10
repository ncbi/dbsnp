
#!/opt/python-3.4/bin/python
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
# Script name: hadoop_json_annotation.py 
# Description: a demo script to parse dbSNP RS JSON object and extract clinical 
# rs data.  The script will produce tab-delimited output containing 
# accession_version, allele_id,measure_set_id,organization, accession, snp_id, 
# create_date,update_date,last_evaluated_date,review_status,disease_names,clinical_significances,
# disease_ids_organization,disease_ids_accession, origins, collection_method, citations, and gene_ids.
#
# Sample use:
# python hadoop_json_clinical.py -r hadoop hdfs:///path/to/input -o hdfs:///path/to/output --no-output --jobconf mapreduce.job.name=test.mrjob --jobconf mapreduce.job.reduces=100
#
# Author:  Qiang Wang  wangq2@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------

import os
import json

from mrjob.job import MRJob

        
class MRJsonProcessor(MRJob):

    def mapper(self, _, line):

        data = json.loads(line)

        snp_id=data["refsnp_id"]

        annotations=data["primary_snapshot_data"]["allele_annotations"]
        
        for a in annotations:
            clinical=a["clinical"]

            for c in clinical:
                accession_version=c["accession_version"]
                allele_id=str(c["allele_id"])
                measure_set_id=str(c["measure_set_id"])
                variant_identifiers=c["variant_identifiers"]
                organization=";".join([vi["organization"] for vi in variant_identifiers])
                accession=";".join([vi["accession"] for vi in variant_identifiers])
                snp_id=c["refsnp_id"]
                create_date=c["create_date"]
                update_date=c["update_date"]

                if "last_evaluated_date" in c:
                    last_evaluated_date=c["last_evaluated_date"]
                else:
                    last_evaluated_date=""

                review_status=c["review_status"]
                disease_names=";".join(c["disease_names"])
                clinical_significances=";".join(c["clinical_significances"])
                disease_ids=c["disease_ids"]
                disease_ids_organization=";".join([di["organization"] for di in disease_ids])
                disease_ids_accession=";".join([di["accession"] for di in disease_ids])
                origins=";".join(c["origins"])
                collection_method=";".join(c["collection_method"])
                citations=";".join([str(i) for i in c["citations"]])
                gene_ids=";".join(c["gene_ids"])
                values=[accession_version, allele_id,measure_set_id,organization, accession, str(snp_id), create_date,update_date,last_evaluated_date,review_status,disease_names,clinical_significances,disease_ids_organization,disease_ids_accession,origins, collection_method,citations,gene_ids]
                print("\t".join(values))


if __name__ == '__main__':
    MRJsonProcessor.run()
