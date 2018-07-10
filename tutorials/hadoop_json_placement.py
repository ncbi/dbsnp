
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
# Script name: hadoop_json_placement.py 
# Description: a demo script to parse dbSNP RS JSON object and extract the records 
# of rs placement.  The script will produce tab-delimited output containing 
# snp_id, seq_id, is_ptlp, is_aln_opposite_orientation, is_mismatch, position, deleted_sequence, inserted_sequence, and hgvs.
#
# Sample use:
# python hadoop_json_placement.py -r hadoop hdfs:///path/to/input -o hdfs:///path/to/output --no-output --jobconf mapreduce.job.name=test.mrjob --jobconf mapreduce.job.reduces=100
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


        placements=data["primary_snapshot_data"]["placements_with_allele"]

        for p in placements:
            is_ptlp=str(int(p["is_ptlp"]=="true"))
            is_aln_opposite_orientation=str(int(p["placement_annot"]["is_aln_opposite_orientation"]=="true"))
            is_mismatch=str(int(p["placement_annot"]["is_mismatch"]=="true"))
            for a in p["alleles"]:
                if "spdi" in a["allele"]:
                    spdi=a["allele"]["spdi"]
                    hgvs=a["hgvs"]
                    values=(snp_id, p["seq_id"],is_ptlp, is_aln_opposite_orientation, is_mismatch, str(spdi["position"]), spdi["deleted_sequence"], spdi["inserted_sequence"], hgvs)
                else:
                    frameshift=a["allele"]["frameshift"]
                    hgvs=a["hgvs"]
                    values=(snp_id, p["seq_id"],is_ptlp, is_aln_opposite_orientation, is_mismatch, str(frameshift["position"]), "", "", hgvs)
                print("\t".join(values))


if __name__ == '__main__':
    MRJsonProcessor.run()
