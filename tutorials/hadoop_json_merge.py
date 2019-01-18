
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
# Script name: hadoop_json_merge.py
# Description: a demo script to parse dbSNP RS JSON object and obtain the
# records of rs merge history.  The script will produce tab-delimited output
# containing current snp_id, merged snp_id, build_id, and merge date.
#
# Sample use:
# python hadoop_json_merge.py                 \
#     -r hadoop hdfs:///path/to/input         \
#     -o hdfs:///path/to/output               \
#     --no-output                             \
#     --jobconf mapreduce.job.name=test.mrjob \
#     --jobconf mapreduce.job.reduces=100
#
# Author:  Qiang Wang  wangq2@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------


import json

from mrjob.job import MRJob


class MRJsonProcessor(MRJob):

    def mapper(self, _, line):

        data = json.loads(line)

        snp_id = data["refsnp_id"]

        merges = data["dbsnp1_merges"]

        for m in merges:
            values = (snp_id, m["merged_rsid"], m["revision"], m["merge_date"])
            print("\t".join(values))


if __name__ == '__main__':
    MRJsonProcessor.run()
