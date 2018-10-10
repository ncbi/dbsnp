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
# Script name: rsjson_allele_info_demo.py
# Description: a demo script to parse dbSNP RS JSON object.  The script will
# produce tab-delimited output containing rs number, handle (if available),
# type (subsnp|clinvar), and id (dbSNP ss|clinvar rcv)
# Author:  Lon Phan  lonphan@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------


import sys
import json


def getSsInfo(rs, obj):
    for ss in obj:
        id = ss['id']
        print("\t".join([str(rs), ss['submitter_handle'],
                        id['type'], id['value']]))


for line in sys.stdin:
    rs_obj = json.loads(line)
    if 'primary_snapshot_data' in rs_obj:
        print("\t".join(["rs", "handle", "type", "ss or RCV"]))
        getSsInfo(rs_obj['refsnp_id'],
                  rs_obj['primary_snapshot_data']['support'])
