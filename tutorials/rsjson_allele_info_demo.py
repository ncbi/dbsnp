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
# produce tab-delimited output containing tthe assembly version, sequence ID,
# position, reference allele, variant allele and ClinVar clinical significance
# if available.
# Author:  Lon Phan  lonphan@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------


import sys
import json
import re

rs = {}

def printAllele_annotations(primary_refsnp):
    '''
    rs clinical significance
    '''
    for annot in primary_refsnp['allele_annotations']:
        for clininfo in annot['clinical']:
            print(",".join(clininfo['clinical_significances']))
            

def getPlacements(info):
    '''
    rs genomic positions
    '''
    rs['alleles'] =  [] # holder for one or more variant alleles
    for alleleinfo in info:
        # has top level placement (ptlp) and assembly info
        if alleleinfo['is_ptlp'] and len(alleleinfo['placement_annot']['seq_id_traits_by_assembly']) > 0:  #get genomic placement and alleles
            assembly_name = alleleinfo['placement_annot'] \
                                      ['seq_id_traits_by_assembly'] \
                                      [0]['assembly_name']

            for a in alleleinfo['alleles']:
                spdi = a['allele']['spdi']
                #print(spdi)
                if spdi['inserted_sequence'] == spdi['deleted_sequence']:
                    rs['alleles'].append({'allele':spdi['deleted_sequence']})
                    rs['seq_id'] = spdi['seq_id']
                    rs['position'] =  spdi['position']
                else: #spdi['inserted_sequence'] != spdi['deleted_sequence']:
                    rs['alleles'].append({'allele':spdi['inserted_sequence']})

def getRefSeqAnnot(info):
    '''
    rs refseq info
    '''
    idx = 0
    for allele in rs['alleles']:

        #print(info[cnt])
        allele_annotation = info[idx]['assembly_annotation'][0]
        if (re.match('^NC_', allele_annotation['seq_id'])): #get only RefSeq annotation on NC
                #print (allele_annotation['seq_id'])
                for g in allele_annotation['genes']:
                    rs['alleles'][idx]['refseq_annot'] = g  #allele and annotation have same ordering
        idx = idx + 1
                       
    


for line in sys.stdin:
    rs_obj = json.loads(line)
    rs['id'] = rs_obj['refsnp_id']
    if 'primary_snapshot_data' in rs_obj:
        getPlacements(rs_obj['primary_snapshot_data']['placements_with_allele'])
        getRefSeqAnnot(rs_obj['primary_snapshot_data']['allele_annotations'])
        idx = 0
        for a in rs['alleles']:
            #print(a)
            rnas = a['refseq_annot']['rnas']
            gene_symbol = a['refseq_annot']['locus']
            gene_name = a['refseq_annot']['name']
            for r in rnas:
                mrna = r['transcript_change']
                protein = r['protein']['variant']['spdi']
                print( "\t".join([rs['id'], a['allele'], gene_name, gene_symbol, mrna['seq_id'], mrna['deleted_sequence'], str(mrna['position']),mrna['deleted_sequence'] ,protein['seq_id'], protein['deleted_sequence'], str(protein['position']),protein['deleted_sequence']]))
 