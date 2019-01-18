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
# Script name: spdi_batch.py
# Description: a demo script to perform batch query to ncbi spdi service using
# VCF, HGVS, or dbSNP rs IDs
#
# Sample use:
# ---Annotate VCF with RS ID and INFO
# python spdi_batch.py -i test_vcf.vcf -t VCF
#
# ---Retrieve RS JSON objects
# python spdi_batch.py -i test_rs.txt -t RS
#
# ---Convert HGVS to SPDI
# python spdi_batch.py -i test_hgvs.txt -t HGVS
#
# ---Convert HGVS to RS
# python spdi_batch.py -i test_hgvs.txt -t HGVS_RS
#
# Author:  Lon Phan  lonphan@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
# ---------------------------------------------------------------------------
import requests
import json
import argparse
import re
import sys
from itertools import islice, chain


parser = argparse.ArgumentParser(description='batch process SPDI requests')
parser.add_argument(
    '-i', dest='input_file', required=True,
    help='The name of the input file to parse (VCF, HGVS or rs list, etc.)')
parser.add_argument(
    '-t', dest='input_format', required=True,
    help='The input file format (VCF, HGVS, or RS')
api_rootURL = 'https://api.ncbi.nlm.nih.gov/variation/v0/'


def apiRequest(url):
    try:
        r = requests.get(url)
    except requests.exceptions.Timeout:
        # Maybe set up for a retry, or continue in a retry loop
        print("ERROR: Timeout")
    except requests.exceptions.TooManyRedirects:
        # Tell the user their URL was bad and try a different one
        print("ERROR: bad url =" + url)
    except requests.exceptions.RequestException as e:
        # catastrophic error. bail.
        print(e)
        sys.exit(1)
    if (r.status_code == 200):
        return r
    else:
        print("ERROR: status code = " + str(r.status_code))
        return None


def batch(iterable, n):
    i = iter(iterable)
    piece = list(islice(i, n))
    while piece:
        yield piece
        piece = list(islice(i, n))


def batchRS(infile):
    for rs in infile:
        rs = re.sub('rs', '', rs.rstrip())
        if rs.isdigit():
            url = api_rootURL + 'beta/refsnp/' + rs
            print(url)
            req = requests.get(url)
            print(req.text)


def batchHGVS(infile, handler=0):
    for hgvs in infile:
        hgvs = hgvs.rstrip()
        url = api_rootURL + 'hgvs/' + hgvs + '/contextuals'
        req = apiRequest(url)
        if req and handler:
            handler(hgvs, req)
        elif req:
            spdi = req2spdi(req)
            print(hgvs + "\t" + spdi)


def hgvs2rs(hgvs, req):
    spdi = req2spdi(req)
    rslist = spdi2rs(spdi)
    print("\t".join([hgvs, spdi, ",".join(map(str, rslist))]))


def batchHGVS2RS(infile):
    batchHGVS(infile, hgvs2rs)


def spdi2rs(spdi):
    url = api_rootURL + 'spdi/' + spdi + '/rsids'
    req = apiRequest(url)
    if req:
        return json.loads(req.text)['data']['rsids']
    else:
        return ["no rs found"]


def req2spdi(req):
    reqjson = json.loads(req.text)
    spdiobj = reqjson['data']['spdis'][0]
    spdi = ':'.join([
        spdiobj['seq_id'],
        str(spdiobj['position']),
        spdiobj['deleted_sequence'],
        spdiobj['inserted_sequence']])
    return spdi


def batchVCF(infile):
    vcfbatchsize = 1000
    for batchiter in batch(infile, vcfbatchsize):
        rowcount = 0
        rowdata = ''
        for row in batchiter:
            if not row.startswith("#"):
                rowcount += 1
                rowdata += row
        if rowcount > 0:
            req = requests.post(
                api_rootURL + 'vcf/file/set_rsids?assembly=GCF_000001405.25',
                data=rowdata)
            print(req.text)


batchfunctions = {
    'VCF': batchVCF,
    'RS': batchRS,
    'HGVS': batchHGVS,
    'HGVS_RS': batchHGVS2RS}
args = parser.parse_args()
infile = open(args.input_file, "r")
batchfunctions[args.input_format](infile)
