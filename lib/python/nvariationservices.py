#!/usr/bin/env python3
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
# Module name: dbsnp
# Description: a module facilitating dbSNP data retrieval
#
import requests
import json
import re
import sys
import urllib
from itertools import islice, chain
from collections import defaultdict

api_rootURL = 'https://api.ncbi.nlm.nih.gov/variation/v0/'
spdi_fields = ['seq_id', 'position', 'deleted_sequence', 'inserted_sequence']
vcf_key_fields = ['chrom', 'pos', 'id', 'ref']

class refsnp(object):
    # Description: perform queries to ncbi spdi service and convert
    # between VCF, HGVS, and dbSNP RSIDs

    def __init__(self, variant=''):
        self.rsid = 0
        self.req = None
        self.rs = dict()
        self.spdi = list()
        self.hgvs = list()
        self.vcf = list()
        m_rsid = re.fullmatch('(rs)?([1-9]\d*)', variant)
        m_spdi = re.fullmatch('([^:]+:){3}[^:]+', variant)
        if m_rsid:
        
            self.rsid = int(m_rsid.group(2))
            url = api_rootURL + 'beta/refsnp/' + str(self.rsid)
            #print(url)
            self.req = self._apiRequest(url)
            self.rs = json.loads(self.req.text)
            
            ptlp = self.__find_ptlp(self.rs['primary_snapshot_data']['placements_with_allele'])
            novar_spdi = dict()
            novar_hgvs = None
            novar_vcf_key = None
            vcf = defaultdict(list)
            for a in ptlp['alleles']:
                s = a['allele']['spdi']
                if a['hgvs'].endswith('='):
                    novar_spdi = ':'.join([str(s[i]) for i in spdi_fields])
                    novar_hgvs = a['hgvs']
                    url = api_rootURL + 'spdi/' + urllib.parse.quote(novar_spdi) + '/vcf_fields'
                    n_vcf = json.loads(self._apiRequest(url).text)['data']
                    n_vcf['id'] = 'rs' + str(self.rsid)
                    novar_vcf_key = "\t".join([str(n_vcf[i]) for i in vcf_key_fields])
                else:
                    spdi = ':'.join([str(s[i]) for i in spdi_fields])
                    self.spdi.append(spdi)
                    self.hgvs.append(a['hgvs'])
                    url = api_rootURL + 'spdi/' + urllib.parse.quote(spdi) + '/vcf_fields'
                    a_vcf = json.loads(self._apiRequest(url).text)['data']
                    a_vcf['id'] = 'rs' + str(self.rsid)
                    vcf_key = "\t".join([str(a_vcf[i]) for i in vcf_key_fields])
                    vcf[vcf_key].append(a_vcf['alt'])
            
            if not self.hgvs and novar_hgvs:
                self.spdi.append(novar_spdi)
                self.hgvs.append(novar_hgvs)
                vcf[novar_vcf_key].append('.')
            
            for vcf_key in vcf:
                self.vcf.append("\t".join([vcf_key, ','.join(vcf[vcf_key])]))

        elif m_spdi:
            spdi = variant
            self.spdi.append(spdi)

            
            
    def __str__(self):
        return json.dumps(self.rs)
    
            
            
    def asSpdi(self):
        return "\n".join(self.spdi)

        
        
    def asHgvs(self):
        return "\n".join(self.hgvs)
    
        
        
    def asRsid(self):
        return str(self.id())

        
        
    def id(self):
        return self.rsid
    
    
    
    def asJson(self):
        return self.req.text
    
    
    
    def asVcf(self):
        return "\n".join(self.vcf)

        

    def _apiRequest(self, url):
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



    def __find_ptlp(self, placements):
        for p in placements:
            if p['is_ptlp']:
                return p
                
        return None
    