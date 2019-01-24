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
# Module name: navs (NCBI API Variation Services)
# Description: dbSNP data retrieval using NCBI Variation Services API
#
# Authors:  Eugene M. Shekhtman
#           Lon Phan  lonphan@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
# ---------------------------------------------------------------------------
import requests
import json
import re
import sys
import urllib
from itertools import islice, chain
from collections import defaultdict, OrderedDict, namedtuple

api_rootURL = 'https://api.ncbi.nlm.nih.gov/variation/v0/'
spdi_fields = ['seq_id', 'position', 'deleted_sequence', 'inserted_sequence']
TVcf = namedtuple('TVcf', ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'])
TVcf.__new__.__defaults__ = ('.',) * len(TVcf._fields)
VcfT1 = namedtuple('VcfT1', ['chrom', 'pos', 'id', 'ref'])
VcfT2 = namedtuple('VcfT2', ['chrom', 'pos', 'ref', 'alt'])
TVcfX = namedtuple('TVcfX', ['qual', 'filter', 'info'])

class Variation:
    # Description: perform queries to ncbi spdi service and convert
    # between VCF, HGVS, and dbSNP RSIDs

    def __init__(self, init_val=''):
        self.rsid = OrderedDict()
        self.req = None
        self.rs = dict()
        self.spdi = OrderedDict()
        self.hgvs = OrderedDict()
        self.vcf = OrderedDict()
        self.vcfx = TVcfX(*['.'] * 3)

        rsid = None
        if isinstance(init_val, int):
            rsid = init_val
        else:
            m_rsid = re.fullmatch('(rs)?(?P<rsid>[1-9]\d*)', init_val)
            if m_rsid:
                rsid = int(m_rsid.group('rsid'))
        
        if rsid:
            self.__init_from_rsid(rsid)
        else:
            m_other = re.fullmatch('(?P<spdi>([^:]+:){3}[^:]+)?(?P<hgvs>[^:]+:[gcmnrp]\.[^:]+)?(?P<vcf>([^\t]+\t){4}[^t].*)?', init_val)
            if m_other:
                if m_other.group('spdi'):
                    self.__init_from_spdi(m_other.group('spdi'))
                elif m_other.group('hgvs'):
                    self.__init_from_hgvs(m_other.group('hgvs'))
                elif m_other.group('vcf'):
                    self.__init_from_vcf(m_other.group('vcf'))
                else:
                    raise ValueError('Variation format not recognized: '+init_val)
            else:
                raise ValueError('Variation format not recognized: '+init_val)
                

                    
    def __init_from_rsid(self, rsid):
        self.rsid[rsid] = 1
        url = api_rootURL + 'beta/refsnp/' + str(rsid)
        self.req = requests.get(url)
        self.rs = json.loads(self.req.text)
        
        ptlp = self.__find_ptlp(self.rs['primary_snapshot_data']['placements_with_allele'])
        novar_spdi = dict()
        novar_hgvs = None
        novar_vcf_key = None
        vcf4rs = defaultdict(list)
        for a in ptlp['alleles']:
            s = a['allele']['spdi']
            if a['hgvs'].endswith('='):
                novar_spdi = ':'.join([str(s[i]) for i in spdi_fields])
                novar_hgvs = a['hgvs']
                url = api_rootURL + 'spdi/' + urllib.parse.quote(novar_spdi) + '/vcf_fields'
                n_vcf = json.loads(requests.get(url).text)['data']                   # {chrom, pos, ref, alt}
                n_vcf['id'] = 'rs' + str(rsid)                                       # {chrom, pos, id, ref, alt}
                novar_vcf_key = VcfT1(**{f: str(n_vcf[f]) for f in VcfT1._fields})   # (chrom, pos, id, ref)
            else:
                spdi = ':'.join([str(s[i]) for i in spdi_fields])
                self.spdi[spdi] = 1
                self.hgvs[a['hgvs']] = 1
                url = api_rootURL + 'spdi/' + urllib.parse.quote(spdi) + '/vcf_fields'
                a_vcf = json.loads(requests.get(url).text)['data']             # {chrom, pos, ref, alt}
                a_vcf['id'] = 'rs' + str(rsid)                                 # {chrom, pos, id, ref, alt}
                vcf_key = VcfT1(**{f: str(a_vcf[f]) for f in VcfT1._fields})   # (chrom, pos, id, ref)
                vcf4rs[vcf_key].append(a_vcf['alt'])                              # {(chrom, pos, id, ref):[alt,]}
        
        if not self.hgvs and novar_hgvs:
            self.spdi[novar_spdi] = 1
            self.hgvs[novar_hgvs] = 1
            vcf4rs[novar_vcf_key].append('.')    # creates pair (chrom, pos, id, ref):['.']
        
        for vcf_key in vcf4rs:
            alt = ','.join(sorted(vcf4rs[vcf_key]))
            vcf_key2 = VcfT2(**{f: getattr(vcf_key, f) for f in ['chrom', 'pos', 'ref']}, alt = alt)
            vcf = TVcf(**vcf_key._asdict(), alt = alt, **self.vcfx._asdict())
            self.vcf[vcf_key2] = vcf



    def __init_from_spdi(self, spdi):
        self.spdi[spdi] = 1
        
        url = api_rootURL + 'spdi/' + urllib.parse.quote(spdi) + '/contextual'
        cona = json.loads(requests.get(url).text)['data']
        cona_spdi = ':'.join([str(cona[i]) for i in spdi_fields])
        self.spdi[cona_spdi] = 1
        
        url = api_rootURL + 'spdi/' + urllib.parse.quote(spdi) + '/canonical_representative'
        cana = json.loads(requests.get(url).text)['data']
        cana_spdi = ':'.join([str(cana[i]) for i in spdi_fields])
        self.spdi[cana_spdi] = 1
        
        for spdi in self.spdi.keys():
            self.hgvs[self._spdi2hgvs(spdi)] = 1
        url = api_rootURL + 'spdi/' + urllib.parse.quote(cona_spdi) + '/rsids'
        rsid_rsp = json.loads(requests.get(url).text)
        if 'data' in rsid_rsp:
            for rs in rsid_rsp['data']['rsids']:
                self.__init_from_rsid(rs)
        else:
            for spdi in self.spdi.keys():
                url = api_rootURL + 'spdi/' + urllib.parse.quote(spdi) + '/vcf_fields'
                vcfd = json.loads(requests.get(url).text)['data']
                if vcfd['alt'] == vcfd['ref']:
                    vcfd['alt'] = '.'
                vcf_key = VcfT2(**{f: str(vcfd[f]) for f in vcfd})
                vcft = TVcf(**vcf_key._asdict(), **self.vcfx._asdict())
                self.vcf[vcf_key] = vcft


    def __init_from_hgvs(self, hgvs):
        self.hgvs[hgvs] = 1
        
        url = api_rootURL + 'hgvs/' + urllib.parse.quote(hgvs) + '/contextuals'
        conas = json.loads(requests.get(url).text)['data']['spdis']
        for a in conas:
            cona_spdi = ':'.join([str(a[i]) for i in spdi_fields])
            self.spdi[cona_spdi] = 1
            self.__init_from_spdi(cona_spdi)



    def __init_from_vcf(self, vcf):
        #self.vcf[vcf] = 1
        vcft = TVcf(*vcf.split("\t"))
        vcf_key = VcfT2(**{f: getattr(vcft, f) for f in VcfT2._fields})
        self.vcf[vcf_key] = vcft
        self.vcfx = TVcfX(**{f: getattr(vcft, f) for f in TVcfX._fields})
        url = api_rootURL + 'vcf/' + '/'.join(urllib.parse.quote(v) for v in list(vcf_key)) + '/contextuals'
        conas_rsp = json.loads(requests.get(url).text)
        if 'data' in conas_rsp:
            conas = conas_rsp['data']['spdis']
            for a in conas:
                cona_spdi = ':'.join([str(a[i]) for i in spdi_fields])
                self.spdi[cona_spdi] = 1
                self.__init_from_spdi(cona_spdi)

            

    def _spdi2hgvs(self, spdi):
        url = api_rootURL + 'spdi/' + urllib.parse.quote(spdi) + '/hgvs'
        return json.loads(requests.get(url).text)['data']['hgvs']



    def __str__(self):
        return json.dumps(self.rs)
    
            
            
    def asSpdiList(self):
        return self.spdi.keys()

        
        
    def asHgvsList(self):
        return self.hgvs.keys()
    
        
        
    def asRsidList(self):
        return self.rsid.keys()

        
        
    def asJson(self):
        return self.req.text
    
    
    
    def asVcfList(self):
        return ["\t".join(list(self.vcf[v])) for v in self.vcf]

        

    def __find_ptlp(self, placements):
        for p in placements:
            if p['is_ptlp']:
                return p
                
        return None
    