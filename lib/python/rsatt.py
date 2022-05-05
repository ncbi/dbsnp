import json
import re
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../lib/python'))
from navs import *
#get RS attribute (rsat) from JSON 
class rsatt(object): 

    acc_chr = {
        'NC_012920': 'MT',
        'NC_000024': 'Y',
        'NC_000023': 'X',
        'NC_000022': '22',
        'NC_000021': '21',
        'NC_000020': '20',
        'NC_000019': '19',
        'NC_000018': '18',
        'NC_000017': '17',
        'NC_000016': '16',
        'NC_000015': '15',
        'NC_000014': '14',
        'NC_000013': '13',
        'NC_000012': '12',
        'NC_000011': '11',
        'NC_000010': '10',
        'NC_000009': '9',
        'NC_000008': '8',
        'NC_000007': '7',
        'NC_000006': '6',
        'NC_000005': '5',
        'NC_000004': '4',
        'NC_000003': '3',
        'NC_000002': '2',
        'NC_000001': '1',
        }
    

    def pmids(self, rs_obj):

        return(rs_obj['citations'])
    
    
    def ss(self, rs_obj):
        
        if 'primary_snapshot_data' in rs_obj:
            # columns: rs, handle, type, ss_or_RCV
            ss_set = []
            for ss in rs_obj['primary_snapshot_data']['support']:
                id = ss['id']
                ss_set.append([rs_obj['refsnp_id'], ss['submitter_handle'], id['type'], id['value']])

            return(ss_set)


    def gene_allele_annot(self, rs_obj):

        rs = {}
        rs['id'] = rs_obj['refsnp_id']
        allele_info = []
        if 'primary_snapshot_data' in rs_obj:
            self.genomic_placements(rs_obj)

            rsa =self.refseq_annot(rs_obj)

            for a in rsa['alleles']:
                if 'refseq_annot' in a:
                    rnas = a['refseq_annot']['rnas']
                    gene_symbol = a['refseq_annot']['locus']
                    gene_name = a['refseq_annot']['name']
                    for r in rnas:
                        if 'codon_aligned_transcript_change' in r:
                            mrna = r['codon_aligned_transcript_change']
                            protein = r['protein']['variant']['spdi']
                            allele_info.append([rs['id'], a['allele'], gene_name,
                                             gene_symbol, mrna['seq_id'],
                                             mrna['deleted_sequence'],
                                             str(mrna['position']),
                                             mrna['deleted_sequence'],
                                             protein['seq_id'],
                                             protein['deleted_sequence'],
                                             str(protein['position']),
                                             protein['deleted_sequence']])

        return(allele_info)
        
            
     def clinical_significance(self, rs_obj):
        '''
        rs clinical significance
        '''
        allele_annot = []
        primary_refsnp = rs_obj['primary_snapshot_data']
        for annot in primary_refsnp['allele_annotations']:
            for clininfo in annot['clinical']:
                allele_annot.append(clininfo['clinical_significances'])

        return(allele_annot)


    def genomic_placements(self, info):
        '''
        rs genomic positions
        '''
        rs = {}
        rs['alleles'] = []  # holder for one or more variant alleles
        for alleleinfo in info['primary_snapshot_data']['placements_with_allele']:
            # has top level placement (ptlp) and assembly info
            if alleleinfo['is_ptlp'] and \
               len(alleleinfo['placement_annot']['seq_id_traits_by_assembly']) > 0:
                # get genomic placement and alleles
                for a in alleleinfo['alleles']:
                    spdi = a['allele']['spdi']
                    if spdi['inserted_sequence'] == spdi['deleted_sequence']:
                        rs['alleles'].append({'allele': spdi['deleted_sequence']})
                        rs['seq_id'] = spdi['seq_id']
                        rs['position'] = spdi['position']
                    else:
                        # spdi['inserted_sequence'] != spdi['deleted_sequence']:
                        rs['alleles'].append({'allele': spdi['inserted_sequence']})
        return rs

    def refseq_annot(self, rsobj):
        '''
        rs refseq info
        '''
        rs = self.genomic_placements(rsobj)
        info = rsobj['primary_snapshot_data']['allele_annotations']
        for idx in range(0, len(rs['alleles'])):
            allele_annotation = info[idx]['assembly_annotation'][0]
            # get only RefSeq annotation on NC
            if (re.match('^NC_', allele_annotation['seq_id'])):
                for g in allele_annotation['genes']:
                    # allele and annotation have same ordering
                    rs['alleles'][idx]['refseq_annot'] = g
        return rs


    def mafs(self, rs_obj):
        
        mafs = {}
        if 'primary_snapshot_data' in rs_obj:        
            for allele in rs_obj['primary_snapshot_data']['allele_annotations']:
                for freq in allele['frequency']:
                    if freq['study_name'] not in mafs:
                        mafs[freq['study_name']] = {}
                        mafs[freq['study_name']]['an'] = freq['total_count']
                        mafs[freq['study_name']]['ac'] = {}

                    mafs[freq['study_name']]['ac'][freq['observation']['inserted_sequence']] \
                        = freq['allele_count']

            for study, maf in mafs.items():
                sorted_ac = sorted(maf['ac'].items(), key=lambda kv: kv[1])
                idx = 0
                # minor allele is the 2nd least abundant allele
                if len(sorted_ac) > 2:
                    idx = 1
                maf['maf_count'] = sorted_ac[idx][1]
                maf['maf'] = float(sorted_ac[idx][1])/maf['an']
                maf['maf_allele'] = sorted_ac[idx][0]
                maf.pop('ac')

        return(mafs)
            
                    
    def gene_consequence(self, rs_obj):
        
        pass

    
    def variation_type(self, rs_obj):

        if 'primary_snapshot_data' in rs_obj:
            return(rs_obj['primary_snapshot_data']['variant_type'])



    def alleles(self, rs_obj):
        alleles = []
        if 'primary_snapshot_data' in rs_obj:
            ptlp = self.__find_ptlp(rs_obj['primary_snapshot_data']['placements_with_allele'])
            if ptlp is not None:
                for allele in ptlp['alleles']:
                    alleles.append(allele['allele']['spdi']['inserted_sequence'])

        return(alleles)


    def chr_pos(self, rs_obj):

        pos = {}
        if 'primary_snapshot_data' in rs_obj:
            ptlp = self.__find_ptlp(rs_obj['primary_snapshot_data']['placements_with_allele'])
            if ptlp is not None:
                for seq_id_trait in ptlp['placement_annot']['seq_id_traits_by_assembly']:
                    if seq_id_trait['is_top_level'] and seq_id_trait['is_chromosome']:
                        pos['assembly'] = seq_id_trait['assembly_name']
                        pos['assembly_accession'] = seq_id_trait['assembly_accession']
                        pos['pos'] = ptlp['alleles'][0]['allele']['spdi']['position']
                        seq_id = ptlp['alleles'][0]['allele']['spdi']['seq_id']
                        acc = seq_id.split('.')[0]
                        pos['chr'] = seq_id
                        if acc in self.acc_chr:
                            pos['chr'] = self.acc_chr[acc]

        return(pos)
                            

    def __find_ptlp(self, placements):
        for p in placements:
            if p['is_ptlp']:
                return p

        return None    
