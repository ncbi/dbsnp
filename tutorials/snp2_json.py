import urllib.request
import json
import re

class SnpJsonParser(object):

    api_base = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/'

    def get_json(self, rs):
        
        with urllib.request.urlopen(self.api_base + str(rs)) as response:
            return(json.loads(response.read().decode()))


    def get_ss_info(self, rs_obj):
        
        if 'primary_snapshot_data' in rs_obj:
            print("\t".join(["rs", "handle", "type", "ss or RCV"]))
            self._getSsInfo(rs_obj['refsnp_id'], \
                            rs_obj['primary_snapshot_data']['support'])
            
        
    def _getSsInfo(self, rs, obj):
        for ss in obj:
            id = ss['id']
            print("\t".join([str(rs), ss['submitter_handle'],
                             id['type'], id['value']]))


    def get_allele_info(self, rs_obj):

        rs = {}
        rs['id'] = rs_obj['refsnp_id']
        if 'primary_snapshot_data' in rs_obj:
            self.getPlacements(
                rs_obj['primary_snapshot_data']['placements_with_allele'], rs)

            self.getRefSeqAnnot(
                rs_obj['primary_snapshot_data']['allele_annotations'], rs)

            for a in rs['alleles']:
                if 'refseq_annot' in a:
                    rnas = a['refseq_annot']['rnas']
                    gene_symbol = a['refseq_annot']['locus']
                    gene_name = a['refseq_annot']['name']
                    for r in rnas:
                        if 'codon_aligned_transcript_change' in r:
                            mrna = r['codon_aligned_transcript_change']
                            protein = r['protein']['variant']['spdi']
                            print("\t".join([rs['id'], a['allele'], gene_name,
                                             gene_symbol, mrna['seq_id'],
                                             mrna['deleted_sequence'],
                                             str(mrna['position']),
                                             mrna['deleted_sequence'],
                                             protein['seq_id'],
                                             protein['deleted_sequence'],
                                             str(protein['position']),
                                             protein['deleted_sequence']]))
        
            
    def printAllele_annotations(self, primary_refsnp):
        '''
        rs clinical significance
        '''
        for annot in primary_refsnp['allele_annotations']:
            for clininfo in annot['clinical']:
                print(",".join(clininfo['clinical_significances']))


    def getPlacements(self, info, rs):
        '''
        rs genomic positions
        '''
        rs['alleles'] = []  # holder for one or more variant alleles
        for alleleinfo in info:
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


    def getRefSeqAnnot(self, info, rs):
        '''
        rs refseq info
        '''
        for idx in range(0, len(rs['alleles'])):
            allele_annotation = info[idx]['assembly_annotation'][0]
            # get only RefSeq annotation on NC
            if (re.match('^NC_', allele_annotation['seq_id'])):
                for g in allele_annotation['genes']:
                    # allele and annotation have same ordering
                    rs['alleles'][idx]['refseq_annot'] = g
