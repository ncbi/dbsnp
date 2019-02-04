import unittest
import json
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../tutorials'))
from snp2_json import *
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../tutorials'))

class SnpJsonParserTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(SnpJsonParserTestCase, self).__init__(*args, **kwargs)
        self.a_var = Variation(328)
        self.rs_obj = json.loads(self.a_var.asJson())        
        self.driver = SnpJsonParser()

    def test_get_ss(self):

        ss_set = self.driver.get_ss_info(self.rs_obj)
        self.assertTrue(len(ss_set) == 102, 'parse ss record failed.')


    def test_get_pubmed(self):
        pmid = self.driver.get_pubmed(self.rs_obj)
        self.assertTrue(len(pmid) == 113, 'parse pubmed id failed.')


    def test_get_allele_info(self):
        allele_info = self.driver.get_allele_info(self.rs_obj)
        self.assertTrue(len(allele_info) == 2, 'get allele info test failed.')

        
    def test_get_Allele_annotations(self):
        allele_annot = self.driver.get_Allele_annotations(self.rs_obj['primary_snapshot_data'])
        expected = [['benign'], ['likely-benign']]

        self.assertTrue(allele_annot == expected, 'parse allele info failed.')
        

    def test_getPlacements(self):
        rs = {}
        self.driver.getPlacements(self.rs_obj['primary_snapshot_data']['placements_with_allele'], rs)
        self.assertTrue('alleles' in rs, 'getPlacements test failed.')

        # def test_getRefseqAnnot(self):
        self.driver.getRefSeqAnnot(self.rs_obj['primary_snapshot_data']['allele_annotations'], rs)
        self.assertTrue('alleles' in rs, 'getRefseqAnnot test failed.')


    def test_get_mafs(self):
        mafs = self.driver.get_mafs(self.rs_obj)
        self.assertTrue(mafs['1000Genomes']['maf_count'] == 463, 'get mafs failed.')


    def test_get_variation_type(self):
        self.assertEqual(self.driver.get_variation_type(self.rs_obj), 'snv', 'get variant type failed.')


    def test_get_alleles(self):
        alleles = self.driver.get_alleles(self.rs_obj)
        expected = ['C', 'G']
        self.assertEqual(alleles, expected, 'get alleles failed.')


    def test_get_pos(self):
        pos = self.driver.get_chr_pos(self.rs_obj)
        self.assertEqual(int(pos['pos']), 19962212, 'get position failed.')
        
        
        
if __name__ == '__main__':

    unittest.main()
