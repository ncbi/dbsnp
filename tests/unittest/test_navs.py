import unittest
import sys
import os
import json

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../lib/python'))
from navs import *

class NvasTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(NvasTestCase, self).__init__(*args, **kwargs)

    
    def test_init_from_rsid(self):

        tester = Variation(328)
        self.assertTrue(tester, 'initiation from rsid test failed.')
        

    def test_init_from_hgvs(self):
        tester = Variation('NC_000007.14:g.8644051C>G')
        self.assertTrue(tester, 'initiation from HGVS test failed.')

        
    def test_init_from_spdi(self):
        tester = Variation('NC_000008.10:19813528:1:G')
        self.assertTrue(tester, 'initiation from SPDI test failed.')

        
    def test_init_from_vcf(self):
        tester = Variation('NC_000007.14\t8644051\t.\tC\tG,T\t.\t.\tINFO')
        self.assertTrue(tester, 'initiation from VCF test failed.')


    def test_invalid_spdi(self):
        
        with self.assertRaises(ValueError) as cm:
            Variation('NC_000008.11:19956017:1.G')
        
        
    def test_unknown_variation(self):

        # spdi NC_000008.10:19813529:1:G doesn't have a known rs
        tester = Variation('NC_000008.10:19813529:1:G')
        self.assertFalse(tester.asJson(), 'unknown variation test failed.')


    def test_json(self):
        tester = Variation(328)
        parsed = json.loads(tester.asJson())
        expected = '328'
        self.assertEqual(parsed['refsnp_id'], expected, 'retrieving json test failed.')
        
        
    def test_hgvs(self):
        tester = Variation(328)
        parsed = set(tester.asHgvsList())
        expected = set(['NC_000008.11:g.19962213C>G'])
        self.assertEqual(parsed, expected, 'retrieving hgvs test failed.')
        

    def test_spdi(self):

        tester = Variation(328)
        parsed = set(tester.asSpdiList())
        expected = set(['NC_000008.11:19962212:C:G'])
        self.assertEqual(parsed, expected, 'retrieving spdi test failed.')


    def test_vcf(self):
        tester = Variation(328)
        parsed = set(tester.asVcfList())
        expected = set(['NC_000008.11	19962213	rs328	C	G	.	.	.'])
        self.assertEqual(parsed, expected, 'retrieving hgvs test failed.')
    
            
if __name__ == '__main__':

    unittest.main()
                
