#!/usr/bin/env python3
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../../lib/python"))

from navs import Variation

test_cases = [
    'rs328',
    338,
    "NC_000007.14\t8644051\t.\tC\tG,T\t.\t.\tINFO",
    "NC_000007.14\t8644051\t.\tC\tT\t.\t.\tINFO",
    'NC_000007.14:g.8644051C>G',
    'NC_000007.14:g.8644052C>G',
    'NC_000008.10:19813528:1:G',
    'NC_000008.10:19813529:1:G',
    'NC_000008.10:19813529:1:T',
    'NC_000008.11:19956017:1.G',
]


for tc in test_cases:
    print()
    print('Input: ' + str(tc))
    print('-------------------------------------------')

    v = Variation(tc)
    print("RSID:\n" + "\n".join([str(rsid) for rsid in v.asRsidList()]))
    print()
    #print(v.asJson())
    #print()
    #print(v)
    print("SPDI:\n" + "\n".join(v.asSpdiList()))
    print()
    print("HGVS:\n" + "\n".join(v.asHgvsList()))
    print()
    print("VCF:\n" + "\n".join(v.asVcfList()))
    print()


