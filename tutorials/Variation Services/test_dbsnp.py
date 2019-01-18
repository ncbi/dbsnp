#!/usr/bin/env python3
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))

from dbsnp import refsnp

variation = 'rs328'
rs = refsnp(variation)
print('Input: ' + variation)
print('-------------------------------------------')
print('RSID: ' + rs.asRsid())
print()
# print(rs.id())
# print()
#print(rs.asJson())
#print(rs)
print("SPDI:\n" + rs.asSpdi())
print()
print("HGVS:\n" + rs.asHgvs())
print()
print("VCF:\n" + rs.asVcf())
print()
print()

variation = '338'
rs = refsnp(variation)
print('Input: ' + variation)
print('-------------------------------------------')
print('RSID: ' + rs.asRsid())
print()
print("SPDI:\n" + rs.asSpdi())
print()
print("HGVS:\n" + rs.asHgvs())
print()
print("VCF:\n" + rs.asVcf())
print()
print()

variation = 'NC_000008.11:19956017:1:G'
rs = refsnp(variation)
print('Input: ' + variation)
print('-------------------------------------------')
print('RSID: ' + rs.asRsid())
print()
print("SPDI:\n" + rs.asSpdi())
print()
print("HGVS:\n" + rs.asHgvs())
print()
print("VCF:\n" + rs.asVcf())
print()
print()


