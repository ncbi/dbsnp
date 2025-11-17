#ifndef VCF_GRAFANC_SNP_GENO_H
#define VCF_GRAFANC_SNP_GENO_H

#include <zlib.h>
#include <errno.h>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "AncestrySnps.h"

class VcfGrafAncSnpGeno
{
private:
    string vcfFile;
    AncestrySnps *ancSnps;
    
    vector<vector<unsigned char>> vcfAncSnpGtVals;
    
    vector<int> vcfAncSnpChrs;      // chr value from The CHROM string
    vector<int> vcfAncSnpPoss;      // pos value from  POS string
    vector<string> vcfAncSnpSnps;   // The ID string
    vector<string> vcfAncSnpRefs;   // The REF string
    vector<string> vcfAncSnpAlts;   // The ALT string
    vector<int> vcfRsIdAncSnpIds;   // Ancestry SNP ID derived using RS ID
    vector<int> vcfGb37AncSnpIds;   // Ancestry SNP ID derived using Build 37 chr + pos
    vector<int> vcfGb38AncSnpIds;   // Ancestry SNP ID derived using Build 38 chr + pos
  
    int totAncSnps;
    int numVcfSamples; // All samples included in the vcf file
    int numSamples;    // Samples to be analyzed
    int numGenoChars;
  
    int totVcfSnps;
    int putativeAncSnps;
    int numRsIdAncSnps;
    int numGb37AncSnps;
    int numGb38AncSnps;
    int numVcfAncSnps;    
  
    AncestrySnpType ancSnpType;
    
    void CompareAncestrySnpAlleles(const string, const string, const char, const char, int*, int*);
    int RecodeGenotypeGivenIntegers(const int, const int, const int, const int);
    
public:
    // Each enotype (per SNP and sample) is coded with number of alts, i.e., 0 = RR, 1 = RA, 2 = AA, 3 = unknown
    // Each row in the array is the genotypes of one SNP, coded as  an array of integers for all samples
    // Ancestry SNP ID of each row is saved in the array of vcfAncSnpIds
    vector<string> vcfSamples;
    vector<int> vcfAncSnpIds;
    vector<unsigned char*> vcfAncSnpCodedGenos; // Use char, instead of int, to save space
    
    VcfGrafAncSnpGeno(string, AncestrySnps*);
    ~VcfGrafAncSnpGeno();
    
    int GetNumVcfSnps() { return totVcfSnps; };
    int GetNumSamples() { return numSamples; };
    int GetTotalVcfSamples() { return numVcfSamples; };
    int GetNumVcfAncestrySnps() { return numVcfAncSnps; };
    bool ReadHeaderFromFile();
    bool ReadDataFromFile();
    bool ReadDataFromFile(int, int, bool=true);
    void RecodeSnpGenotypes();
    
    void ShowSummary();
    void DeleteAncSnpGtValues();
    void DeleteAncSnpCodedGenos();
};


#endif
