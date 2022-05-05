#ifndef BIM_FILE_ANCESTRY_SNPS_H
#define BIM_FILE_ANCESTRY_SNPS_H

#include "Util.h"
#include "AncestrySnps.h"

class BimFileAncestrySnps
{
    int totAncSnps;

    string filename;
    int numBimSnps;
    int numBimAncSnps;
    int numGoodAncSnps;
    int numDupAncSnps = 0;

    int numRsAncSnps;
    int numPos37Snps;
    int numPos38Snps;

    AncestrySnpType ancSnpType;

    // For each bim SNP, if it is an Ancestry SNP, the SNP ID is saved here.
    // Ancestry SNP ID is 0-based. If a bim SNP is not an Ancestry SNP, the SNP ID is -1
    vector<int> bimSnpAncSnpIds;
    vector<int> bimSnpAlleleMatches; // SNP allele matches: 0 = not match; -1, -2: flip; 2, -2: swap

private:
    char FlipAllele(char);

public:
    BimFileAncestrySnps();
    BimFileAncestrySnps(int);
    ~BimFileAncestrySnps();
    void SetTotalAncestrySnps(int totSnps) { totAncSnps = totSnps; };
    char* RecodeBedSnpGeno(char*, bool);
    int ReadAncestrySnpsFromFile(string, AncestrySnps*);
    int CompareAncestrySnpAlleles(const char, const char, const char, const char);
    int GetNumBimSnps() { return numBimSnps; };
    int GetNumBimAncestrySnps() { return numBimAncSnps; };
    int GetAncSnpIdGivenBimSnpPos(int bimSnpPos) {
        return bimSnpPos >= 0 && bimSnpPos < numBimSnps ? bimSnpAncSnpIds[bimSnpPos] : -1;
    };
    int GetAlleleMatchGivenBimSnpPos(int bimSnpPos) {
        return bimSnpPos >= 0 && bimSnpPos < numBimSnps ? bimSnpAlleleMatches[bimSnpPos] : 1000;
    };
    void ShowSummary();
};

#endif
