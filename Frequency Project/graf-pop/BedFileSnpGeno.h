#ifndef BED_FILE_SNP_GENO_H
#define BED_FILE_SNP_GENO_H

#include <fstream>
#include "Util.h"
#include "AncestrySnps.h"
#include "BimFileAncestrySnps.h"
#include "FamFileSamples.h"
#include "SampleGenoDist.h"

static const int BYTE1_IN_BED_FILE = 108;
static const int BYTE2_IN_BED_FILE = 27;
static const int BYTE_OF_SNP_MODE  = 1;

class BedFileSnpGeno
{
public:
    unsigned long baseNums[64];    // bits 1, 10, 100 ... for decoding genos in bed file

    int numAncSnps;
    int numSamples;
    int numBimSnps;
    int numBimAncSnps;

    string bedFile;
    AncestrySnps *ancSnps;
    BimFileAncestrySnps *bimSnps;
    FamFileSamples *famSmps;
    SampleGenoDist *vtxExpGd0;    // Genetic distances from 3 vertices to ref populations when all SNPs have genotypes

public:
    vector<char*> ancSnpSmpGenos; // Genotypes of Ancestry SNPs in an array of chars (0 = AA, 1 = AB; 2 = BB) of chars
    vector<int> ancSnpSnpIds;     // Genotypes of Ancestry SNPs in an array of SNP IDs

    BedFileSnpGeno(string, AncestrySnps*, BimFileAncestrySnps*, FamFileSamples*);
    ~BedFileSnpGeno();
    bool ReadGenotypesFromBedFile();
    void ShowSummary();
    void InitPopPvalues();

private:
    int genoFileLineLen;         // Max length of one sample geno line (with sample info)
    int smpNameLen;

    char GetCompAllele(char);
    int  GetSnpGenoInt(bool, bool);
    char* RecodeBedSnpGeno(char*, int, bool);
};

#endif
