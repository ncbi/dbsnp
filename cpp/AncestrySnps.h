#ifndef ANCESTRY_SNPS_H
#define ANCESTRY_SNPS_H

#include "Util.h"

static const int numAllAncSnps = 282424;
static const int numRefPops = 5;
static const int numVtxPops = 3;
static const int numSubPops = 26;
static const int ancSnpFileOthCols = 6;

class AncestrySnp
{
public:
    int snpId;
    int rs;
    int chr;
    int posG37;
    int posG38;
    char ref;
    char alt;
    float vtxPopAfs[numVtxPops];    // E, F, A, or EUR, AFR, EAS
    float refPopAfs[numRefPops];    // EUR, AFA, ASN, LAT, SAS
    float refSubPopAfs[numSubPops]; // Other reference subcontinental populations

public:
    AncestrySnp(int, int, int, int, int, char, char, float*, float*);
  
    void SetRefSubPopAfs(float*);
};

class AncestrySnps
{
    map<int, int> rsToAncSnpId;
    map<long int, int> pos37ToAncSnpId;
    map<long int, int> pos38ToAncSnpId;

public:
    AncestrySnps();
    ~AncestrySnps();
    vector<AncestrySnp> snps;
    // For each SNP, keeps the expected genetic distance from the 3 vertices to the 3 ref population
    double vtxExpGenoDists[numVtxPops][numVtxPops][numAllAncSnps];
    // Vertex genetic distances summed up using all ancestry SNPs
    GenoDist vtxPopExpGds[numVtxPops];

    string refPopNames[numRefPops];
    string refSubPopNames[numSubPops];
    string nomSubPopNames[numSubPops];
    
    int ReadAncestrySnpsFromFile(string);
    int ReadRefSubPopSnpsFromFile(string);
    int ReadNomSubPopSnpsFromFile(string);
    int FindSnpIdGivenRs(int);
    int FindSnpIdGivenChrPos(int, int, int, const char*, const char*);
    AncestrySnp GetAncestrySnp(int);
    void SetVertexExpecteGeneticDists();
    int GetNumAncestrySnps() { return snps.size(); };
    void ShowAncestrySnps();
};

#endif
