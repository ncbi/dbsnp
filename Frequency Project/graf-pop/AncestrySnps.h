#ifndef ANCESTRY_SNPS_H
#define ANCESTRY_SNPS_H

#include "Util.h"

static const int numAncSnps = 100437;
static const int numRefPops = 5;
static const int numVtxPops = 3;

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
    float vtxPopAfs[numVtxPops]; // E, F, A, or EUR, AFR, EAS
    float refPopAfs[numRefPops]; // EUR, AFA, ASN, LAT, SAS

public:
    AncestrySnp(int, int, int, int, int, char, char, float*, float*);
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
    double vtxExpGenoDists[numVtxPops][numVtxPops][numAncSnps];
    // Vertex genetic distances summed up using all ancestry SNPs
    GenoDist vtxPopExpGds[numVtxPops];

    string refPopNames[numRefPops];

    int ReadAncestrySnpsFromFile(string);
    int FindSnpIdGivenRs(int);
    int FindSnpIdGivenChrPos(int, int, int);
    AncestrySnp GetAncestrySnp(int);
    void SetVertexExpecteGeneticDists();
    int GetNumAncestrySnps() { return snps.size(); };
    void ShowAncestrySnps();
};

#endif
