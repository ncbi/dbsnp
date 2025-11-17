#ifndef SAMPLE_GENO_ANCESTRY_H
#define SAMPLE_GENO_ANCESTRY_H

#include <fstream>
#include "Util.h"
#include "AncestrySnps.h"
#include "FamFileSamples.h"
#include "SampleGenoDist.h"

static const int numSubPopScores = 15;

class GenoSample
{
public:
    string name;
    string father;
    string mother;
    int sex;

    // Ancestry results calculated from genotypes
    int numAncSnps;
    int numHetSnps;
    float hetRate;
    bool ancIsSet;
    float gd1, gd2, gd3;
    float ePct, fPct, aPct;    // Ancestry (EUR, AFR, EAS) components of the sample
    float rawPe, rawPf, rawPa; // Raw scores of above percentages
    int ancGroupId;
    
    float subPopScores[numSubPopScores];

public:
    GenoSample(string);
    void SetAncestryScores(int, float, float, float, float, float, float, float, float, float, bool);
    void SetSubPopGdScores(float*);
    void SetGrafAncGroups();
};


class SampleGenoAncestry
{
private:
    int numSamples;
    int numAncSmps;

    int minAncSnps;
    int totAncSnps;
    int numAncSnps;
    int numThreads;                // Number of threads for parallel computing

    AncestrySnps *ancSnps;
    SampleGenoDist *vtxExpGd0;    // Genetic distances from 3 vertices to ref populations when all SNPs have genotypes

public:
    // Values of the 8 bits, for coding 4 genotypes in one byte. The first geno is stored in the first two bits
    const unsigned int charGenoBitVals[8] = {128, 64, 32, 16, 8, 4, 2, 1};    
  
    vector<GenoSample> samples;
    vector<int> *ancSnpIds;
    vector<unsigned char*> *ancSnpCodedGenos; // Use char, instead of int, to save space

    //                                             0          1        2           3          4        5           6        7         8   
    const string ancSnpFileCols[numSubPops] = {"UKBBEUR", "UKBBAFR", "UKBBEAS", "Nigeria", "Ghana", "Zimbabwe", "Uganda", "Iran", "Barbados", 
                                     
    //        9         10            11        12       13        14            15          16          17                                              
          "China", "Philippines", "Thailand", "Japan", "Nepal", "Pakistan", "Bangladesh", "SriLanka", "India2", 

    //       18       19         20       21        22        23        24        25               
          "Irish", "Finland", "Italy", "TGPPUR", "TGPSAS", "TGPMEX", "France", "Poland"};
          

    const string popScoreNames[numSubPopScores] = {"EA1", "EA2", "EA3", "EA4", "AF1", "AF2", "AF3", "EU1", "EU2", "EU3", "SA1", "SA2", "IC1", "IC2", "IC3"};
    
    // Pair of ref pops to use for each score   
    //                                         EA1 EA2 EA3 EA4 AF1 AF2 AF3 EU1 EU2 EU3 SA1 SA2 IC1 IC2 IC3
    
    const int scorePopIdx1[numSubPopScores] = { 9,  9,  9,  9,  3,  3,  6, 18, 18, 24, 16, 17, 23, 20, 20};
    const int scorePopIdx2[numSubPopScores] = {11, 12, 10, 13,  5,  4,  8, 20, 19, 25, 15, 14, 22, 21,  7};

    // Expected GD scores of subpopulations when all SNPs are included
    float subPopGd0P1[numSubPopScores];
    float subPopGd0P2[numSubPopScores];
    
    // Normalizing GD scores to subPopGd0P1 = -1 and subPopGd0P2 = 1
    const float subPopGdNormP1 = -1.0;
    const float subPopGdNormP2 = 1.0;
    
    SampleGenoAncestry(AncestrySnps*, int=100);
    ~SampleGenoAncestry();

    void SetGenoSamples(const vector<string>&);
    void SetGenoSamples(const vector<FamSample>&);
    int SaveAncestryResults(string, bool=false);
    void SetAncestryPvalues(int);
    void SetSnpGenoData(vector<int>*, vector<unsigned char*>*);
    void SetNumThreads(int);
    void InitPopPvalues();
    void CalculateSubPopGd0Values();

    int GetNumSamples() { return numSamples; };
    int GetNumAncSamples() { return numAncSmps; };
    bool HasEnoughAncestrySnps(int numSnps) { return numSnps >= minAncSnps; }
    void ResetSamples() { samples.clear(); }
    
    void ShowSummary();
};

#endif
