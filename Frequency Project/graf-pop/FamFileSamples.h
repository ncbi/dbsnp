#ifndef FAM_FILE_SAMPLES_H
#define FAM_FILE_SAMPLES_H

#include "Util.h"

class FamSample
{
public:
    // Info from the fam file
    string name;
    string father;
    string mother;
    int sex;

public:
    FamSample(string, string, string, int);
};

class FamFileSamples
{
    string filename;
    int numFamSmps;
    int numMales;
    int numFemales;

    bool Summarize();

private:
    int ReadSamplesFromFile();

public:
    vector<FamSample> samples;

    FamFileSamples(string);
    int GetNumFamSamples() {return numFamSmps;};
    void ShowSummary();
};


#endif
