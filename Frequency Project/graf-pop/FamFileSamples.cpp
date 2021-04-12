#include "FamFileSamples.h"

using namespace std;


FamSample::FamSample(string smp, string dad, string mom, int gender)
{
    name = smp;
    father = dad;
    mother = mom;
    sex = gender;
}


FamFileSamples::FamFileSamples(string file)
{
    filename = file;
    numFamSmps = 0;
    numMales = 0;
    numFemales = 0;

    ReadSamplesFromFile();
}

bool FamFileSamples::Summarize()
{
    bool success = false;
    int arraySize = samples.size();

    if (arraySize == numFamSmps) {
        for (int i = 0; i < numFamSmps; i++) {
            FamSample smp = samples[i];
            if (smp.sex == 1) {
                numMales++;
            }
            else if (smp.sex == 2) {
                numFemales++;
            }
        }

        success = true;
    }

    return success;
}

void FamFileSamples::ShowSummary()
{
    bool success = false;
    cout << "Total " << numFamSmps << " samples in fam file " << filename << ".\n";
    cout << "\t" << numMales << " males\n";
    cout << "\t" << numFemales << " females\n";
    cout << "\n";
}

int FamFileSamples::ReadSamplesFromFile()
{
    int numFileSmps = 0;

    ASSERT(FileExists(filename.c_str()), "File " << filename << " does not exist.");

    int lineLen = 300;
    char fpLine[lineLen];

    FILE *ifp = fopen(filename.c_str(), "r");
    ASSERT(ifp, "Couldn't open file " << filename << ".\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int numSmps = 0;
    int smpSex, pheno;
    char famId[80], smpId[80], dadId[80], momId[80];

    while (fgets(fpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        sscanf(fpLine, "%s %s %s %s %d %d", famId, smpId, dadId, momId, &smpSex, &pheno);

        if (!smpSex) smpSex = 0;

        if (smpId) {
            FamSample sample(smpId, dadId, momId, smpSex);
            samples.push_back(sample);
            numFileSmps++;
        }

        lineNo++;
    }
    fclose(ifp);

    numFamSmps = numFileSmps;

    Summarize();

    return numFileSmps;
}

