#include "GrafPop.h"

SampleGenoAncestry *smpGenoAnc = NULL;

int main(int argc, char* argv[])
{
    string usage = "Usage: grafpop <Binary PLINK set or VCF file> <output file>\n";

    string disclaimer =
    "\n *==========================================================================="
    "\n *  GrafPop: Software to infer subject ancestry from genotypes quickly"
    "\n *  Yumi (Jimmy) Jin, PhD"
    "\n *  jinyu@ncbi.nlm.nih.gov"
    "\n *  04/05/2021"
    "\n *"
    "\n *                            PUBLIC DOMAIN NOTICE"
    "\n *               National Center for Biotechnology Information"
    "\n *"
    "\n *  This software/database is a \"United States Government Work\" under the"
    "\n *  terms of the United States Copyright Act.  It was written as part of"
    "\n *  the author's official duties as a United States Government employee and"
    "\n *  thus cannot be copyrighted.  This software/database is freely available"
    "\n *  to the public for use. The National Library of Medicine and the U.S."
    "\n *  Government have not placed any restriction on its use or reproduction."
    "\n *"
    "\n *  Although all reasonable efforts have been taken to ensure the accuracy"
    "\n *  and reliability of the software and data, the NLM and the U.S."
    "\n *  Government do not and cannot warrant the performance or results that"
    "\n *  may be obtained by using this software or data. The NLM and the U.S."
    "\n *  Government disclaim all warranties, express or implied, including"
    "\n *  warranties of performance, merchantability or fitness for any particular"
    "\n *  purpose."
    "\n *"
    "\n *  Please cite the author in any work or product based on this material."
    "\n *"
    "\n *===========================================================================";

    if (argc < 3) {
        cout << disclaimer << "\n\n";
        cout << usage << "\n";
        exit(0);
    }

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    string genoDs, outputFile;
    genoDs = argv[1];
    outputFile = argv[2];

    string fileBase = "";
    GenoDatasetType fileType = CheckGenoDataFile(genoDs, &fileBase);

    if (fileType == GenoDatasetType::NOT_EXISTS) {
        cout << "\nERROR: Genotype file " << genoDs << " doesn't exist!\n\n";
    	return 0;
    }
    else if (fileType == GenoDatasetType::IS_PLINK_GZ) {
        cout << "\nERROR: PLINK set " << genoDs << " is zipped. Please unzip it.\n\n";
        return 0;
    }
    else if (fileType == GenoDatasetType::IS_OTHER) {
        cout << "\nERROR: Genotype file " << genoDs << " should be a binary PLINK set or vcf or vcf.gz file..\n\n";
        return 0;
    }

    string ancSnpFile = FindFile("AncInferSNPs.txt");
    if (ancSnpFile == "") {
        cout << "\nERROR: didn't find file AncInferSNPs.txt. Please put the file under 'data' directory.\n\n";
        return 0;
    }
    AncestrySnps *ancSnps = new AncestrySnps();
    ancSnps->ReadAncestrySnpsFromFile(ancSnpFile);
    //ancSnps->ShowAncestrySnps();

    int totAncSnps = ancSnps->GetNumAncestrySnps();
    int minAncSnps = 100;

    int numThreads = thread::hardware_concurrency();
    numThreads--;

    smpGenoAnc = new SampleGenoAncestry(ancSnps, minAncSnps);

    if (fileType == GenoDatasetType::IS_VCF || fileType == GenoDatasetType::IS_VCF_GZ) {
        VcfSampleAncestrySnpGeno *vcfGeno = new VcfSampleAncestrySnpGeno(genoDs, ancSnps);
        bool dataRead = vcfGeno->ReadDataFromFile();
        if (!dataRead) {
            cout << "\nFailed to read genotype data from " << genoDs << "\n\n";
            return 0;
        }
        vcfGeno->ShowSummary();
        vcfGeno->RecodeSnpGenotypes();

        int numAncSnps = vcfGeno->vcfAncSnpIds.size();
        int numVcfSmps = vcfGeno->GetNumSamples();

        if (smpGenoAnc->HasEnoughAncestrySnps(numAncSnps)) {
            smpGenoAnc->SetGenoSamples(vcfGeno->vcfSamples);
            smpGenoAnc->SetSnpGenoData(&vcfGeno->vcfAncSnpIds, &vcfGeno->vcfAncSnpCodedGenos);
        }
        else {
            cout << "\nWARNING: Ancestry inference not done due to lack of genotyped ancestry SNPs "
             << "(at least " << minAncSnps << " ancestry SNPs are needed).\n\n";
            return 0;
        }
    }
    else if (fileType == GenoDatasetType::IS_PLINK) {
        string bedFile = fileBase + ".bed";
        string bimFile = fileBase + ".bim";
        string famFile = fileBase + ".fam";

        if ( !FileExists(bedFile.c_str()) ||
             !FileExists(bimFile.c_str()) ||
             !FileExists(famFile.c_str())    ) {
            if (!FileExists(bedFile.c_str())) cout << "\nERROR: didn't find " << bedFile << "\n";
            if (!FileExists(bimFile.c_str())) cout << "\nERROR: didn't find " << bimFile << "\n";
            if (!FileExists(famFile.c_str())) cout << "\nERROR: didn't find " << famFile << "\n";
            cout << "\n";
            return 0;
        }

        FamFileSamples *famSmps = new FamFileSamples(famFile);
        famSmps->ShowSummary();

        smpGenoAnc->SetGenoSamples(famSmps->samples);
        int numSmps = smpGenoAnc->GetNumSamples();

        BimFileAncestrySnps *bimSnps = new BimFileAncestrySnps(totAncSnps);
        bimSnps->ReadAncestrySnpsFromFile(bimFile, ancSnps);
        int numBimAncSnps = bimSnps->GetNumBimAncestrySnps();
        bimSnps->ShowSummary();

        if (smpGenoAnc->HasEnoughAncestrySnps(numBimAncSnps)) {
            BedFileSnpGeno *bedGenos = new BedFileSnpGeno(bedFile, ancSnps, bimSnps, famSmps);
            bool hasErr = bedGenos->ReadGenotypesFromBedFile();
            if (hasErr) return 0;
            bedGenos->ShowSummary();

            smpGenoAnc->SetSnpGenoData(&bedGenos->ancSnpSnpIds, &bedGenos->ancSnpSmpGenos);
        }
        else {
            cout << "Ancestry inference not done due to lack of genotyped ancestry SNPs.\n\n";
            return 0;
        }
    }

    cout << "\nLaunching " << numThreads << " threads to calculate ancestry scores.\n";
    smpGenoAnc->SetNumThreads(numThreads);

    mutex iomutex;
    vector<thread> threads(numThreads);

    for (unsigned i = 0; i < numThreads; ++i) {
        threads[i] = thread([&iomutex, i] {
            {
                lock_guard<mutex> iolock(iomutex);
            }

	        smpGenoAnc->SetAncestryPvalues(i);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    smpGenoAnc->SaveAncestryResults(outputFile);

    gettimeofday(&t2, NULL);
    cout << "\n";
    ShowTimeDiff(t1, t2);

    return 1;
}

string GetExecutablePath()
{
    char rawPathName[PATH_MAX];
    realpath(PROC_SELF_EXE, rawPathName);

    string exePath = string(rawPathName);
    size_t slashPos = exePath.find_last_of("/\\");
    string exeDir = exePath.substr(0, slashPos);

    return exeDir;
}

string FindFile(string filename)
{
    string fullFile = filename;

    if (FileExists(fullFile.c_str())) return fullFile;

    string exeDir = GetExecutablePath();
    fullFile = exeDir + "/data/" + filename;
    if (FileExists(fullFile.c_str())) return fullFile;

    fullFile = exeDir + "/" + filename;
    if (FileExists(fullFile.c_str())) return fullFile;

    if(const char* grafPath = getenv("GRAFPATH")) {
        string grafDir = string(grafPath);
        fullFile = grafDir + "/data/" + filename;
        if (FileExists(fullFile.c_str())) return fullFile;

        fullFile = grafDir + "/" + filename;
        if (FileExists(fullFile.c_str())) return fullFile;
    }

    return "";
}
