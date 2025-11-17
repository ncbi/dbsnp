#include "BedFileSnpGeno.h"

BedFileSnpGeno::BedFileSnpGeno(string bFile, AncestrySnps *aSnps, BimFileAncestrySnps *bSnps, FamFileSamples *fSmps)
{
    bedFile = bFile;
    ancSnps = aSnps;
    bimSnps = bSnps;
    famSmps = fSmps;

    unsigned long base = 1;
    for (int i = 0; i < 64; i++) {
        baseNums[i] = base;
        base = base << 1;
    }

    numAncSnps = ancSnps->GetNumAncestrySnps();
    numBimSnps = bimSnps->GetNumBimSnps();
    numBimAncSnps = bimSnps->GetNumBimAncestrySnps();
    numFamSamples = famSmps->GetNumFamSamples();
    
    numSamples = numFamSamples;
    startSmpNo = 0;
    
    ancSnpSmpGenos = {};
    ancSnpSnpIds = {};

    vtxExpGd0 = new SampleGenoDist(&aSnps->vtxPopExpGds[0], &aSnps->vtxPopExpGds[1],
    &aSnps->vtxPopExpGds[2], &aSnps->vtxPopExpGds[0]);
    vtxExpGd0->TransformAllDists();
    vtxExpGd0->CalculateBaryCenters();
}

BedFileSnpGeno::~BedFileSnpGeno()
{
    for (int i = 0; i < ancSnpSmpGenos.size(); i++) {
        delete ancSnpSmpGenos[i];
    }
    ancSnpSmpGenos.clear();
    ancSnpSnpIds.clear();
}

bool BedFileSnpGeno::SelectFamSampleIds(int stSmpNo, int numSmps)
{
    assert(stSmpNo % 4 == 0); // Read 4 sample genos each time and save them to one byte
    assert(stSmpNo >= 0);
    if (stSmpNo + numSmps > numFamSamples) numSmps = numFamSamples - stSmpNo;

    startSmpNo = stSmpNo;
    numSamples = numSmps;
    
    return true;
}

char BedFileSnpGeno::GetCompAllele(char a)
{
    char c = '0';
    if      (a == 'A') c = 'T';
    else if (a == 'T') c = 'A';
    else if (a == 'G') c = 'C';
    else if (a == 'C') c = 'G';

    return c;
}

unsigned char* BedFileSnpGeno::RecodeBedSnpGeno(char *snpBedGenos, bool swap)
{
    int numBytes = (numSamples - 1) / 4 + 1;
  
    unsigned char *snpGenos = new unsigned char[numBytes]; // char only takes one byte
    for (int i = 0; i < numBytes; i++) snpGenos[i] = 255;

    int smpNo = 0;
    int byteNo = 0;

    // So that 01 01 11 10 = 1 * 64 + 1 * 16 + 3 * 4 + 2 = 94
    unsigned int charBaseInts[4] = {64, 16, 4, 1};
    
    for (byteNo = 0; byteNo < numBytes; byteNo++) {
        char genoByte = snpBedGenos[byteNo];

        unsigned int chrIntGeno = 0; // One char stores genotypes of 4 samples
        
        for (int byteSmpNo = 0; byteSmpNo < 4; byteSmpNo++) {
            int bit1Pos = byteSmpNo * 2;
            int bit2Pos = bit1Pos + 1;;

            int bit1 = genoByte & baseNums[bit1Pos];
            int bit2 = genoByte & baseNums[bit2Pos];

            unsigned int intGeno = 3;
            if      ( bit1 &&  bit2) intGeno = 2;
            else if (!bit1 &&  bit2) intGeno = 1;
            else if (!bit1 && !bit2) intGeno = 0;

            if (swap) {
                if      (intGeno == 0) intGeno = 2;
                else if (intGeno == 2) intGeno = 0;
            }

            chrIntGeno += charBaseInts[byteSmpNo] * intGeno;
            smpNo++;
        }
        
        snpGenos[byteNo] = chrIntGeno;
    }

    return snpGenos;
}

bool BedFileSnpGeno::ReadGenotypesFromBedFile()
{
    bool hasErr = false;

    char header[2];
    char mode[1];

    ifstream bedFilePtr (bedFile, ios::in | ios::binary);

    int snpTotBytes = (numFamSamples - 1) / 4 + 1; // Bytes read from bed file for each SNP
    int snpNumBytes = (numSamples - 1) / 4 + 1;    // Bytes needed to code sample genos to be analyzed
    
    // Check if not all samples are to be analyzed
    bool notAllSmps = numSamples < numFamSamples ? true : false;
    
    long expFileLen = (long)snpTotBytes * numBimSnps + 3;

    bedFilePtr.seekg (0, bedFilePtr.end);
    long fileLen = bedFilePtr.tellg();
    bedFilePtr.seekg(0, bedFilePtr.beg);

    bedFilePtr.read (header, 2);
    bedFilePtr.read (mode, 1);

    if (header[0] != 108 || header[1] != 27) {
        cout << "ERROR: File " << bedFile << " is not a valid PLINK bed file!\n";
        hasErr = true;
    }
    else if (mode[0] != 1) {
        cout << "ERROR: File " << bedFile << " is not in SNP mode!\n";
        hasErr = true;
    }

    if (hasErr) return hasErr;

    char buff[snpTotBytes];             // Reusable memory to keep the genotypes
    int bimAncSnpNo = 0;
    
    for (int i = 0; i < numBimSnps; i++) {
        bedFilePtr.read (buff, snpTotBytes);
        int ancSnpId = bimSnps->GetAncSnpIdGivenBimSnpPos(i);
        int match = bimSnps->GetAlleleMatchGivenBimSnpPos(i);
        bool swap = match ==  2 || match == -2 ? true : false;

        if (ancSnpId >= 0) {
            char* snpGenoStr = new char[snpNumBytes];
            if (notAllSmps) {
                int startByteNo = startSmpNo / 4;
                for (int k = 0; k < snpNumBytes; k++) snpGenoStr[k] = buff[k+startByteNo];
            }
            else 
                for (int j = 0; j < snpTotBytes; j++) snpGenoStr[j] = buff[j];
            
            ASSERT(bimAncSnpNo < numAncSnps, "bim ancestry SNP ID " << bimAncSnpNo << " not less than " << numAncSnps << "\n");

            unsigned char *snpSmpGeno = RecodeBedSnpGeno(snpGenoStr, swap);

            ancSnpSmpGenos.push_back(snpSmpGeno);
            ancSnpSnpIds.push_back(ancSnpId);
            
            delete snpGenoStr;

            bimAncSnpNo++;
        }
    }

    bedFilePtr.close();
    numBimAncSnps = bimAncSnpNo;

    cout << "Bim file has " << numBimSnps << " SNPs. " << bimAncSnpNo << " SNPs are GrafAnc ancestry SNPs.\n";
    cout << "Read genotypes of " << numSamples << " samples from bed file.\n";
    cout << "\n";
    
    if (fileLen != expFileLen) {
        cout << "\n**************************************** WARNING ****************************************\n";
        cout << "Bed file has " << fileLen << " bytes, but is expected to have " << expFileLen << " bytes based on fam and bim files.\n";
        cout << "Please check the PLINK set to make sure nothing is wrong.\n";
        cout << "*****************************************************************************************\n\n";
    }
    
    return 0;
}


int BedFileSnpGeno::GetSnpGenoInt(bool b1, bool b2)
{
    int g = 3;

    if (!b1 && !b2) {
        g = 0;
    }
    else if (!b1 && b2) {
        g = 1;
    }
    else if (b1 && b2) {
        g = 2;
    }

    return g;
}

void BedFileSnpGeno::ShowSummary()
{
    cout << "\n";
    cout << "Total " << numSamples << " samples\n";
    cout << "Total " << numAncSnps << " Ancestry SNPs\n";
    cout << "Total " << numBimAncSnps << " Ancestry SNPs in bim file\n";
}
