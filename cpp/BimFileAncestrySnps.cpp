#include "BimFileAncestrySnps.h"

BimFileAncestrySnps::BimFileAncestrySnps()
{
    numDupAncSnps = 0;
    filename = "";
    numBimSnps = 0;
}

BimFileAncestrySnps::BimFileAncestrySnps(int totSnps)
{
    totAncSnps = totSnps;
    numDupAncSnps = 0;
    filename = "";
    numBimSnps = 0;
    numBimAncSnps = 0;
    numGoodAncSnps = 0;
    numRsAncSnps = 0;
    numPos37Snps = 0;
    numPos38Snps = 0;
}

BimFileAncestrySnps::~BimFileAncestrySnps()
{
    bimSnpAncSnpIds.clear();
}

char BimFileAncestrySnps::FlipAllele(char allele)
{
    char flipAllele = '0';

    switch(allele) {
        case 'A': flipAllele = 'T'; break;
        case 'T': flipAllele = 'A'; break;
        case 'G': flipAllele = 'C'; break;
        case 'C': flipAllele = 'G'; break;
    }

    return flipAllele;
}

int BimFileAncestrySnps::CompareAncestrySnpAlleles(const char a1, const char a2, const char expA1, const char expA2)
{
    int match = 0;

    char fa1 = FlipAllele(a1);
    char fa2 = FlipAllele(a2);

    if (a1 == expA1 && a2 == expA2) {
        match = 1;
    }
    else if (a1 == expA2 && a2 == expA1) {
        match = 2;  // swap
    }
    else if (fa1 == expA1 && fa2 == expA2) {
        match = -1; // flip
    }
    else if (fa1 == expA2 && fa2 == expA1) {
        match = -2; // swap and flip
    }

    return match;
}

int BimFileAncestrySnps::ReadAncestrySnpsFromFile(string bimFile, AncestrySnps* ancSnps)
{
    ASSERT(FileExists(bimFile.c_str()), "File " << bimFile << " does not exist.");

    filename = bimFile;

    int lineLen = 65536;
    char fpLine[lineLen];

    FILE *ifp = fopen(bimFile.c_str(), "r");
    ASSERT(ifp, "Could not open " << bimFile << "\n");

    bool fileIsValid = true;

    int i;
    int rs, pos;
    char chrStr[128], rsStr[128], cm[64];
    char refStr[32768], altStr[32768]; // In case there are very, very long refs or alts

    // Read the bed file and save lines with potential ancestry SNPs into memory
    vector<int> bimSnpIds;
    vector<int> rsAncSnpIds;
    vector<int> pos37SnpIds;
    vector<int> pos38SnpIds;
    vector<char> refs;
    vector<char> alts;

    int numSaveSnps = 0;
    int numRsAncSnps = 0;
    int numPos37Snps = 0;
    int numPos38Snps = 0;

    numBimSnps = 0;
    while (fgets(fpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        sscanf(fpLine, "%s %s %s %d %s %s", chrStr, rsStr, cm, &pos, refStr, altStr);

        int chr = GetChromosomeFromString(chrStr);
        int rsNum = GetRsNumFromString(rsStr);

        int rsAncSnpId = ancSnps->FindSnpIdGivenRs(rsNum);
        int pos37SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 37, refStr, altStr);
        int pos38SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 38, refStr, altStr);

        if (rsAncSnpId > -1 || pos37SnpId > -1 || pos38SnpId > -1) {
            char ref = 0, alt = 0;
            if (strlen(refStr) == 1) ref = refStr[0];
            if (strlen(altStr) == 1) alt = altStr[0];

            bimSnpIds.push_back(numBimSnps);
            rsAncSnpIds.push_back(rsAncSnpId);
            pos37SnpIds.push_back(pos37SnpId);
            pos38SnpIds.push_back(pos38SnpId);
            refs.push_back(ref);
            alts.push_back(alt);

            if (rsAncSnpId > -1) numRsAncSnps++;
            if (pos37SnpId > -1) numPos37Snps++;
            if (pos38SnpId > -1) numPos38Snps++;

            numSaveSnps++;
        }

        numBimSnps++;
    }

    fclose(ifp);

    // Rs ID, GB37, or GB38, use whichever returns the most ancestry SNPs to find these SNPs
    ancSnpType = AncestrySnpType::RSID;
    int maxBimAncSnps = numRsAncSnps;

    if (numPos37Snps > maxBimAncSnps) {
        ancSnpType = AncestrySnpType::GB37;
        maxBimAncSnps = numPos37Snps;
    }

    if (numPos38Snps > maxBimAncSnps) {
        ancSnpType = AncestrySnpType::GB38;
        maxBimAncSnps = numPos38Snps;
    }

    for (i = 0; i < numBimSnps; i++) {
        bimSnpAncSnpIds.push_back(-1);
        bimSnpAlleleMatches.push_back(0);
    }

    // Avoid adding same SNP more than once
    numDupAncSnps = 0;
    bool ancIdAdded[totAncSnps];
    for (i = 0; i < totAncSnps; i++) ancIdAdded[i] = false;

    numBimAncSnps = 0;
    int numSwaps = 0;
    for (i = 0; i < numSaveSnps; i++) {
        int bimSnpId = bimSnpIds[i];
        char ref = refs[i];
        char alt = alts[i];

        int ancSnpId = -1;
        if (ancSnpType == AncestrySnpType::RSID) {
            ancSnpId = rsAncSnpIds[i];
        }
        else if (ancSnpType == AncestrySnpType::GB37) {
            ancSnpId = pos37SnpIds[i];
        }
        else if (ancSnpType == AncestrySnpType::GB38) {
            ancSnpId = pos38SnpIds[i];
        }

        if (ancSnpId > -1) {
            AncestrySnp ancSnp = ancSnps->GetAncestrySnp(ancSnpId);
            int match = CompareAncestrySnpAlleles(ref, alt, ancSnp.ref, ancSnp.alt);

            // Only save SNPs with expected alleles
            if (match) {
                if (ancIdAdded[ancSnpId]) {
                    numDupAncSnps++;
                }
                else {
                    ancIdAdded[ancSnpId] = true;
                    bimSnpAncSnpIds[bimSnpId] = ancSnpId;
                    bimSnpAlleleMatches[bimSnpId] = match;
                    if (match == -2 || match ==  2) numSwaps++;
                    numGoodAncSnps++;
                }
            }
            numBimAncSnps++;
        }
    }

    return numBimSnps;
}

void BimFileAncestrySnps::ShowSummary()
{
    int numBadAncSnps = numBimAncSnps - numGoodAncSnps;

    string showSnpType = "RS IDs";
    if      (ancSnpType == AncestrySnpType::GB37) showSnpType = "GRCh 37 chromosome positions";
    else if (ancSnpType == AncestrySnpType::GB38) showSnpType = "GRCh 38 chromosome positions";

    cout << "Total " << numBimSnps << " SNPs in bim file. " << numBimAncSnps << " SNPs are ancestry SNPs.\n";
    cout << "\t" << showSnpType << " are used to find ancestry SNPs.\n";
    cout << "\t" << numGoodAncSnps << " SNPs have expected alleles and will be used for ancestry inference.\n";
    if (numDupAncSnps > 0) cout << "\t" << numDupAncSnps << " ancestry SNPs have multiple entries.\n";

    if (numBadAncSnps > 0) {
        cout << "\t" << numBadAncSnps << " ancestry SNPs do not have expected alleles.\n";
    }
    cout << "\n";

    //  int snpId = GetAncSnpIdGivenBimSnpPos(i*100);
    //  int match = GetAlleleMatchGivenBimSnpPos(i*100);
    //  cout << "No. " << i << ": snpID " << snpId << " match " << match << "\n";
    //}
}
