#include "AncestrySnps.h"
#include <fstream>

AncestrySnp::AncestrySnp(int id, int rsNum, int ch, int g37, int g38, char a1, char a2, float* refPops, float *vtxPops)
{
    snpId = id;
    rs = rsNum;
    chr = ch;
    posG37 = g37;
    posG38 = g38;
    ref = a1;
    alt = a2;

    for (int i = 0; i < numRefPops; i++) refPopAfs[i] = refPops[i];
    for (int i = 0; i < numVtxPops; i++) vtxPopAfs[i] = vtxPops[i];
}

void AncestrySnp::SetRefSubPopAfs (float* popAfs)
{
    // By default, AFs for normalization is the same as reference AFs
    for (int i = 0; i < numSubPops; i++) {
        refSubPopAfs[i] = popAfs[i];
    }
}

AncestrySnps::AncestrySnps()
{
  
  refPopNames[0] = "European";
  refPopNames[1] = "African";
  refPopNames[2] = "Asian";
  refPopNames[3] = "Mexican";
  refPopNames[4] = "Indian-Pakistani";

  string subPopNames[numSubPops];
  
  subPopNames[0] = "chn"; // Chinease, especially Northern Chinese
  subPopNames[1] = "jpn"; // Japanese
  subPopNames[2] = "sea"; // Southeast Asian: Thai, Vietnamese, ...
  
  for (int popId = 0; popId < numSubPops; popId++) {
      refSubPopNames[popId] = "ref_" + subPopNames[popId];
      nomSubPopNames[popId] = "nom_" + subPopNames[popId];
  }
  
  snps = {};
}


AncestrySnps::~AncestrySnps()
{
    snps.clear();
    rsToAncSnpId.clear();
    pos37ToAncSnpId.clear();
    pos38ToAncSnpId.clear();
}

int AncestrySnps::ReadAncestrySnpsFromFile(string ancSnpFile)
{
    ASSERT(FileExists(ancSnpFile.c_str()), "File " << ancSnpFile << " does not exist!\n");

    double popExpPfSums[numRefPops];
    double popExpPaSums[numRefPops];
    double popExpPeSums[numRefPops];

    for (int popId = 0; popId < numRefPops; popId++) {
        popExpPeSums[popId] = 0;
        popExpPfSums[popId] = 0;
        popExpPaSums[popId] = 0;
    }

    int lineLen = 500;
    char snpLine[lineLen];

    FILE *ifp = fopen(ancSnpFile.c_str(), "r");
    ASSERT(ifp, "Couldn't open file " << ancSnpFile << ".\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int snpId = 0;
    int rsNum, chr, g37, g38;
    char a1, a2;
    float rfEur, rfAfr, rfEas, rfLat, rfSas, vtEur, vtAfr, vtEas;

    while (fgets(snpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        if (lineNo == 0) {
            if (snpLine[0] != 'c' || snpLine[1] != 'h' || snpLine[2] != 'r') {
                fileIsValid = false;
            }
        }
        else {
            sscanf(snpLine, "%d %d %d %d %c %c %f %f %f %f %f",
            &chr, &g37, &g38, &rsNum, &a1, &a2, &rfEur, &rfAfr, &rfEas, &rfLat, &rfSas);

            float* refPopAfs = new float[numRefPops];
            float* vtxPopAfs = new float[numVtxPops];

            refPopAfs[0] = rfEur;
            refPopAfs[1] = rfAfr;
            refPopAfs[2] = rfEas;
            refPopAfs[3] = rfLat;
            refPopAfs[4] = rfSas;

            vtxPopAfs[0] = rfEur;
            vtxPopAfs[1] = rfAfr;
            vtxPopAfs[2] = rfEas;

            AncestrySnp ancSnp(snpId, rsNum, chr, g37, g38, a1, a2, refPopAfs, vtxPopAfs);

            snps.push_back(ancSnp);

            long int chrPos37 = (long)chr * 1000000000 + g37;
            long int chrPos38 = (long)chr * 1000000000 + g38;

            rsToAncSnpId[rsNum] = snpId;
            pos37ToAncSnpId[chrPos37] = snpId;
            pos38ToAncSnpId[chrPos38] = snpId;

            double pev = refPopAfs[0];
            double pfv = refPopAfs[1];
            double pav = refPopAfs[2];

            double qev = 1 - pev;
            double qfv = 1 - pfv;
            double qav = 1 - pav;

            for (int vtxId = 0; vtxId < 3; vtxId++) {
                double pv  = vtxPopAfs[vtxId];
                double qv  = 1 - pv;

                double aaPev = log(pev) * 2;
                double bbPev = log(qev) * 2;
                double abPev = log(pev) + log(qev) + log(2);

                double aaPfv = log(pfv) * 2;
                double bbPfv = log(qfv) * 2;
                double abPfv = log(pfv) + log(qfv) + log(2);

                double aaPav = log(pav) * 2;
                double bbPav = log(qav) * 2;
                double abPav = log(pav) + log(qav) + log(2);

                double eGd = aaPev * pv * pv + bbPev * qv * qv + abPev * 2 * pv * qv;
                double fGd = aaPfv * pv * pv + bbPfv * qv * qv + abPfv * 2 * pv * qv;
                double aGd = aaPav * pv * pv + bbPav * qv * qv + abPav * 2 * pv * qv;

                vtxExpGenoDists[vtxId][0][snpId] = eGd;
                vtxExpGenoDists[vtxId][1][snpId] = fGd;
                vtxExpGenoDists[vtxId][2][snpId] = aGd;

                popExpPeSums[vtxId] += eGd;
                popExpPfSums[vtxId] += fGd;
                popExpPaSums[vtxId] += aGd;
            }

            delete refPopAfs;
            delete vtxPopAfs;

            snpId++;
        }

        lineNo++;
    }
    fclose(ifp);

    ASSERT(snpId == numAllAncSnps, "snpId = " << snpId << " all anc SNPs = " << numAllAncSnps << ".\n");

    for (int vtxId = 0; vtxId < 3; vtxId++) {
        vtxPopExpGds[vtxId].e = -1 * popExpPeSums[vtxId]/snpId;
        vtxPopExpGds[vtxId].f = -1 * popExpPfSums[vtxId]/snpId;
        vtxPopExpGds[vtxId].a = -1 * popExpPaSums[vtxId]/snpId;
    }

    if (0) {
        cout << "Expected vertex genetic distances (#SNPs = " << numAllAncSnps << ")\n";
        for (int vtxId = 0; vtxId < 3; vtxId++) {
            cout << "\tVertex " << vtxId << "\n";
            cout << "\t\tEUR: " << vtxPopExpGds[vtxId].e << "\n";
            cout << "\t\tAFR: " << vtxPopExpGds[vtxId].f << "\n";
            cout << "\t\tEAS: " << vtxPopExpGds[vtxId].a << "\n";
        }
    }
    
    return snpId;
}

int AncestrySnps::ReadRefSubPopSnpsFromFile(string refPopFile)
{
    ASSERT(FileExists(refPopFile.c_str()), "File " << refPopFile << " does not exist!\n");
    
    int snpId = 0;
    
    FILE* fp = fopen(refPopFile.c_str(), "r");
    if (fp == NULL) exit(EXIT_FAILURE);
    
    char* line = NULL;
    size_t len = 0;
    const string delim = "\t";
    int rsCol = 3;
    
    int numExpCols = 7;
    string expCols[7] = {"chr", "pos_37", "pos_38", "rs", "ref", "alt", "UKBBEUR"};
    
    // Make sure the file has correct columns
    bool isRightFile = true;
    if ((getline(&line, &len, fp)) != -1) {
        const string header(line);
        vector<string> vars = SplitString(header, delim);
        int numVars = vars.size();
        
        if (numVars <= numSubPops) {
            isRightFile = false;
        }
        else {
            for (int i = 0; i < numExpCols; i++) {
                string var = TrimString(vars[i]);

                if (var != expCols[i]) {
                    isRightFile = false;
                    break;
                }
            }
        }
    }
    
    string expColStr = expCols[0];
    for (int i = 1; i < numExpCols; i++) {
        expColStr += ", " + expCols[i];
    }
    expColStr += " ...";
    ASSERT(isRightFile, "File " << refPopFile << " does not have expected columns (" << expColStr << ")!\n");
    
    while ((getline(&line, &len, fp)) != -1) {
        ASSERT(snpId < numAllAncSnps, "File " << refPopFile << " has too many rows!\n");
        
        vector<string> vals = SplitString(string(line), delim);
        int rs = stoi(vals[rsCol]); 
        ASSERT(rs == snps[snpId].rs, "ERROR in " << refPopFile << ": Row #" << snpId << " has rs" << rs << " not " << snps[snpId].rs << "!\n");
        
        float popAfs[numSubPops];
        for (int i = 0; i < numSubPops; i++) {
            // TODO: Add code to make sure AF value is a number
            popAfs[i] = stof(vals[i+ancSnpFileOthCols]);
        }
        
        snps[snpId].SetRefSubPopAfs(popAfs);
        
        snpId++;
    }
    fclose(fp);
    
    if (line) free(line);
    
    return snpId;
}

int AncestrySnps::FindSnpIdGivenRs(int rsNum)
{
    int snpId = -1;

    if (rsToAncSnpId.find(rsNum) != rsToAncSnpId.end()) {
        snpId = rsToAncSnpId[rsNum];
    }

    return snpId;
}

int AncestrySnps::FindSnpIdGivenChrPos(int chr, int pos, int build, const char* ref, const char* alt)
{
    int snpId = -1;

    long int chrPos = long(chr) * 1000000000 + pos;
    
    if (sizeof(ref) == 8 && sizeof(alt) == 8) {
        if (build == 37) {
            if (pos37ToAncSnpId.find(chrPos) != pos37ToAncSnpId.end()) {
                snpId = pos37ToAncSnpId[chrPos];
            }
        }
        else if (build == 38) {
            if (pos38ToAncSnpId.find(chrPos) != pos38ToAncSnpId.end()) {
                snpId = pos38ToAncSnpId[chrPos];
            }
        }
        
        if (snpId > -1) {
            char eRef = snps[snpId].ref;
            char eAlt = snps[snpId].alt;
            char sRef = ref[0];
            char sAlt = alt[0];
            char cRef = FlipAllele(sRef);          
            char cAlt = FlipAllele(sAlt);
            
            if (!( (sRef == eRef && sAlt == eAlt) ||
                   (sRef == eAlt && sAlt == eRef) ||
                   (cRef == eRef && cAlt == eAlt) ||
                   (cRef == eAlt && cAlt == eRef)   ) ) {
                snpId = -1;
            }
        }
    }
    
    return snpId;
}

AncestrySnp AncestrySnps::GetAncestrySnp(int snpId)
{
    return snps[snpId];
}

void AncestrySnps::ShowAncestrySnps()
{
    int numAncSnps = snps.size();

    cout << "Total " << numAncSnps << " Ancestry SNPs.\n";
    bool debug = 0;

    if (debug) {
        for (int i = 0; i < 20; i++) {
            int snpId = i * 5000;
            AncestrySnp snp = snps[snpId];
            cout << "SNP " << snpId << " rs " << snp.rs << " chr " << snp.chr << " pos " << snp.posG37
            << " ref " << snp.ref << " alt " << snp.alt;
            for (int j = 0; j < numRefPops; j++) cout << " Ref " << j << " = " << snp.refPopAfs[j];
            for (int j = 0; j < numVtxPops; j++) cout << " Vtx " << j << " = " << snp.vtxPopAfs[j];
            cout << "\n";
        }

        cout << "Positions (x, y, z coordinates) of the three vertices when all SNPs have genotypes:\n";
        printf("\tE:  %5.4f  %5.4f  %5.4f\n", vtxPopExpGds[0].e, vtxPopExpGds[0].f, vtxPopExpGds[0].a);
        printf("\tF:  %5.4f  %5.4f  %5.4f\n", vtxPopExpGds[1].e, vtxPopExpGds[1].f, vtxPopExpGds[1].a);
        printf("\tA:  %5.4f  %5.4f  %5.4f\n", vtxPopExpGds[2].e, vtxPopExpGds[2].f, vtxPopExpGds[2].a);
    
        cout << "\nExpected genetic distances\n";
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cout << i << "-" << j << ": ";
                for (int k = 0; k < 5; k++) {
                    cout << vtxExpGenoDists[i][j][k] << " ";
                }
                cout << "\n";
            }
        }
    }
}

