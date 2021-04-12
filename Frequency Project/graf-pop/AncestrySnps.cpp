#include "AncestrySnps.h"

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

AncestrySnps::AncestrySnps()
{

    refPopNames[0] = "African";
    refPopNames[1] = "European";
    refPopNames[2] = "Asian";
    refPopNames[3] = "Mexican";
    refPopNames[4] = "Indian-Pakistani";

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

    int lineLen = 300;
    char snpLine[lineLen];

    FILE *ifp = fopen(ancSnpFile.c_str(), "r");
    ASSERT(ifp, "Couldn't open file " << ancSnpFile << ".\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int numSnps = 0;
    int rsNum, chr, g37, g38;
    char a1, a2;
    float rfEur, rfAfa, rfAsn, rfLat, rfSas, vtEur, vtAfr, vtEas;

    while (fgets(snpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        if (lineNo == 0) {
            if (snpLine[0] != 'c' || snpLine[1] != 'h' || snpLine[2] != 'r') {
                fileIsValid = false;
            }
        }
        else {
            sscanf(snpLine, "%d %d %d %d %c %c %f %f %f %f %f %f %f %f",
            &chr, &g37, &g38, &rsNum, &a1, &a2, &rfEur, &rfAfa, &rfAsn, &rfLat, &rfSas, &vtEur, &vtAfr, &vtEas);

            float* refPopAfs = new float[numRefPops];
            float* vtxPopAfs = new float[numVtxPops];

            refPopAfs[0] = rfEur;
            refPopAfs[1] = rfAfa;
            refPopAfs[2] = rfAsn;
            refPopAfs[3] = rfLat;
            refPopAfs[4] = rfSas;

            vtxPopAfs[0] = vtEur;
            vtxPopAfs[1] = vtAfr;
            vtxPopAfs[2] = vtEas;

            AncestrySnp ancSnp(numSnps, rsNum, chr, g37, g38, a1, a2, refPopAfs, vtxPopAfs);

            snps.push_back(ancSnp);

            long int chrPos37 = (long)chr * 1000000000 + g37;
            long int chrPos38 = (long)chr * 1000000000 + g38;

            rsToAncSnpId[rsNum] = numSnps;
            pos37ToAncSnpId[chrPos37] = numSnps;
            pos38ToAncSnpId[chrPos38] = numSnps;

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

                vtxExpGenoDists[vtxId][0][numSnps] = eGd;
                vtxExpGenoDists[vtxId][1][numSnps] = fGd;
                vtxExpGenoDists[vtxId][2][numSnps] = aGd;

                popExpPeSums[vtxId] += eGd;
                popExpPfSums[vtxId] += fGd;
                popExpPaSums[vtxId] += aGd;
            }

            delete refPopAfs;
            delete vtxPopAfs;

            numSnps++;
        }

        lineNo++;
    }
    fclose(ifp);

    ASSERT(numSnps == numAncSnps, "numSnps = " << numAncSnps << ".\n");

    for (int vtxId = 0; vtxId < 3; vtxId++) {
        vtxPopExpGds[vtxId].e = -1 * popExpPeSums[vtxId]/numSnps;
        vtxPopExpGds[vtxId].f = -1 * popExpPfSums[vtxId]/numSnps;
        vtxPopExpGds[vtxId].a = -1 * popExpPaSums[vtxId]/numSnps;
    }

    cout << "Read " << numSnps << " ancestry SNPs from file " << ancSnpFile << "\n\n";
    if (0) {
        cout << "Expected vertex genetic distances\n";
        for (int vtxId = 0; vtxId < 3; vtxId++) {
            cout << "\tVertex " << vtxId << "\n";
            cout << "\t\tEUR: " << vtxPopExpGds[vtxId].e << "\n";
            cout << "\t\tAFR: " << vtxPopExpGds[vtxId].f << "\n";
            cout << "\t\tEAS: " << vtxPopExpGds[vtxId].a << "\n";
        }
    }
}

int AncestrySnps::FindSnpIdGivenRs(int rsNum)
{
    int snpId = -1;

    if (rsToAncSnpId.find(rsNum) != rsToAncSnpId.end()) {
        snpId = rsToAncSnpId[rsNum];
    }

    return snpId;
}

int AncestrySnps::FindSnpIdGivenChrPos(int chr, int pos, int build)
{
    int snpId = -1;

    long int chrPos = long(chr) * 1000000000 + pos;

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
