#include "SampleGenoAncestry.h"

bool debug = false; /////////////////////////////////////////////// for debugging

GenoSample::GenoSample(string smp)
{
    name = smp;
    father = "";
    mother = "";
    sex = 0;

    numHetSnps = 0;
    hetRate = -1;

    numAncSnps = 0;
    ancIsSet = false;
}

void GenoSample::SetAncestryScores(int numSnps, float d1, float d2, float d3, float e, float f, float a, float rPe, float rPf, float rPa, bool isAnc)
{
    numAncSnps = numSnps;
    ancIsSet = isAnc;
    gd1 = d1;
    gd2 = d2;
    gd3 = d3;
    ePct = e;
    fPct = f;
    aPct = a;
    rawPe = rPe;
    rawPf = rPf;
    rawPa = rPa;
}

void GenoSample::SetSubPopGdScores(float* scores)
{
    for (int i = 0; i < numSubPopScores; i++)
        subPopScores[i] = scores[i];
}

SampleGenoAncestry::SampleGenoAncestry(AncestrySnps *aSnps, int minSnps)
{
    ancSnps = aSnps;
    if (minSnps) minAncSnps = minSnps;
    else         minAncSnps = 100;
    numSamples = 0;
    numAncSnps = 0;
    totAncSnps = ancSnps->GetNumAncestrySnps();

    ancSnpIds = NULL;
    ancSnpCodedGenos = NULL;

    samples = {};

    numThreads = 1;
    vtxExpGd0 = new SampleGenoDist(&aSnps->vtxPopExpGds[0], &aSnps->vtxPopExpGds[1],
        &aSnps->vtxPopExpGds[2], &aSnps->vtxPopExpGds[0]);

    vtxExpGd0->TransformAllDists();
    vtxExpGd0->CalculateBaryCenters();
}

SampleGenoAncestry::~SampleGenoAncestry()
{
    delete vtxExpGd0;
    samples.clear();
}

void SampleGenoAncestry::SetNumThreads(int threads)
{
    numThreads = threads;
}

void SampleGenoAncestry::SetGenoSamples(const vector<string> &smps)
{
    if (!smps.empty()) {
        numSamples = smps.size();

        for (int i = 0; i < numSamples; i++) {
            GenoSample genoSmp = GenoSample(smps[i]);
            samples.push_back(genoSmp);
        }
    }
    numAncSmps = 0;
}

void SampleGenoAncestry::SetGenoSamples(const vector<FamSample> &smps)
{
    if (!smps.empty()) {
        numSamples = smps.size();

        for (int i = 0; i < numSamples; i++) {
            GenoSample genoSmp = GenoSample(smps[i].name);
            samples.push_back(genoSmp);
        }
    }
    numAncSmps = 0;
}

void SampleGenoAncestry::SetSnpGenoData(vector<int> *snpIds, vector<unsigned char*> *snpCodedGenos)
{
    ancSnpIds = snpIds;
    ancSnpCodedGenos = snpCodedGenos;
    numAncSnps = ancSnpIds->size();
}

int SampleGenoAncestry::SaveAncestryResults(string outFile, bool isAppend)
{
    int numSaveSmps = 0;
    for (int i = 0; i < samples.size(); i++) {
        if (samples[i].ancIsSet) numSaveSmps++;
    }

    if (numSaveSmps < 1) {
        cout << "\nNOTE: None of the " << numSamples << " samples have enough genotypes for ancestry inference."
        <<  " No ancestry results were generated.\n";
        return numSaveSmps;
    }

    string vtxTitle = "Positions (x, y, z coordinates) of the three vertices";
    // vtxExpGd0->ShowPositions(vtxTitle);

    FILE *ifp = NULL;
    if (isAppend) ifp = fopen(outFile.c_str(), "a");
    else          ifp = fopen(outFile.c_str(), "w");
    
    if(ifp) {
        char line[512];
        if (!isAppend) {
            sprintf(line, "%s\t%s\tGD1\tGD2\tGD3", "Sample", "SNPs");
            for (int i = 0; i < numSubPopScores; i++) {
                strcat(line, "\t");
                strcat(line, popScoreNames[i].c_str());
            }
            strcat(line, "\tPe\tPf\tPa\tRawPe\tRawPf\tRawPa\tAncGroupID");

            fprintf(ifp, "%s\n", line);
        }

        for (int i = 0; i < numSaveSmps; i++) {
            GenoSample smp = samples[i];
            if (!smp.ancIsSet) continue;

            sprintf(line, "%s\t%d", smp.name.c_str(), smp.numAncSnps);
            
            sprintf(line, "%s\t%7.6f\t%7.6f\t%7.6f", line, smp.gd1, smp.gd2, smp.gd3);
            for (int i = 0; i < numSubPopScores; i++) {
                sprintf(line, "%s\t%7.6f", line, smp.subPopScores[i]);
            }
            sprintf(line, "%s\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%d",
                    line, smp.ePct, smp.fPct, smp.aPct, smp.rawPe, smp.rawPf, smp.rawPa, smp.ancGroupId);

            fprintf(ifp, "%s\n", line);
        }
    }
    else {
        cout << "ERROR: Can't open " << outFile << " for writing!\n";
        return 0;
    }

    fclose(ifp);
    if (isAppend) cout << "Added population results of " << numSaveSmps << " samples to " << outFile << ".\n";
    else          cout << "Saved population results of " << numSaveSmps << " samples to " << outFile << ".\n";
        
    
    return numSaveSmps;
}

void SampleGenoAncestry::CalculateSubPopGd0Values()
{
    double gdScoreSumP1[numSubPopScores];
    double gdScoreSumP2[numSubPopScores];
    int numScoreSnps[numSubPopScores];
    
    for (int i = 0; i < numSubPopScores; i++) {
        gdScoreSumP1[i] = 0;
        gdScoreSumP2[i] = 0;
        numScoreSnps[i] = 0;
    }
    
    for (int snpNo = 0; snpNo < numAllAncSnps; snpNo++) {
        if (debug && snpNo % 50000 == 0) cout << "SNP No " << snpNo << "\n";

        for (int i = 0; i < numSubPopScores; i++) {
            int refPopId1 = scorePopIdx1[i];
            int refPopId2 = scorePopIdx2[i];
            
            float p1 = ancSnps->snps[snpNo].refSubPopAfs[refPopId1];
            float p2 = ancSnps->snps[snpNo].refSubPopAfs[refPopId2];
            
            // Only count those SNPs for which there are both ref and normalization pop freqs
            if (p1 > 0 && p1 < 1 && p2 > 0 && p2 < 1) {
            double logp = log(p2/p1);
                double logq = log((1-p2)/(1-p1));
                
                // Score from only one of the two alleles is added              
                gdScoreSumP1[i] += (p1 * logp + (1-p1) * logq) * 2;
                gdScoreSumP2[i] += (p2 * logp + (1-p2) * logq) * 2;
                
                numScoreSnps[i]++;
            }
            
            if (debug && snpNo % 50000 == 0)
                printf("\tPop %d: p1 %5.4f p2 %5.4f; n = %d S1 %7.4f S2 %7.4f\n",
                       i, p1, p2, numScoreSnps[i], gdScoreSumP1[i], gdScoreSumP2[i]);
        }
    }
    
    for (int i = 0; i < numSubPopScores; i++) {
        // Subject GD scores are calculated for two alleles on each SNP
        subPopGd0P1[i] = gdScoreSumP1[i] / numScoreSnps[i];
        subPopGd0P2[i] = gdScoreSumP2[i] / numScoreSnps[i];
    }
}

void SampleGenoAncestry::SetAncestryPvalues(int thNo)
{
    int meanThSmps = int(numSamples / numThreads);
    int rmSmps = numSamples % numThreads;
    int chkThSmps = meanThSmps;
    if (thNo+1 <= rmSmps) chkThSmps = meanThSmps + 1;

    int stSmp = thNo * chkThSmps;
    if (thNo >= rmSmps) stSmp = thNo * chkThSmps + rmSmps;
    int edSmp = stSmp + chkThSmps - 1;

    int smpCnt = 0;
    for (int smpNo = stSmp; smpNo <= edSmp; smpNo++) {
        // Calculate 14 scores for each sample, based on the SNPs with genotypes, i.e.,
        // 9 expected genetics distances from the 3 vertices to the first 3 reference populations, and
        // 5 genetic distances from the sample to the 5 referene populations

        double popPvalues[numRefPops];
        double popMeanPvals[numRefPops];

        for (int popId = 0; popId < numRefPops; popId++) {
            popPvalues[popId] = 0;   // The raw log p-value
            popMeanPvals[popId] = 0; // Genetic distancs from the sample to each ref population
        }

        // Assuming the three vertex populations have allele freqs for all Ancestry SNPs
        int numGenoSnps = 0;
        double vtxExpPeSums[3];
        double vtxExpPfSums[3];
        double vtxExpPaSums[3];

        // Initialize expected values for each vertex population (E, F, A)
        for (int vtxId = 0; vtxId < 3; vtxId++) {
            vtxExpPfSums[vtxId] = 0;
            vtxExpPeSums[vtxId] = 0;
            vtxExpPaSums[vtxId] = 0;
        }

        // Declare and initialize values for calculating subpopulation delta GD scores
        double smpGdScoreSum[numSubPopScores];    // Sum of sample GD over SNPs
        double nomGdScoreSumP1[numSubPopScores];  // Sum of norm pop1 GD over SNPs
        double nomGdScoreSumP2[numSubPopScores];  // Sum of norm pop2 GD over SNPs
        double smpGdScoreSnpNum[numSubPopScores];
        
        for (int i = 0; i < numSubPopScores; i++) {
            smpGdScoreSum[i] = 0;
            nomGdScoreSumP1[i] = 0;
            nomGdScoreSumP2[i] = 0;
            smpGdScoreSnpNum[i] = 0;
        }
        
        float smpGdRawScore[numSubPopScores];
        float smpGdScore[numSubPopScores];
        float nomGdScoreP1[numSubPopScores];
        float nomGdScoreP2[numSubPopScores];

        // In case some reference populations might not have allele freqs for all Ancestry SNPs
        int refPopSnps[numRefPops]; // Counts of SNPs with freqs for each ref population
        for (int popId = 0; popId < numRefPops; popId++) {
            refPopSnps[popId] = 0;
        }

        bool hasRefPv = false;
        int ancSnpNo = 0;
        int genoChrNo = smpNo / 4; // Each char stores genotypes of 4 samples
        int chrGenoIdx = smpNo % 4;

        int bit1Pos = chrGenoIdx * 2;
        int bit2Pos = bit1Pos + 1;
        
        // Count heterozygous genotypes for each sample
        int hetSnps = 0;
        
        for (ancSnpNo = 0; ancSnpNo < numAncSnps; ancSnpNo++) {
            unsigned char* chrGenos = (*ancSnpCodedGenos)[ancSnpNo];
            unsigned char chrGeno = chrGenos[genoChrNo];
            int gsize = sizeof(chrGenos);

            unsigned int geno = 0;
            if (chrGeno & charGenoBitVals[bit2Pos]) geno++;
            if (chrGeno & charGenoBitVals[bit1Pos]) geno += 2;

            int snpNo = (*ancSnpIds)[ancSnpNo];

            // Alt allele freq p of the 3 vertices
            double v0p = ancSnps->snps[snpNo].vtxPopAfs[0];
            double v1p = ancSnps->snps[snpNo].vtxPopAfs[1];
            double v2p = ancSnps->snps[snpNo].vtxPopAfs[2];

            if (geno < 3) {
                if (0) {
                    cout << "    Anc snp No. " << ancSnpNo << " Snp No. " << snpNo
                    << " rs " << ancSnps->snps[snpNo].rs
                    << " ref " << ancSnps->snps[snpNo].ref
                    << " alt " << ancSnps->snps[snpNo].alt
                    << " v0 " << ancSnps->snps[snpNo].vtxPopAfs[0]
                    << " p0 " << ancSnps->snps[snpNo].refPopAfs[0]
                    << " geno " << geno << "\n";
                }

                for (int popId = 0; popId < numRefPops; popId++) {
                    float pv = ancSnps->snps[snpNo].refPopAfs[popId]; // reference population allele freq p value

                    if (pv > 0 && pv < 1) {
                        float qv = 1 - pv;

                        if      (geno == 0) {
                            popPvalues[popId] += log(qv) * 2;
                        }
                        else if (geno == 1) {
                            popPvalues[popId] += log(pv * qv * 2);
                        }
                        else if (geno == 2) {
                            popPvalues[popId] += log(pv) * 2;
                        }

                        hasRefPv = true;
                        refPopSnps[popId]++;
                    }
                }

                // Calculate expected P values for each of the vertices E, F, A
                for (int vtxId = 0; vtxId < 3; vtxId++) {
                    vtxExpPeSums[vtxId] += ancSnps->vtxExpGenoDists[vtxId][0][snpNo];
                    vtxExpPfSums[vtxId] += ancSnps->vtxExpGenoDists[vtxId][1][snpNo];
                    vtxExpPaSums[vtxId] += ancSnps->vtxExpGenoDists[vtxId][2][snpNo];
                }

                // Calculate extra GD scores from subpopulation references
                for (int sId = 0; sId < numSubPopScores; sId++) {
                    int refPopId1 = scorePopIdx1[sId];
                    int refPopId2 = scorePopIdx2[sId];

                    float p1 = ancSnps->snps[snpNo].refSubPopAfs[refPopId1];
                    float p2 = ancSnps->snps[snpNo].refSubPopAfs[refPopId2];
                    
                    // Calculate three values for each SNP locus: sample GD, and the two norm pop GDs
                    if (p1 > 0 && p1 < 1 && p2 > 0 && p2 < 1) {
                        float logp = log(p2/p1);
                        float logq = log((1-p2)/(1-p1));
                        
                        if      (geno == 2) smpGdScoreSum[sId] += logp * 2;
                        else if (geno == 1) smpGdScoreSum[sId] += logp + logq;
                        else if (geno == 0) smpGdScoreSum[sId] += logq * 2;

                        nomGdScoreSumP1[sId] += (p1 * logp + (1-p1) * logq) * 2;
                        nomGdScoreSumP2[sId] += (p2 * logp + (1-p2) * logq) * 2;
                        
                        smpGdScoreSnpNum[sId]++;
                    }                      
                }

                if (geno == 1) hetSnps++;
                numGenoSnps++;
            }
        }

        for (int i = 0; i < numSubPopScores; i++) {
            smpGdRawScore[i] = smpGdScoreSum[i]   / smpGdScoreSnpNum[i];
            nomGdScoreP1[i]  = nomGdScoreSumP1[i] / smpGdScoreSnpNum[i];
            nomGdScoreP2[i]  = nomGdScoreSumP2[i] / smpGdScoreSnpNum[i];
            // Normalize the GD score
            smpGdScore[i] = subPopGdNormP1 + (subPopGdNormP2-subPopGdNormP1)*(smpGdRawScore[i]-nomGdScoreP1[i])/(nomGdScoreP2[i]-nomGdScoreP1[i]);
        }

        for (int popId = 0; popId < numRefPops; popId++) {
            if (refPopSnps[popId] > 0) {
                popMeanPvals[popId] = -1 * popPvalues[popId]/refPopSnps[popId];
            }
        }

        float gd1 = 0, gd2 = 0, gd3 = 0;
        float ePct = 0, fPct = 0, aPct = 0;
        float rPe = 0, rPf = 0, rPa = 0;    // Raw E, F, A percentages
        bool hasAncGeno = false;

        if (numGenoSnps >= minAncSnps) {
            GenoDist smpDist;
            GenoDist vtxExpDists[numVtxPops];

            smpDist.e = popMeanPvals[0];
            smpDist.f = popMeanPvals[1];
            smpDist.a = popMeanPvals[2];

            for (int vtxId = 0; vtxId < numVtxPops; vtxId++) {
                vtxExpDists[vtxId].e  = -1 * vtxExpPeSums[vtxId]/numGenoSnps;
                vtxExpDists[vtxId].f  = -1 * vtxExpPfSums[vtxId]/numGenoSnps;
                vtxExpDists[vtxId].a  = -1 * vtxExpPaSums[vtxId]/numGenoSnps;
            }

            // Calculate GD and ancestry components using the raw scores
            SampleGenoDist *smpGd = new SampleGenoDist(&vtxExpDists[0], &vtxExpDists[1], &vtxExpDists[2], &smpDist);
            smpGd->TransformAllDists();
            smpGd->CalculateBaryCenters();
            
            // Show rotated x, y, z values as GD1, GD2, GD3
            gd1 = smpGd->eWt * vtxExpGd0->ePt.x + smpGd->fWt * vtxExpGd0->fPt.x + smpGd->aWt * vtxExpGd0->aPt.x;
            gd2 = smpGd->eWt * vtxExpGd0->ePt.y + smpGd->fWt * vtxExpGd0->fPt.y + smpGd->aWt * vtxExpGd0->aPt.y;
            gd3 = smpGd->sPt.z;

            float ejWt = smpGd->eWt > 0 ? smpGd->eWt : 0;
            float fjWt = smpGd->fWt > 0 ? smpGd->fWt : 0;
            float ajWt = smpGd->aWt > 0 ? smpGd->aWt : 0;
            float totWt = fjWt + ejWt + ajWt;
            ePct = ejWt * 100 / totWt;
            fPct = fjWt * 100 / totWt;
            aPct = ajWt * 100 / totWt;
            rPe = smpGd->eWt * 100;
            rPf = smpGd->fWt * 100;
            rPa = smpGd->aWt * 100;

            hasAncGeno = true;
            numAncSmps++;
            delete smpGd;
        }

        samples[smpNo].SetAncestryScores(numGenoSnps, gd1, gd2, gd3, ePct, fPct, aPct, rPe, rPf, rPa, hasAncGeno);
        samples[smpNo].SetSubPopGdScores(smpGdScore);
        samples[smpNo].SetGrafAncGroups();
        
        samples[smpNo].numHetSnps = hetSnps;
        samples[smpNo].hetRate = numGenoSnps > 0 ? hetSnps * 1.0 / numGenoSnps : -1;
        
        smpCnt++;
        if (thNo == 0 && smpCnt % 100 == 0)
            cout  << "\tCalculated scores for " << smpCnt << " of total " << chkThSmps << " samples\n";
    }
}

void SampleGenoAncestry::ShowSummary()
{
    cout << "Num samples: " << numSamples << "\n";

    cout << "Tot ancestry SNPs: " << totAncSnps << "\n";
    cout << "Num ancestry SNPs in dataset: " << numAncSnps << "\n";

    for (int i = 0; i < numAncSnps; i++) {
        int snpNo = (*ancSnpIds)[i];
        cout << "No. " << i << ": " << snpNo << ": ";

        for (int j = 0; j < numSamples; j++) {
            if  (j < 20) cout << int((*ancSnpCodedGenos)[i][j]) << " ";
        }
        cout << "\n";

        if (i > 20) break;
    }
}

void GenoSample::SetGrafAncGroups()
{
    float ea1 = subPopScores[0];
    float ea2 = subPopScores[1];
    // float ea3 = subPopScores[2];
    float ea4 = subPopScores[3];
    float af1 = subPopScores[4];
    float af2 = subPopScores[5];
    float af3 = subPopScores[6];
    float eu1 = subPopScores[7];
    float eu2 = subPopScores[8];
    float eu3 = subPopScores[9];
    float sa1 = subPopScores[10];
    float sa2 = subPopScores[11];
    float ic1 = subPopScores[12];
    float ic2 = subPopScores[13];
    float ic3 = subPopScores[14];
    
    int ancId = 0;
    
    if (aPct > 50 && fPct > 10 && ic1 > 0.4 && gd3 > 0.035) {
        ancId = 700; // Oceania
    }                
    else if (ic1 > 0.5 && fPct < 15) {
        // South Asia
        
        if (sa2 < -0.1) {
            ancId = 402;  // Gujarati
        }
        else if (sa2 > 0.7) {
            ancId = 403;  // Pakistan
        }
        else if (sa1 < -0.6) {
            ancId = 404;  // Sri Lanka
        }
        else if (sa1 > 0.1) {
            ancId = 405;  // Bangaladesh
        }
        else {
            ancId = 401;  // India
        }
    }
    else if (aPct > 15 && fPct > 15) {
        ancId = 800;      // Multiple ancestry
    }
    else if (ic1 > -0.3 && aPct > 40) {
        // East Asia
        
        if (ea2 > 1.4) {
            ancId = 501; // Ryukyu
        }
        else if (ea2 > 0.4) {
            ancId = 502; // Japan
        }
        else if (ea4 > -0.7 && ea4 < -0.2 && ea1 < -0.3) {
            if (ea4 < -0.48) {
                ancId = 506; // Northern China 2 (Daur, Xibo, Hezhe) 
            }
            else {
                ancId = 508; // Southern China 2 (Naxi, Tu, Yi)
            }
        }
        else if (ea2 > -0.3 && ea1 < -0.3) {
            ancId = 503; // Korea
        }
        else if (ea4 > -0.2 && ea4 > ea1 + 0.1) {
            ancId = 504; // Northern Asian
        }
        else if (ea1 < -0.6) {
            ancId = 505; // Northern China 1
        }
        else if (ea1 < -0.1) {
            ancId = 507; // Southern China 1
        }
        else if (ea1 < 0.5) {
            ancId = 509; // Southeast Asian
        }
        else if (ea2 > -1) {
            ancId = 510; // Thailand
        }
        else {
            ancId = 511; // Other East Asia
        }
    }
    else if (eu1 < 1.6 && ic2 < -0.25) {
        // Europe and part of MENA
      
        if (ic3 < 3.5 - 3.1 * eu1) {
            // Europe
            
            if (ic3 > -0.3) {
                ancId = 308; // Other Europe
            }
            else if (eu2 > 0.3) {
                ancId = 301; // Finland
            }
            else if (eu3 > -0.4 - eu1) {
                if (eu3 > eu1 + 1.1) {
                    ancId = 305; // Northeast Europe
                }                            
                else if (eu3 < 1.3 * eu1 - 1.3) {
                    ancId = 307; // Balkans
                }
                else {
                    ancId = 306; // Southeast Europe
                }                            
            }
            else {
                if (eu3 < eu1 - 1.0) {
                    ancId = 304; // South Europe
                }                            
                else if (eu2 > -0.75 &&  eu1 < -0.75) {
                    ancId = 302; // North Europe
                }                            
                else {
                    ancId = 303; // West Europe
                }                            
            }
        }
        else {
            // MENA
          
            if (ic3 > 0.3) {
                ancId = 203; // Middle East 1
            }
            else {
                ancId = 202; // Middle East 2
            }
        }
    }
    else if (eu1 > 1.6 && ic2 < 0.4) {
        if (ic3 < -0.3) {
            ancId = 201; // North Africa
        }
        else {
            ancId = 202; // Middle East 2
        }
    }
    else if (eu1 > 2.2 && ic2 < 1.4) {
          ancId = 106; // Northeast Africa
    }
    else if (gd1  > 1.4758) {
        // Latin American 1
        
        if (ic1 < -0.2) {
            if (gd3 > 0.05) {
                ancId = 603; // Native American
            }
            else {
                ancId = 602; // Other Latin American
            }
        }
        else {
            ancId = 308; // Other Europe
        }
    }
    else {
        // Black backgound
        
        if (ic1 < 0.9 * gd1 - 1.36) {
            ancId = 601; // Latin American 1
        }
        else {
            // West Africa
            
            if (af3 < af1 * 4.0 / 3 + 0.4) {
                if (af3 > 0.5) {
                    ancId = 108; // Other Africa
                }
                else {
                    if (af1 > 0.6) {
                        if (af3 > -0.1) {
                            ancId = 104; // Kenya
                        }
                        else {
                            ancId = 105; // Southern Africa
                        }
                    }
                    else {
                        ancId = 103;     // Central Africa
                    }
                }
            }
            else {
                if (af1 > -0.8 && af2 > -0.6 && af2 < 0.9 * af1 + 0.8) {
                    ancId = 107; // African American
                }
                else {
                    if (af3 > 0.5) {
                        ancId = 108; // Other Africa
                    }
                    else if (af2 < -0.1) {
                        ancId = 101; // Nigeria
                    }
                    else {
                        ancId = 102; // West Africa
                    }
                }
            }
        }
    } 
    
    ancGroupId = ancId;
}
