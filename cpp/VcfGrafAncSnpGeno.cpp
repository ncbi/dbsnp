#include "VcfGrafAncSnpGeno.h"

VcfGrafAncSnpGeno::VcfGrafAncSnpGeno(string file, AncestrySnps *aSnps)
{
    ancSnps = aSnps;
    vcfFile = file;

    vcfSamples = {};
    vcfAncSnpGtVals = {};
    vcfAncSnpChrs = {};
    vcfAncSnpPoss = {};
    vcfAncSnpSnps = {};
    vcfAncSnpRefs = {};
    vcfAncSnpAlts = {};
    vcfRsIdAncSnpIds = {};
    vcfGb37AncSnpIds = {};
    vcfGb38AncSnpIds = {};
    vcfAncSnpCodedGenos = {};
    vcfAncSnpIds = {};

    totAncSnps = ancSnps->GetNumAncestrySnps();

    numSamples = 0;
    totVcfSnps = 0;
    putativeAncSnps = 0;
    numRsIdAncSnps = 0;
    numGb37AncSnps = 0;
    numGb38AncSnps = 0;
    numVcfAncSnps = 0;
}

VcfGrafAncSnpGeno::~VcfGrafAncSnpGeno()
{
    vcfSamples.clear();
    vcfAncSnpChrs.clear();
    vcfAncSnpPoss.clear();
    vcfAncSnpSnps.clear();
    vcfAncSnpRefs.clear();
    vcfAncSnpAlts.clear();
    vcfRsIdAncSnpIds.clear();
    vcfGb37AncSnpIds.clear();
    vcfGb38AncSnpIds.clear();
    vcfAncSnpIds.clear();
    
    DeleteAncSnpCodedGenos();
    DeleteAncSnpGtValues();
}

void VcfGrafAncSnpGeno::DeleteAncSnpCodedGenos()
{
    for (int i = 0; i < vcfAncSnpCodedGenos.size(); i++) {
        if (vcfAncSnpCodedGenos[i]) delete vcfAncSnpCodedGenos[i];
    }
    vcfAncSnpCodedGenos.clear();
}

void VcfGrafAncSnpGeno::DeleteAncSnpGtValues()
{
    for (int i = 0; i < vcfAncSnpGtVals.size(); i++) {
        if (!vcfAncSnpGtVals[i].empty()) vcfAncSnpGtVals[i].clear();
    }
    vcfAncSnpGtVals.clear();
}

bool VcfGrafAncSnpGeno::ReadHeaderFromFile()
{
    htsFile *fp = hts_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening VCF file " << vcfFile << std::endl;
        return false;
    }

    bcf_hdr_t *header = bcf_hdr_read(fp);
    if (!header) {
        std::cerr << "Error reading VCF header from file " << vcfFile << std::endl;
        return false;
    }

    numVcfSamples = bcf_hdr_nsamples(header);

    bcf_hdr_destroy(header);
    hts_close(fp);
    
    return true;
}

bool VcfGrafAncSnpGeno::ReadDataFromFile()
{
    return ReadDataFromFile(-1, 0, true);
}

bool VcfGrafAncSnpGeno::ReadDataFromFile(int stSbjNo, int numReadSbjs, bool verbose)
{
    unsigned int charBaseInts[4] = {64, 16, 4, 1};
    
    htsFile *fp = hts_open(vcfFile.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening VCF file " << vcfFile << std::endl;
        return false;
    }
    
    bcf_hdr_t *header = bcf_hdr_read(fp);
    if (!header) {
        std::cerr << "Error reading VCF header from file " << vcfFile << std::endl;
        return false;
    }

    int totVcfSamples = bcf_hdr_nsamples(header);

    if (stSbjNo >= 0 || numReadSbjs < totVcfSamples) {
        assert(stSbjNo % 4 == 0);
        if (stSbjNo + numReadSbjs > totVcfSamples) numReadSbjs = totVcfSamples - stSbjNo;
        
        numSamples = numReadSbjs;
        for (int i = 0; i < numSamples; i++) vcfSamples.push_back(header->samples[i+stSbjNo]);
    }
    else {
        numSamples = totVcfSamples;
        for (int i = 0; i < numSamples; i++) vcfSamples.push_back(header->samples[i]);
    }
    
    numGenoChars = (numSamples - 1) / 4 + 1; // Save 4 genotypes to one char
    if (numSamples < totVcfSamples) cout << numSamples << " samples were read to memory";

    bcf1_t *record = bcf_init();

    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    
    int snpNo = 0;
    int valSnpNo = 0;
    int numVcfSnps = 0;
    
    while (bcf_read(fp, header, record) == 0) {
        bcf_unpack(record, BCF_UN_STR);
  
        string chrStr = string(bcf_hdr_id2name(header, record->rid));
        int32_t pos = record->pos + 1;
        string snpStr = string(record->d.id);
        const char* gRef = record->d.allele[0];
        const char* gAlt = record->d.allele[1];
        string refStr = string(gRef);
        string altStr = string(gAlt);
        
        ngt = bcf_get_genotypes(header, record, &gt, &ngt_arr);

        int chr = GetChromosomeFromString(chrStr.c_str());
        int rsNum = GetRsNumFromString(snpStr.c_str());
        
        int rsSnpId = ancSnps->FindSnpIdGivenRs(rsNum);
        int gb37SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 37, gRef, gAlt);
        int gb38SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 38, gRef, gAlt);

        // Testing
        if (snpNo < 0) {        
            cout << "chr " << chr << " pos " << pos << " rs " << rsNum << " ref " << gRef << " alt " << gAlt << " GT: " << ngt;
            cout << " rsSnpId " << rsSnpId << " 37Id " << gb37SnpId << " 38Id " << gb38SnpId << "\n";
        }
        
        if (ngt == totVcfSamples * 2) {
            // Testing
            if (snpNo < 0) {        
                for (int i = 0; i < 20; i++) {
                    int g1 = bcf_gt_allele(gt[i*2]);
                    int g2 = bcf_gt_allele(gt[i*2+1]);
                    cout << " " << g1 << "|" << g2;
                }
                cout << "\n";
            }
            
            if (rsSnpId > -1 || gb37SnpId > -1 || gb38SnpId > -1) {
                putativeAncSnps++;
                if (rsSnpId > -1)   numRsIdAncSnps++;
                if (gb37SnpId > -1) numGb37AncSnps++;
                if (gb38SnpId > -1) numGb38AncSnps++;
                
                vcfRsIdAncSnpIds.push_back(rsSnpId);
                vcfGb37AncSnpIds.push_back(gb37SnpId);
                vcfGb38AncSnpIds.push_back(gb38SnpId);
                vcfAncSnpSnps.push_back(snpStr);
                vcfAncSnpRefs.push_back(refStr);
                vcfAncSnpAlts.push_back(altStr);
                
                vector<unsigned char> gtVals;
                for (int charNo = 0; charNo < numGenoChars; charNo++) {
                    unsigned char charGeno = 0;

                    for (int charGenoNo = 0; charGenoNo < 4; charGenoNo++) {
                        int smpNo = charNo * 4 + charGenoNo;
                        if (stSbjNo >= 0) smpNo += stSbjNo;
                        
                        if (smpNo < totVcfSamples) {
                            int refGeno = bcf_gt_allele(gt[smpNo*2]);
                            int altGeno = bcf_gt_allele(gt[smpNo*2+1]);
                            unsigned char numAlts = 3;
                            
                            if (refGeno == 0 && altGeno == 0) {
                                numAlts = 0;
                            }
                            else if (refGeno == 0 && altGeno == 1) {
                                numAlts = 1;
                            }
                            else if (refGeno == 1 && altGeno == 0) {
                                numAlts = 1;
                            }
                            else if (refGeno == 1 && altGeno == 1) {
                                numAlts = 2;
                            }
                            
                            charGeno += numAlts * charBaseInts[charGenoNo];
                        }
                    }
                    gtVals.push_back(charGeno);
                }

                vcfAncSnpGtVals.push_back(gtVals);
            }
            
            numVcfSnps++;

            valSnpNo++;
        }
        else {
            cout << "\nERROR at line #" << snpNo << ": chr " << chr
                 << ", pos " << pos << ", snp " << snpStr
                 << ". #genotypes (" << ngt << ") is different from #samples (" << numSamples << ").\n";
            return false;
        }
        
        snpNo++;
        
        if (verbose && snpNo % 100000 == 0) {
            cout << "\tChecked " << snpNo << " SNPs. Found " << putativeAncSnps << " SNPs with ancestry SNPs.\n";
        }
    }

    cout << "Checked " << snpNo << " SNPs. Found " << putativeAncSnps << " SNPs with ancestry SNPs\n";
    totVcfSnps = snpNo;
    
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    hts_close(fp);

    return true;
}

void VcfGrafAncSnpGeno::RecodeSnpGenotypes()
{
    // Rs ID, GB37, or GB38, use whichever returns the most ancestry SNPs to find these SNPs
    ancSnpType = AncestrySnpType::RSID;
    int maxVcfAncSnps = numRsIdAncSnps;
    
    if (numGb37AncSnps > maxVcfAncSnps) {
        ancSnpType = AncestrySnpType::GB37;
        maxVcfAncSnps = numGb37AncSnps;
    }
    
    if (numGb38AncSnps > maxVcfAncSnps) {
        ancSnpType = AncestrySnpType::GB38;
        maxVcfAncSnps = numGb38AncSnps;
    }

    int saveSnpNo = 0; // Putative SNPs saved after vcf file was read
    int ancSnpNo = 0;  // Final list of ancestry SNPs to be used for ancestry inference
    
    unsigned int charBaseInts[4] = {64, 16, 4, 1};
    unsigned int baseNums[8];
    unsigned int base = 1;
    for (int i = 0; i < 8; i++) {
        baseNums[i] = base;
        base = base << 1;
    }

    map<int, int> allAncSnpIds;
    for (saveSnpNo = 0; saveSnpNo < putativeAncSnps; saveSnpNo++) {
        int ancSnpId = -1;
        if      (ancSnpType == AncestrySnpType::RSID) ancSnpId = vcfRsIdAncSnpIds[saveSnpNo];
        else if (ancSnpType == AncestrySnpType::GB37) ancSnpId = vcfGb37AncSnpIds[saveSnpNo];
        else if (ancSnpType == AncestrySnpType::GB38) ancSnpId = vcfGb38AncSnpIds[saveSnpNo];

        if (ancSnpId > -1 && !allAncSnpIds[ancSnpId]) {
            char eRef = ancSnps->snps[ancSnpId].ref;
            char eAlt = ancSnps->snps[ancSnpId].alt;

            string vcfRef = vcfAncSnpRefs[saveSnpNo];
            string vcfAlt = vcfAncSnpAlts[saveSnpNo];
            
            int expRefIdx = -1;
            int expAltIdx = -1;
            
            CompareAncestrySnpAlleles(vcfRef, vcfAlt, eRef, eAlt, &expRefIdx, &expAltIdx);

            if (expRefIdx > -1 && expAltIdx > -1) {
                unsigned char* smpGenos = new unsigned char[numGenoChars];

                for (int charNo = 0; charNo < numGenoChars; charNo++) {
                    int charGtVal = int(vcfAncSnpGtVals[saveSnpNo][charNo]);
                    unsigned char charGeno = 0;

                    for (int cGenoNo = 0; cGenoNo < 4; cGenoNo++) {
                        int smpNo = charNo * 4 + cGenoNo;
                        if (smpNo < numSamples) {
                            int numVcfAlts = charGtVal / charBaseInts[cGenoNo];
                            charGtVal -= numVcfAlts * charBaseInts[cGenoNo];
                            
                            unsigned char geno = 3;                            
                            
                            if (expAltIdx == 1) {
                                geno = numVcfAlts; 
                            }
                            else if (expAltIdx == 0) {
                                geno = numVcfAlts; 
                                if (numVcfAlts < 3) {
                                    geno = 2 - numVcfAlts;
                                }
                            }
                            
                            charGeno += geno * charBaseInts[cGenoNo];
                        }
                    }
                    
                    smpGenos[charNo] = charGeno;
                }
                
                vcfAncSnpIds.push_back(ancSnpId);
                vcfAncSnpCodedGenos.push_back(smpGenos);
                vcfAncSnpGtVals[saveSnpNo].clear();
                
                allAncSnpIds[ancSnpId] = 1;
                
                ancSnpNo++;
            }
        }
    }
//    cout << "Total " << ancSnpNo << " SNPs in vcf. " << ancSnpNo << " SNPs have expected alleles.\n";
    
    DeleteAncSnpGtValues();
}

void VcfGrafAncSnpGeno::ShowSummary()
{
//    int numBadAncSnps = numBimAncSnps - numGoodAncSnps;
    
    string showSnpType = "RS IDs";
    if      (ancSnpType == AncestrySnpType::GB37) showSnpType = "GRCh 37 chromosome positions";
    else if (ancSnpType == AncestrySnpType::GB38) showSnpType = "GRCh 38 chromosome positions";
    
    cout << "\t" << showSnpType << " are used to find ancestry SNPs.\n";
//    cout << "\t" << numGoodAncSnps << " SNPs have expected alleles and will be used for ancestry inference.\n";
//    if (numDupAncSnps > 0) cout << "\t" << numDupAncSnps << " ancestry SNPs have multiple entries.\n";
    
//    if (numBadAncSnps > 0) {
//      cout << "\t" << numBadAncSnps << " ancestry SNPs do not have expected alleles.\n";
//    }
    cout << "\n";
}

void VcfGrafAncSnpGeno::CompareAncestrySnpAlleles(const string refStr, const string altsStr,
                                                  const char eRef, const char eAlt, int* expRefIdx, int* expAltIdx)
{
    // Let index of the ref alleles in vcf be 0, and indices of alts be 1, 2, 3, ...
    // Check them to find which one is the expected ref and and which is alt for the Ancestry SNP
    *expRefIdx = -1;
    *expAltIdx = -1;
    
    int refLen = refStr.length();
    if (refLen == 1) {
        char ref = refStr[0];
        char fRef = FlipAllele(ref);  // Flipped ref allele
        
        if (ref == eRef || fRef == eRef) {
            *expRefIdx = 0;
        }
        else if (ref == eAlt || fRef == eAlt) {
            *expAltIdx = 0;
        }

                
        // If ref seems to be flipped, then all the alts are flipped
        bool flip = false;
        if (fRef == eRef || fRef == eAlt) {
            flip = true;
        }
        
        if (*expRefIdx > -1 || *expAltIdx > -1) {
            vector<string> altWords = SplitString(altsStr, ",");
            int numAlts = altWords.size();
            
            for (int altNo = 0; altNo < numAlts; altNo++) {
                string altWord = altWords[altNo];
                int alleleIdx = altNo + 1;  // allele index starts from 0 = ref, then 1 = first alt, ...
                if (altWord.length() == 1) {
                    char alt = altWord[0];
                    if (flip)  alt = FlipAllele(alt);
                    
                    if      (alt == eRef && *expRefIdx < 0) *expRefIdx = alleleIdx;
                    else if (alt == eAlt && *expAltIdx < 0) *expAltIdx = alleleIdx;
                }
            }
        } // end if (*expRefIdx ...)
    }
}

int VcfGrafAncSnpGeno::RecodeGenotypeGivenIntegers(const int expRefIdx, const int expAltIdx, const int g1Num, const int g2Num)
{
    int genoInt = 3; // number of alt alleles, valid counts are 0, 1, 2
    
    // Genotype is valid only if both alleles are valid, i.e., same as one of the expected alleles
    int numValidGenos = 0;
    int numAlts = 0;
    
    if (g1Num == expRefIdx || g1Num == expAltIdx) {
        numValidGenos++;
        if (g1Num == expAltIdx) numAlts++;
    }
    
    if (g2Num == expRefIdx || g2Num == expAltIdx) {
        numValidGenos++;
        if (g2Num == expAltIdx) numAlts++;
    }
    
    if (numValidGenos == 2) genoInt = numAlts;
    
    return genoInt;
}
