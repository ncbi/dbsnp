#include "Util.h"

bool FileExists (const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (fp) {
        fclose(fp);
        return true;
    }
    else {
        return false;
    }
}

bool FileWriteable (const char* filename)
{
    FILE* fp = fopen(filename, "a");
    if (fp) {
        fclose(fp);
        return true;
    }
    else {
        return false;
    }
}

string LowerString(const string& inStr)
{
    int strLen = inStr.length();
    string outStr = "";
    for (int i = 0; i < strLen; i++) {
        char lowCh = tolower(inStr[i]);
        outStr.push_back(lowCh);
    }

    return outStr;
}

string UpperString(const string& inStr)
{
    int strLen = inStr.length();
    string outStr = "";
    for (int i = 0; i < strLen; i++) {
        char upCh = toupper(inStr[i]);
        outStr.push_back(upCh);
    }

    return outStr;
}

GenoDatasetType CheckGenoDataFile(const string& file, string *baseName)
{
    GenoDatasetType fileType = GenoDatasetType::IS_OTHER;

    // Check if this is a plink set basename
    bool plinkExists = false;
    string bedFile = file + ".bed";
    string bimFile = file + ".bim";
    string famFile = file + ".fam";

    if (FileExists(bedFile.c_str()) && FileExists(bimFile.c_str()) && FileExists(famFile.c_str())) {
        *baseName = file;
        fileType = GenoDatasetType::IS_PLINK;
        return fileType;
    }

    // Check if file exists
    bool fileExists = false;
    if (FileExists(file.c_str())) {
        fileExists = true;
    }
    else {
        fileType = GenoDatasetType::NOT_EXISTS;
    }

    // Check if it is a .gz file
    int fileLen = file.length();
    string gzFileBase = file; // File name. If it is .gz, stripped off ".gz"
    bool isGz = false;
    if (fileLen > 3 && file.substr(fileLen-3, 3).compare(".gz") == 0) {
        isGz = true;
        gzFileBase = file.substr(0, fileLen-3);
        fileLen -= 3;
    }

    // Check if it is a vcf or PLINK
    string fileBase = gzFileBase;
    string fileExt = "";
    size_t dotPos = gzFileBase.find_last_of(".");

    if (dotPos != string::npos) {
        fileExt = gzFileBase.substr(dotPos+1, fileLen-dotPos-1);
        fileBase = gzFileBase.substr(0, dotPos);
    }

    bool isVcf = false;
    bool isPlink = false;
    if (fileExt.compare("vcf") == 0) {
        isVcf = true;
    }
    else if (fileExt.compare("bed") == 0 ||
             fileExt.compare("bim") == 0 ||
             fileExt.compare("fam") == 0   ) {
        isPlink = true;
    }
    else {
        fileBase = gzFileBase;
    }

    if (fileExists) {
        if      (isPlink && !isGz) fileType = GenoDatasetType::IS_PLINK;
        else if (isPlink &&  isGz) fileType = GenoDatasetType::IS_PLINK_GZ;
        else if (isVcf   && !isGz) fileType = GenoDatasetType::IS_VCF;
        else if (isVcf   &&  isGz) fileType = GenoDatasetType::IS_VCF_GZ;
    }

    *baseName = fileBase;

    return fileType;
}

int GetChromosomeFromString(const char* chrStr)
{
    int chrNum = 0;
    int i = 0;

    if (strlen(chrStr) > 3 &&
        (chrStr[0] == 'c' || chrStr[0] == 'C') &&
        (chrStr[1] == 'h' || chrStr[1] == 'H') &&
        (chrStr[2] == 'r' || chrStr[2] == 'R') ) {
        i = 3;
    }

    while(chrStr[i] != 0) {
        int num = chrStr[i] - '0';

        if (num >= 0 && num < 10) {
            chrNum = chrNum * 10 + num;
        }
        else {
            chrNum = 0;
            break;
        }

        i++;
    }

    if (chrNum < 1 || chrNum > 22) chrNum = 0;

    return chrNum;
}

int GetRsNumFromString(const char* rsStr)
{
    int rsNum = 0;

    if (strlen(rsStr) > 2 &&
        (rsStr[0] == 'r' || rsStr[0] == 'R') &&
        (rsStr[1] == 's' || rsStr[1] == 'S') ) {
        int i = 2;

        while(rsStr[i] != 0) {
            int num = rsStr[i] - '0';

            if (num >= 0 && num < 10) {
                rsNum = rsNum * 10 + num;
            }
            else {
                rsNum = 0;
                break;
            }

            i++;
        }
    }

    return rsNum;
}

char FlipAllele(char allele)
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

vector<string> SplitString(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;

    while (pos < str.length() && prev < str.length()) {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }

    return tokens;
}

string TrimString(const string &str) 
{
    string newStr = str;

    newStr.erase(0, newStr.find_first_not_of("\t\n\v\f\r "));
    newStr.erase(newStr.find_last_not_of("\t\n\v\f\r ") + 1);

    return newStr;
}

void ShowTimeDiff(const struct timeval &t1, const struct timeval &t2)
{
    int usec = t2.tv_usec - t1.tv_usec;
    int sec = t2.tv_sec - t1.tv_sec;
    if (usec < 0) {
        usec += 1000000;
        sec += 1;
    }
    int min = 0;
    if (sec > 60) {
        min = sec / 60;
        sec = sec % 60;
    }
    int hour = 0;
    if (min > 60) {
        hour = min / 60;
        min = min % 60;
    }

    printf("Time used: ");
    if (hour > 0) { printf("%d hours ", hour); }
    if (min > 0) { printf("%d minutes ", min); }
    if (sec > 0) { printf("%d seconds ", sec); }
    printf("%d microseconds\n\n", usec);
}
