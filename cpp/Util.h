#ifndef NDEBUG
#   define ASSERT(condition, message) \
    if (!(condition)) {	\
        cerr << "Assertion `" #condition "` failed in " << __FILE__  << " line " << __LINE__ << ": " << message << "\n"; \
        terminate(); \
    }
#else
#   define ASSERT(condition, message) {}
#endif

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <map>
#include <vector>
#include <unistd.h>
#include <fstream>
  
const double pi = 3.1415926;

using namespace std;

enum class AncestrySnpType
{
    RSID = 0,
    GB37 = 1,
    GB38 = 2
};

enum class GenoDatasetType
{
    NOT_EXISTS = 0,
    IS_PLINK = 1,
    IS_PLINK_GZ = 2,
    IS_VCF = 3,
    IS_VCF_GZ = 4,
    IS_OTHER = 5
};

// Define Genetic Distances to the three reference populations
struct GenoDist
{
    double e; // To European
    double f; // To African
    double a; // To East Asian
};

// A 3-D point in space
struct Point
{
    double x;
    double y;
    double z;
};

bool FileExists (const char*);
bool FileWriteable (const char*);
int GetChromosomeFromString(const char*);
int GetRsNumFromString(const char*);
void ShowTimeDiff(const struct timeval&, const struct timeval&);
char FlipAllele(char);
vector<string> SplitString(const string&, const string&);
string TrimString(const string &str); 
string LowerString(const string&);
string UpperString(const string&);
GenoDatasetType CheckGenoDataFile(const string&, string*);


#endif
