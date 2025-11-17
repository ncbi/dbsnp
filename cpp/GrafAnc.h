using namespace std;

#ifndef GRAFANC_H
#define GRAFANC_H

#include <getopt.h>
#include "Util.h"
#include "AncestrySnps.h"
#include "VcfGrafAncSnpGeno.h"
#include "FamFileSamples.h"
#include "BimFileAncestrySnps.h"
#include "BedFileSnpGeno.h"
#include "SampleGenoDist.h"
#include "SampleGenoAncestry.h"
#include <thread>
#include <mutex>

#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "sys/types.h"
#include "sys/sysinfo.h"

#if defined(__sun)
#define PROC_SELF_EXE "/proc/self/path/a.out"
#else
#define PROC_SELF_EXE "/proc/self/exe"
#endif

string GetExecutablePath(void);
string FindFile(string);

int GetAvailableMemoryInMb();
int GetAllocatableMemoryInMb();

typedef struct parameters {
    int samplesPerRound; // Number of samples analyzed in each round 
    int numThreads;      // Maximum threads to use
    int maxMemoryInMb;   // Maximum memory to use
} parameters;

parameters GetParameters(int, char**);

#endif
