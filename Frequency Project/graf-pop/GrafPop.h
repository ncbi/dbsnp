using namespace std;

#ifndef GRAFPOP_H
#define GRAFPOP_H

#include "Util.h"
#include "AncestrySnps.h"
#include "VcfSampleAncestrySnpGeno.h"
#include "FamFileSamples.h"
#include "BimFileAncestrySnps.h"
#include "BedFileSnpGeno.h"
#include "SampleGenoDist.h"
#include "SampleGenoAncestry.h"
#include <thread>
#include <mutex>

#if defined(__sun)
#define PROC_SELF_EXE "/proc/self/path/a.out"
#else
#define PROC_SELF_EXE "/proc/self/exe"
#endif

string GetExecutablePath(void);
string FindFile(string);

#endif
