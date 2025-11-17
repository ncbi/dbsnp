BASEDIR = ./
INCLUDES = -I. -I$(BASEDIR)

#------ Compiler and options -----------------
CXX = /usr/bin/g++
HTSLIB = ~/htslib/1.20
CXXFLAGS = -std=c++11 -pthread -g -lm -lz $(INCLUDES) -I $(HTSLIB)/include -L $(HTSLIB)/lib -lhts -lstdc++

HDIR = ./
SRCDIR = ./

#----- Suffix Rules ---------------------------
.SUFFIXES: .cpp .C .cc

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

#----- File Dependencies ----------------------

SRC = Util.cpp AncestrySnps.cpp VcfGrafAncSnpGeno.cpp FamFileSamples.cpp BimFileAncestrySnps.cpp BedFileSnpGeno.cpp SampleGenoDist.cpp SampleGenoAncestry.cpp  GrafAnc.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

grafanc: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

Util.o: $(HDIR)Util.h
	$(CXX) $(CXXFLAGS) -c Util.cpp
AncestrySnps.o: $(HDIR)AncestrySnps.h
	$(CXX) $(CXXFLAGS) -c AncestrySnps.cpp
VcfGrafAncSnpGeno.o: $(HDIR)VcfGrafAncSnpGeno.h
	$(CXX) $(CXXFLAGS) -c VcfGrafAncSnpGeno.cpp
FamFileSamples.o: $(HDIR)FamFileSamples.h
	$(CXX) $(CXXFLAGS) -c FamFileSamples.cpp
BimFileAncestrySnps.o: $(HDIR)BimFileAncestrySnps.h
	$(CXX) $(CXXFLAGS) -c BimFileAncestrySnps.cpp
BedFileSnpGeno.o: $(HDIR)BedFileSnpGeno.h
	$(CXX) $(CXXFLAGS) -c BedFileSnpGeno.cpp
SampleGenoDist.o: $(HDIR)SampleGenoDist.h
	$(CXX) $(CXXFLAGS) -c SampleGenoDist.cpp
SampleGenoAncestry.o:$(HDIR)SampleGenoAncestry.h
	$(CXX) $(CXXFLAGS) -c SampleGenoAncestry.cpp

depend:
	makedepend $(CXXFLAGS) -Y $(SRC)

clean:
	rm -f $(OBJ) *~

