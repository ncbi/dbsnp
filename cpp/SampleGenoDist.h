#ifndef SAMPLE_GENO_DIST_H
#define SAMPLE_GENO_DIST_H
#include "Util.h"

// Given positions of the 3 vertices a, e, f and that of sample s:
// 1. Rotate all points so that the triangle is on the z=0 plane with f-a side parallel to x-axis
// 2. Calulate the Barycentric weights of the sample
class SampleGenoDist
{
private:
    // Genetic distances to the three reference populaton of the three vertex populations and the sample
    GenoDist aDist; // EAS
    GenoDist eDist; // EUR
    GenoDist fDist; // AFR
    GenoDist sDist; // Sample being checked

public:
    // The positions in the space of:the three vertex populations and the sample
    Point aPt;
    Point ePt;
    Point fPt;
    Point sPt;

    // The barycentric weights of the sample to the three vertices
    double aWt;
    double eWt;
    double fWt;

    // Position of fPt after transformation
    Point afrPosition;

    SampleGenoDist(GenoDist*, GenoDist*, GenoDist*, GenoDist*);

    void CopyGenoDist(GenoDist*, GenoDist*);
    void SetPointWithDist(Point*, GenoDist*);
    void CopyPoint(Point*, Point*);
    void MovePoint(Point*, double, double, double);
    void RotatePointOnx(Point*, double);
    void RotatePointOny(Point*, double);
    void RotatePointOnz(Point*, double);
    Point TransformGenoDist(GenoDist);
    void TransformAllDists();
    void CalculateBaryCenters();
    void ShowPositions(string, bool=false);
};

#endif
