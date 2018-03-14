//lightweight minimiser

#ifndef DUAL_SWITCH_MINIMISE_H
#define DUAL_SWITCH_MINIMISE_H

#include "crd.h"

class HarmonicMinimiser {
private:
    //potential model
    int nHarmonicPairs, nFixedPnts, nIntersectionPairs;
    vector<int> fixedPoints; //points fixed in place
    vector<Pair> harmonicPairs; //points linked by harmonic pairs
    vector<double> harmonicR0; //minimum in harmonic potential
    double harmonicK; //force constant for harmonic potential
    double convergenceCriteria, lineSearchIncrement;
    int maxIterations;
    vector<DoublePair> intersectionPairs; //lines to check for intersection

    //calculation variables
    int nPnts; //number of points
    vector<Crd2d> coordinates, forces;

    //periodic boundary condition functions
    double pbcX, pbcY, pbcRX, pbcRY; //periodic boundary conditions and reciprocals

    //minimisation variables
    bool complete, converged; //flags for minimisation completion status
    int iterations, deltaEZeroCount; //number of minimisation iterations
    double previousEnergy, currentEnergy; //energy of previous iteration
    Crd2d zeroForce; //force of zero
    vector<Crd2d> zeroForces; //vector where all forces set to zero

    //minimisation functions
    void calculateForces(); //force on each point
    double calculateEnergy(); //energy of configuration
    Crd2d harmonicForce(Crd2d &c1, Crd2d &c2, double &r0, double &k); //force of harmonic crds
    double harmonicEnergy(Crd2d &c1, Crd2d &c2, double &r0, double &k); //energy of harmonic crds
    void lineSearch(); //move along lines of force to minimum energy
    bool checkIntersections(); //check if any line intersections
    void checkConvergence(); //check energy convergence or maximum iterations

    //periodic minimisation functions
    void calculateForcesPBC(); //force on each point with pbc
    double calculateEnergyPBC(); //energy of configuration with pbc
    Crd2d harmonicForcePBC(Crd2d &c1, Crd2d &c2, double &r0, double &k); //force of harmonic crds with pbc
    double harmonicEnergyPBC(Crd2d &c1, Crd2d &c2, double &r0, double &k); //energy of harmonic crds with pbc
    void lineSearchPBC(); //move along lines of force to minimum energy with pbc
    bool checkIntersectionsPBC(); //check if any line intersections
    void wrapAroundCoordinates(); //apply periodic boundary to wrap coordinates
    //###### tests ######
    void runTests();

public:
    HarmonicMinimiser(); //default constructor invokes tests
    HarmonicMinimiser(vector<Crd2d> crds, vector<Pair> harmPairs, vector<int> fixedPnts,
                      vector<double> harmR0, double harmK, double cc, double inc, int maxIt, vector<DoublePair> lineInt);
    void steepestDescent(); //steepest descent minimisation
    void steepestDescent(double x, double y, double rx, double ry); //steepest descent minimisation with periodic boundary conditions
    vector<Crd2d> getMinimisedCoordinates(); //return minimised coordinates
};



#endif //DUAL_SWITCH_MINIMISE_H
