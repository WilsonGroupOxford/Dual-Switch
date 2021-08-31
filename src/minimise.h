//lightweight minimiser

#ifndef DUAL_SWITCH_MINIMISE_H
#define DUAL_SWITCH_MINIMISE_H

#include <map>
#include "crd.h"
#include "vectorManip.h"
#include "customIO.h"

using namespace std;

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

    //overlap resolution functions
    bool resolveInitialIntersections(); //attempt to remove starting intersections
    vector<DoublePair> getIntersections(); //get intersecting lines
    Pair findMajorIntersection(vector<DoublePair> &intersectingLines, int &majorCount); //find intersection common to all
    vector<int> findUniquePoints(vector<DoublePair> &intersectingLines, Pair &majorLine); //get list of points which could be moved
    bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections); //move subset of points to resolve intersections

    //periodic minimisation functions
    void calculateForcesPBC(); //force on each point with pbc
    double calculateEnergyPBC(); //energy of configuration with pbc
    Crd2d harmonicForcePBC(Crd2d &c1, Crd2d &c2, double &r0, double &k); //force of harmonic crds with pbc
    double harmonicEnergyPBC(Crd2d &c1, Crd2d &c2, double &r0, double &k); //energy of harmonic crds with pbc
    void lineSearchPBC(); //move along lines of force to minimum energy with pbc
    bool checkIntersectionsPBC(); //check if any line intersections
    void wrapAroundCoordinates(); //apply periodic boundary to wrap coordinates

    //periodic overlap resolution functions
    bool resolveInitialIntersectionsPBC(); //attempt to remove starting intersections
    vector<DoublePair> getIntersectionsPBC(); //get intersecting lines
    bool moveIntersectingPointsPBC(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections); //move subset of points to resolve intersections

    //###### tests ######
    void runTests();

public:
    HarmonicMinimiser(); //default constructor invokes tests
    HarmonicMinimiser(vector<Crd2d> crds, vector<Pair> harmPairs, vector<int> fixedPnts,
                      vector<double> harmR0, double harmK, double cc, double inc, int maxIt, vector<DoublePair> lineInt);
    int steepestDescent(); //steepest descent minimisation
    int steepestDescent(double x, double y, double rx, double ry); //steepest descent minimisation with periodic boundary conditions
    vector<Crd2d> getMinimisedCoordinates(); //return minimised coordinates
    double getEnergy(); //return final energy
    int getIterations(); //number of iterations before optimisation completed
};

class KeatingMinimiser {
private:
//     //potential model
    int nBonds, nAngles, nIntersectionPairs;
//    vector<int> fixedPoints; //points fixed in place
    vector<Pair> bonds; //points forming bonds
    vector<Trio> angles; //points forming angles
    vector<DoublePair> intersectionPairs; //lines to check for intersection
//    vector<double> harmonicR0; //minimum in harmonic potential
//    double harmonicK; //force constant for harmonic potential
    double convergenceCriteria, lineSearchIncrement;
    int maxIterations;
    double a, alpha, beta; //bond/angle params
    double aSq, aSq_2; //square and half of square of a
    double kBondForce, kAngleForce, kBondEnergy, kAngleEnergy; //constants for foce + angle calculations


    //calculation variables
    int nPnts; //number of points
    vector<Crd2d> coordinates, forces;
//
//    //periodic boundary condition functions
//    double pbcX, pbcY, pbcRX, pbcRY; //periodic boundary conditions and reciprocals
//
    //minimisation variables
    bool complete, converged; //flags for minimisation completion status
    int iterations, deltaEZeroCount; //number of minimisation iterations
    double previousEnergy, currentEnergy; //energy of previous iteration
    Crd2d zeroForce; //force of zero
    vector<Crd2d> zeroForces; //vector where all forces set to zero
//
//    //minimisation functions
    void calculateForces(); //force on each point
    double calculateEnergy(); //energy of configuration
    Crd2d bondForce(Crd2d &c1, Crd2d &c2); //force from bonds
    void angleForce(Crd2d &cI, Crd2d &cJ, Crd2d &cK, Crd2d &fI, Crd2d &fJ, Crd2d &fK); //force from angles
    double bondEnergy(Crd2d &c1, Crd2d &c2); //energy of bonds
    double angleEnergy(Crd2d &cI, Crd2d &cJ, Crd2d &cK); //energy of angles
    void lineSearch(); //move along lines of force to minimum energy
    bool checkIntersections(); //check if any line intersections
    void checkConvergence(); //check energy convergence or maximum iterations
//
//    //overlap resolution functions
//    bool resolveInitialIntersections(); //attempt to remove starting intersections
//    vector<DoublePair> getIntersections(); //get intersecting lines
//    Pair findMajorIntersection(vector<DoublePair> &intersectingLines, int &majorCount); //find intersection common to all
//    vector<int> findUniquePoints(vector<DoublePair> &intersectingLines, Pair &majorLine); //get list of points which could be moved
//    bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections); //move subset of points to resolve intersections
//
//    //periodic minimisation functions
//    void calculateForcesPBC(); //force on each point with pbc
//    double calculateEnergyPBC(); //energy of configuration with pbc
//    Crd2d harmonicForcePBC(Crd2d &c1, Crd2d &c2, double &r0, double &k); //force of harmonic crds with pbc
//    double harmonicEnergyPBC(Crd2d &c1, Crd2d &c2, double &r0, double &k); //energy of harmonic crds with pbc
//    void lineSearchPBC(); //move along lines of force to minimum energy with pbc
//    bool checkIntersectionsPBC(); //check if any line intersections
//    void wrapAroundCoordinates(); //apply periodic boundary to wrap coordinates
//
//    //periodic overlap resolution functions
//    bool resolveInitialIntersectionsPBC(); //attempt to remove starting intersections
//    vector<DoublePair> getIntersectionsPBC(); //get intersecting lines
//    bool moveIntersectingPointsPBC(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections); //move subset of points to resolve intersections

public:
    KeatingMinimiser(vector<Crd2d> crds, double cc, double inc, int maxIt); //coordinates and search parameters
    void setParameters(double kpA, double kpAlpha, double kpBeta);
    void setInteractions(vector<Pair> kpBonds, vector<Trio> kpAngles, vector<DoublePair> ints); //lists of atoms for potential
    int steepestDescent(); //steepest descent minimisation
    int steepestDescent(double x, double y, double rx, double ry); //steepest descent minimisation with periodic boundary conditions
    vector<Crd2d> getMinimisedCoordinates(); //return minimised coordinates
    double getEnergy(); //return final energy
    int getIterations(); //number of iterations before optimisation completed
};


#endif //DUAL_SWITCH_MINIMISE_H
