#ifndef DUAL_SWITCH_GEOMETRYOPT_H
#define DUAL_SWITCH_GEOMETRYOPT_H

#include <map>
#include "crd.h"
#include "vectorManip.h"
#include "customIO.h"

class GeometryOptimiser {
//base class for optimisers
protected:
    //optimisation parameters
    int maxIterations; //maximum possible iterations for optimisation
    double convergenceCriteria, lineSearchIncrement; //for calculation end and optimisation line search

    //potential parameters

    //system parameters
    int nPnts, nBonds, nAngles, nFixed, nIntersectionLines; //number of points, bonds, angles, fixed and intersecting lines
    vector<Crd2d> coordinates; //list of x,y coordinates
    vector<Pair> bonds; //points forming bonds
    vector<Trio> angles; //points forming angles
    vector<int> fixed; //points fixed in place
    vector<DoublePair> intersectionLines; //lines forming intersections

    //calculation variables
    bool complete, converged; //optimisation status
    int iterations, deltaEZeroCount; //of optimisation, no energy change counter
    double previousEnergy, currentEnergy; //energies of previous and current iteration
    Crd2d zeroForce; //single force of zero
    vector<Crd2d> zeroForces; //vector where all forces set to zero
    vector<Crd2d> forces; //force for each coordinate

    //steepest descent functions
    virtual int checkInitialStructure()=0; //ensure starting structure is valid
    void calculateForces(); //from all contributions
    virtual void calculateBondForces()=0; //forces due to bonds
    virtual void calculateAngleForces()=0; //forces due to angles
    void lineSearch(); //find energy minimum along direction of force
    double calculateEnergy(); //of system
    virtual bool checkIntersections()=0; //check for line overlaps
    virtual double calculateBondEnergies()=0; //energies due to bonds
    virtual double calculateAngleEnergies()=0; //energies due to angles
    void checkConvergence(); //check if optimisation converged or hit limits

public:
    //constructors
    GeometryOptimiser();

    //setters
    void setCoordinates(vector<Crd2d> crds);
    void setOptisationParameters(double convCriteria, double lineInc, int maxIt);
    void setSystemParameters(vector<Pair> sysBonds, vector<Trio> sysAngles, vector<int> sysFixed, vector<DoublePair> sysInt);
    //set potential parameters in concrete class

    //getters
    vector<Crd2d> getMinimisedCoordinates(); //return minimised coordinates
    double getEnergy(); //return final energy
    int getIterations(); //number of iterations before optimisation completed

    //optimisation functions
    int steepestDescent(); //main geometry optimisation function
};

class IntersectionResolver {
    //base class for finding and resolving line intersections
protected:
    virtual vector<DoublePair> getIntersections()=0;
    Pair findMajorIntersection(vector<DoublePair> &intersectingLines, int &majorCount);
    vector<int> findUniquePoints(vector<DoublePair> &intersectingLines, Pair &majorLine);
    virtual bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections)=0;

public:
    IntersectionResolver(); //default constructor
    bool resolveIntersections(); //attempt to resolve line intersections
};

class HarmonicAperiodicGO: public GeometryOptimiser, public IntersectionResolver {
private:
    //virtual functions to define
    int checkInitialStructure() override;
    void calculateBondForces() override;
    void calculateAngleForces() override;
    bool checkIntersections() override;
    double calculateBondEnergies() override;
    double calculateAngleEnergies() override;
    vector<DoublePair> getIntersections() override;
    bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections) override;

    //additional functions
    Crd2d bondForce(Crd2d &cI, Crd2d &cJ, double &r0);

    //potential parameters
    double forceK, energyK; //force constant and energy constant
    vector<double> r0; //potential minimum

public:
    HarmonicAperiodicGO();
    void setPotentialParameters(double harmK, vector<double> harmR0);
};


#endif //DUAL_SWITCH_GEOMETRYOPT_H
