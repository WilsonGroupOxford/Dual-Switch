#ifndef DUAL_SWITCH_GEOMETRYOPT_H
#define DUAL_SWITCH_GEOMETRYOPT_H

#include <map>
#include "crd.h"
#include "vectorManip.h"
#include "customIO.h"

class GeometryOptimiser{
//base class for geometry optimisers
protected:
    //optimisation parameters
    bool resolveStructure; //whether to attempt to resolve initial structure
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
    int checkInitialStructure(); //ensure starting structure is valid
    void calculateForces(); //from all contributions
    virtual void calculateBondForces()=0; //forces due to bonds
    virtual void calculateAngleForces()=0; //forces due to angles
    void lineSearch(); //find energy minimum along direction of force
    double calculateEnergy(); //of system
    virtual bool checkIntersections()=0; //check for line overlaps
    virtual double calculateBondEnergies()=0; //energies due to bonds
    virtual double calculateAngleEnergies()=0; //energies due to angles
    void checkConvergence(); //check if optimisation converged or hit limits

    //line overlap resolution
    bool resolveIntersections(); //attempt to resolve line overlaps
    virtual vector<DoublePair> getIntersections()=0;
    Pair findMajorIntersection(vector<DoublePair> &intersectingLines, int &majorCount);
    vector<int> findUniquePoints(vector<DoublePair> &intersectingLines, Pair &majorLine);
    virtual bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections)=0;

public:
    //constructors
    GeometryOptimiser();

    //setters
    void setCoordinates(vector<Crd2d> crds);
    void setOptisationParameters(double convCriteria, double lineInc, int maxIt, bool res);
    void setSystemParameters(vector<Pair> sysBonds, vector<Trio> sysAngles, vector<int> sysFixed, vector<DoublePair> sysInt);
    //set potential parameters in concrete class

    //getters
    vector<Crd2d> getMinimisedCoordinates(); //return minimised coordinates
    double getEnergy(); //return final energy
    int getIterations(); //number of iterations before optimisation completed

    //optimisation functions
    virtual int steepestDescent()=0; //main geometry optimisation function
};

class AperiodicGeometryOptimiser: public GeometryOptimiser{
//base class for aperiodic geometry optimisers
protected:
    //virtual functions to define
    bool checkIntersections() override;
    vector<DoublePair> getIntersections() override;
    bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections) override;

public:
    //constructor
    AperiodicGeometryOptimiser();

    //optimisation functions
    int steepestDescent() override; //main geometry optimisation function
};

class PeriodicGeometryOptimiser: public GeometryOptimiser{
//base class for periodic geometry optimisers
protected:
    //virtual functions to define
    bool checkIntersections() override;
    vector<DoublePair> getIntersections() override;
    bool moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections) override;

    //additional functions
    void wrapCoordinates(); //apply periodic boundary to coordinates

    //additional variables
    double pbcX, pbcY, pbcRX, pbcRY; //periodic boundary and reciprocals

public:
    //constructor
    PeriodicGeometryOptimiser();

    //setters
    void setPeriodicBoundary(double x, double y, double rx, double ry);

    //optimisation functions
    int steepestDescent() override; //main geometry optimisation function
};

class HarmonicAperiodicGO: public AperiodicGeometryOptimiser {
//harmonic, aperiodic
private:
    //virtual functions to define
    void calculateBondForces() override;
    void calculateAngleForces() override;
    double calculateBondEnergies() override;
    double calculateAngleEnergies() override;

    //additional functions
    Crd2d bondForce(Crd2d &cI, Crd2d &cJ, double &r0);

    //potential parameters
    double forceK, energyK; //force constant and energy constant
    vector<double> r0; //potential minimum

public:
    HarmonicAperiodicGO();
    void setPotentialParameters(double harmK, vector<double> harmR0);
};

class HarmonicPeriodicGO: public PeriodicGeometryOptimiser {
//harmonic, periodic
private:
    //virtual functions to define
    void calculateBondForces() override;
    void calculateAngleForces() override;
    double calculateBondEnergies() override;
    double calculateAngleEnergies() override;

    //additional functions
    Crd2d bondForce(Crd2d &cI, Crd2d &cJ, double &r0);

    //potential parameters
    double forceK, energyK; //force constant and energy constant
    vector<double> r0; //potential minimum

public:
    HarmonicPeriodicGO();
    void setPotentialParameters(double harmK, vector<double> harmR0);
};

class KeatingAperiodicGO: public AperiodicGeometryOptimiser {
//keating, aperiodic
private:
    //virtual functions to define
    void calculateBondForces() override;
    void calculateAngleForces() override;
    double calculateBondEnergies() override;
    double calculateAngleEnergies() override;

    //additional functions
    Crd2d bondForce(Crd2d &cI, Crd2d &cJ);
    void angleForce(Crd2d &cI, Crd2d &cJ, Crd2d &cK, Crd2d &fI, Crd2d &fJ, Crd2d &fK);

    //potential parameters
    double aSq, aSq_2, kBondForce, kAngleForce, kBondEnergy, kAngleEnergy;

public:
    KeatingAperiodicGO();
    void setPotentialParameters(double a, double alpha, double beta);
};

class KeatingPeriodicGO: public PeriodicGeometryOptimiser {
//keating, periodic
private:
    //virtual functions to define
    void calculateBondForces() override;
    void calculateAngleForces() override;
    double calculateBondEnergies() override;
    double calculateAngleEnergies() override;

    //additional functions
    Crd2d bondForce(Crd2d &cI, Crd2d &cJ);
    void angleForce(Crd2d &cI, Crd2d &cJ, Crd2d &cK, Crd2d &fI, Crd2d &fJ, Crd2d &fK);

    //potential parameters
    double aSq, aSq_2, kBondForce, kAngleForce, kBondEnergy, kAngleEnergy;

public:
    KeatingPeriodicGO();
    void setPotentialParameters(double a, double alpha, double beta);
};

class KeatingAperiodicBondOnlyGO: public AperiodicGeometryOptimiser {
//keating, aperiodic, no angle term
private:
    //virtual functions to define
    void calculateBondForces() override;
    void calculateAngleForces() override;
    double calculateBondEnergies() override;
    double calculateAngleEnergies() override;

    //additional functions
    Crd2d bondForce(Crd2d &cI, Crd2d &cJ, int bond);

    //potential parameters
    vector<double> aSq, aSq_2, kBondForce, kBondEnergy;

public:
    KeatingAperiodicBondOnlyGO();
    void setPotentialParameters(double alpha, vector<double> a);
};

#endif //DUAL_SWITCH_GEOMETRYOPT_H
