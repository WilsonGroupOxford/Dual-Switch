//network of single random seed

#ifndef DUAL_SWITCH_NETWORK_H
#define DUAL_SWITCH_NETWORK_H

#include <map>
#include <random>
#include "customIO.h"
#include "node.h"
#include "vectorManip.h"
#include "geometryopt.h"

using namespace std;

class Network {
private:

    //###### Input variables ######

    //IO variables
    string inPrefix, outPrefix; //for input and output files

    //Network properties
    bool periodic, load, overridePbc; //periodic network, load existing network, override default pbc
    double targetAlpha; //target for aboav-weaire alpha parameter
    string initialLatticeType; //starting lattice for mc
    vector<int> initialLatticeDimensions, ringSizeLimits; //starting lattice size, max/min ring size values
    vector<double> targetPVector; //target ring statistics
    vector<double> overridelatticePbc; //overriden periodc boundary conditions

    //Potential model
    bool localGeomOpt, globalGeomOpt; //geometry optimisation local after each accepted step of global at end
    int localGeomOptMaxIt, globalGeomOptMaxIt; //maximum number of steepest descent iterations
    double localGeomOptCC, globalGeomOptCC; //convergence criteria
    double atomicSeparation, geomOptLineSearchInc; //ideal distance between atoms, line search increment
    double harmonicK; //force constant

    //Monte carlo dual switching properties
    int mcSeed, mcMaxMoves; //random seed, maximum proposed moves
    double mcTemperature, mcConvergence, mcAlphaScaleFactor; //temperature, properties convergence criteria, scale factor for energy functional

    //Analysis variables
    bool convertToAtomic, periodicVisualisation, spatialRdf, topoRdf, assortativeMix, areaLaw, clustering; //convert dual to atomic network, make periodic image, spatial partial rdfs, topological rdfs, assortative mixing, ring areas, clustering
    double spatialRdfBinwidth, spatialRdfExtent, topoRdfExtent; //size of bin rdf, max distance for rdf

    //Atomic Geometry Optimisation
    bool atomicGeomOpt; //geometry optimisation of atomic network at end
    double keatingA, keatingAlpha, keatingBeta; //keating potential parameters
    int atomicGeomOptMaxIt; //maximum number of steepest descent iterations
    double atomicGeomOptCC, atomicLineSearchInc; //convergence criteria, line search increment

    //###### Construction variables ######

    //Further network properties
    bool consistent, noIntersections; //network and p vector/matrix consistent, check for edge intersections
    string name; //name of network for writing to files
    int minRingSize, maxRingSize, nRingSizes, index6, nNodes, nRings; //ring size limits, index of 6mr, number of nodes in network, number of rings in network
    double periodicBoxX, periodicBoxY, rPeriodicBoxX, rPeriodicBoxY; //dimensions of periodic box and reciprocals
    vector<Node> nodes; //nodes in network
    vector<int> pVector; //number of nodes of each size
    vector< vector<int> > pMatrix; //number of nodes adjacent to nodes of each size
    vector<double> aboavWeaireParams; //alpha, mu, rsq
    vector<Ring> nodeRings; //rings in dual graph

    //Further potential model variables
    vector< vector<double> > harmonicR0Matrix; //harmonic minimum for nodes of different sizes

    //Further monte carlo properties
    bool mcTargetReached; //whether targets met within convergence
    int testCounter, geomOptRejectCount, mcProposedMoves, mcAcceptedMoves; //number of moves proposed and accepted during monte carlo
    double mcEnergy, targetMu; //energy of system relative to target properties, mu
    double rMcTemperature, rTargetMu; //reciprocals of key values
    vector<double> rTargetPVector; //reciprocals of key values
    mt19937 randomGenerator; //generator for distributions
    uniform_int_distribution<int> nodeDistribution; //mersenne twister uniform distribution between 0->nNodes
    uniform_real_distribution<double> zeroOneDistribution, connectionPickDistribution; //mersenne twister uniform distribution between 0->1

    //further geometry optimisation variables
    double geomOptEnergy; //energy of global geometry optimisation
    int geomOptIterations; //number of iterations of global geometry optimisation

    //###### Construction functions ######

    //Initialisation
    void initialiseNetworkProperties(); //additional network properites
    void initialisePotentialModel(); //make parameters for potential model
    void initialisePeriodicLattice(); //make initial periodic hexagonal lattice
    void initialiseAperiodicLattice(); //make initial aperiodic hexagonal lattice
    void initialisePeriodicCrystalLattice(string crystal); //make alternative starting crystal lattice
    void initialiseAperiodicCrystalLattice(); //XXX
    void hexPeriodicCrystal(int layer, int offset, int blockLength, int blockSkip, bool reverse); //defined moves in hexagonal symmetry
    void cubicPeriodicCrystal(int layer, int offset, int blockLength, int blockSkip); //defined moves in cubic symmetry
    void initialiseMonteCarlo(); //set up variables and initialise random number generators
    void loadPeriodicLattice(); //read in existing periodic lattice

    //Random numbers
    int pickRandomNode(), pickRandomConnection(int nCnxs); //select node from network, select connection
    double metropolisRandomNum(); //random number for mc metropolis condition

    //Monte carlo
    void monteCarlo(); //main monte carlo process
    void definedMove(vector<int> &switchTriangles); //perform defined dual switch move
    vector<int> pickRandomTrianglePairPeriodic(); //nodes to dual switch
    vector<int> pickRandomTrianglePairAperiodic(); //nodes to dual switch
    vector<int> pickDefinedTrianglePair(int m); //nodes to dual switch - for development
    void calculateTrialPPeriodic(vector<int> &triangles, vector<int> &trialPVector, vector< vector<int> > &trialPMatrix); //calculate trial p vector and matrix
    void calculateTrialPAperiodic(vector<int> &triangles, vector<int> &trialPVector, vector< vector<int> > &trialPMatrix); //calculate trial p vector and matrix
    vector<double> calculateAboavWeaireFit(vector<int> &pVec, vector< vector<int> > &pMat); //calculate Aboav-Weaire parameters
    double mcEnergyFunctional(vector<double> &awParams, vector<int> &pVec); //calculate energy for mc
    bool evaluateMetropolisCondition(double &trialEnergy, double &currEnergy); //accept or reject mc move
    bool acceptDualSwitchAperiodic(vector<int> &switchTriangles, vector<int> &trialPVec, vector< vector<int> > &trialPMat, double &trialMcEnergy, vector<double> &trialAwParams); //enact dual switch and update trial->current variables
    bool acceptDualSwitchPeriodic(vector<int> &switchTriangles, vector<int> &trialPVec, vector< vector<int> > &trialPMat, double &trialMcEnergy, vector<double> &trialAwParams); //enact dual switch and update trial->current variables
    void findNodeRings(); //find rings of nodes
    int localMinimisationPeriodic(vector<int> &switchTriangles); //minimise locally
    int localMinimisationAperiodic(vector<int> &switchTriangles); //minimise locally
    void globalMinimisationPeriodic(); //minimise globally
    void globalMinimisationAperiodic(); //minimise globally
    vector<int> nextTopologicalShell(vector<int> &currentShell); //find two-sided topological shell
    vector<int> nextTopologicalShell(vector<int> &currentShell, vector<int> &prevShell); //find two-sided topological shell
    vector<int> findNodeRing(int n); //find ring around a central node

    //Checking
    void checkFidelity(); //check to ensure consistency
    void checkGeometry(); //check for node edge overlap

    //###### Analysis variables ######
    int nVertices; //number of vertices of atomic graph
    vector<Vertex> vertices; //vertices of atomic graph
    vector<Ring> vertexRings; //rings of atoms
    vector<double> ringStatistics; //ring statistics of network
    vector< vector<double> > piMatrix; //normalised p matrix
    vector<double> spatialRdfDensities; //density of nodes
    vector<Rdf> spatialPartialRdfs; //rdfs for each rings size pair
    vector< vector<int> > topoRdfShellSizes; //number of rings in each shell for each ring size
    vector<Rdf> topoPartialRdfs; //rdfs for each ring size pair
    vector<double> assortativeMixing; //pearson cc via two methods
    vector<double> ringAreas; //areas of rings of different sizes
    int atomicGeomOptStatus, atomicGeomOptIterations; //whether optimised, number of iterations
    double atomicGeomOptEnergy; //final energy of geometry optimisation
    vector<double> atomicGeomBondMean, atomicGeomAngleMean; //mean of atomic bond and angles for each ringsize
    vector<double> atomicGeomBondStd, atomicGeomAngleStd; //standard deviation of atomic bond and angles for each ringsize
    vector<double> atomicBondDistribution, atomicAngleDistribution; //full bond and angle distribution, listed in order of ring size
    vector< vector<int> > clusterValues, clusterDistributions; //sizes of clusters

    //###### Analysis functions ######
    void convertDualToAtomicNetwork(); //triangulate nodes to get atomic network
    void analyseRingStatistics(); //calculate ring statistics and normalise pi matrix
    void analyseRingAreas(); //calculate ring areas
    void analyseAboavWeaire(); //calculate aboav-weaire parameters
    void analysePartialSpatialRdfs(); //calculate partial rdfs for each ring size
    void analysePartialTopologicalRdfs(); //calculate partial topological rdfs for each ring size
    void analyseAssortativity(); //calculate assortativity
    void analyseClusters(); //analyse size of clusters
    void analyseAtomicGeometry(); //calculate average/std of bond lengths and angles
    void geometryOptimiseAtomicNetworkPeriodic(); //globally minimise periodic atomic network
    void geometryOptimiseAtomicNetworkAperiodic(); //globally minimise aperiodic atomic network
    void meanAndStdDeviation(vector<double> &values, double &mean, double &stdev); //calculates mean and std of values

    //###### Write functions ######
    void writeDual(); //write out dual network
    void writeAtomicNetwork(); //write out atomic network
    void writePeriodicDualNetwork(); //calculate and write out network for periodic visualisation
    void writePeriodicAtomicNetwork(vector<Node> &periodicNodes); //calculate and write out atomic network for visualisation
    void writeGeometryOptimisationEnergy(); //write energy of system
    void writeRingStatistics(); //write out ring statistics
    void writeAboavWeaire(); //write out aboav weaire parameters
    void writeSpatialPartialRdfs(); //write out partial rdfs
    void writeTopoPartialRdfs(); //write out partial rdfs
    void writeRingAreas(); //write out mean areas ring sizes
    void writeClusters(); //write out size of clusters of each ring size
    void writeAssortativeMixing(); //write out assortative pearson correlation
    void writeAtomicGeometryOptimisation(); //write energy of atomic system
    void writeAtomicGeometrySummary(); //write bond and angle averages and stdev

public:
    void setIO(string in, string out); //set input/output properties
    void setProperties(bool per, bool readIn, string latType, vector<int> latDim, bool ovrPbc, vector<double> ovrLatPbc, vector<int> ringLim, double alpha, vector<double> p); //set initial and target network properties
    void setPotential(double sep, double k, bool local, int localMaxIt, double localCC,
                      bool global, int globalMaxIt, double globalCC, double lineInc); //set parameters for potential model
    void setMonteCarlo(int seed, double t, int moves, double conv, double asf); //set monte carlo limits
    void setAnalysis(bool convert, bool perVis, bool rdf, double rdfBw, double rdfExt, bool tRdf, double tRdfExt, bool rArea, bool clst, bool aMix); //set analysis tools
    void setAtomicPotential(bool opt, double kA, double kAlpha, double kBeta, int maxIt, double cc, double lineInc); //set parameters for atomic geometry optimisation
    bool getConsistency(); //get whether network constructed is consistent
    bool getTargetStatus(); //get whether target is met
    bool getIntersectionStatus(); //get whether has intersecting edges in dual
    void construct(ofstream &logfile); //aim to make network with supplied properties
    void analyse(ofstream &logfile); //analyse network properties
    void write(); //write out control
};


#endif //DUAL_SWITCH_NETWORK_H
