//network of single random seed

#ifndef DUAL_SWITCH_NETWORK_H
#define DUAL_SWITCH_NETWORK_H

#include <map>
#include "customIO.h"
#include "node.h"
#include "vectorManip.h"
#include "minimise.h"

using namespace std;

class Network {
private:

    //###### Input variables ######

    //IO variables
    string inPrefix, outPrefix; //for input and output files

    //Network properties
    bool periodic, load; //periodic network, load existing network
    double targetAlpha; //target for aboav-weaire alpha parameter
    vector<int> initialLatticeDimensions, ringSizeLimits; //starting hex lattice size, max/min ring size values
    vector<double> targetPVector; //target ring statistics

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
    bool periodicVisualisation, spatialRdf, topoRdf; //make periodic image, spatial partial rdfs, topological rdfs
    double spatialRdfBinwidth, spatialRdfExtent, topoRdfExtent; //size of bin rdf, max distance for rdf

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
    void initialiseMonteCarlo(); //set up variables and initialise random number generators

    //Random numbers
    int pickRandomNode(), pickRandomConnection(int nCnxs); //select node from network, select connection
    double metropolisRandomNum(); //random number for mc metropolis condition

    //Monte carlo
    void monteCarlo(); //main monte carlo process
    vector<int> pickRandomTrianglePairPeriodic(); //nodes to dual switch
    vector<int> pickRandomTrianglePairAperiodic(); //nodes to dual switch
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
    vector<double> ringStatistics; //ring statistics of network
    vector< vector<double> > piMatrix; //normalised p matrix
    vector<double> spatialRdfDensities; //density of nodes
    vector<Rdf> spatialPartialRdfs; //rdfs for each rings size pair
    vector< vector<int> > topoRdfShellSizes; //number of rings in each shell for each ring size
    vector<Rdf> topoPartialRdfs; //rdfs for each ring size pair

    //###### Analysis functions ######
    void analyseRingStatistics(); //calculate ring statistics and normalise pi matrix
    void analyseAboavWeaire(); //calculate aboav-weaire parameters
    void analysePartialSpatialRdfs(); //calculate partial rdfs for each ring size
    void analysePartialTopologicalRdfs(); //calculate partial topological rdfs for each ring size

    //###### Write functions ######
    void writeDual(); //write out dual network
    void writePeriodicNetwork(); //calculate and write out network for periodic visualisation
    void writeRingStatistics(); //write out ring statistics
    void writeAboavWeaire(); //write out aboav weaire parameters
    void writeSpatialPartialRdfs(); //write out partial rdfs
    void writeTopoPartialRdfs(); //write out partial rdfs

public:
    void setIO(string in, string out); //set input/output properties
    void setProperties(bool per, bool readIn, vector<int> latDim, vector<int> ringLim, double alpha, vector<double> p); //set initial and target network properties
    void setPotential(double sep, double k, bool local, int localMaxIt, double localCC,
                      bool global, int globalMaxIt, double globalCC, double lineInc); //set parameters for potential model
    void setMonteCarlo(int seed, double t, int moves, double conv, double asf); //set monte carlo limits
    void setAnalysis(bool perVis, bool rdf, double rdfBw, double rdfExt, bool tRdf, double tRdfExt); //set analysis tools
    bool getConsistency(); //get whether network constructed is consistent
    bool getTargetStatus(); //get whether target is met
    bool getIntersectionStatus(); //get whether has intersecting edges in dual
    void construct(ofstream &logfile); //aim to make network with supplied properties
    void analyse(ofstream &logfile); //analyse network properties
    void write(); //write out control
};


#endif //DUAL_SWITCH_NETWORK_H
