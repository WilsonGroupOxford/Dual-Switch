//network of single random seed

#ifndef DUAL_SWITCH_NETWORK_H
#define DUAL_SWITCH_NETWORK_H

#include <map>
#include "customIO.h"
#include "node.h"
#include "vectorManip.h"

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
    double atomicSeparation; //ideal distance between atoms

    //Monte carlo dual switching properties
    int mcSeed, mcMaxMoves; //random seed, maximum proposed moves
    double mcTemperature, mcConvergence, mcAlphaScaleFactor; //temperature, properties convergence criteria, scale factor for energy functional

    //Analysis variables
    bool periodicVisualisation; //make periodic image

    //###### Construction variables ######

    //Further network properties
    bool consistent; //network and p vector/matrix consistent
    string name; //name of network for writing to files
    int minRingSize, maxRingSize, nRingSizes, index6, nNodes; //ring size limits, index of 6mr, number of nodes in network
    double periodicBoxX, periodicBoxY, rPeriodicBoxX, rPeriodicBoxY; //dimensions of periodic box and reciprocals
    vector<Node> nodes; //nodes in network
    vector<int> pVector; //number of nodes of each size
    vector< vector<int> > pMatrix; //number of nodes adjacent to nodes of each size
    vector<double> aboavWeaireParams; //alpha, mu, rsq

    //Further potential model variables
    vector< vector<double> > harmonicR0Matrix; //harmonic minimum for nodes of different sizes

    //Further monte carlo properties
    bool mcTargetReached; //whether targets met within convergence
    int mcProposedMoves; //number of moves proposed during monte carlo
    double mcEnergy, targetMu; //energy of system relative to target properties, mu
    double rMcTemperature, rTargetMu; //reciprocals of key values
    vector<double> rTargetPVector; //reciprocals of key values
    mt19937 randomGenerator; //generator for distributions
    uniform_int_distribution<int> nodeDistribution; //mersenne twister uniform distribution between 0->nNodes
    uniform_real_distribution<double> zeroOneDistribution, connectionPickDistribution; //mersenne twister uniform distribution between 0->1



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
    bool acceptDualSwitch(vector<int> &switchTriangles, vector<int> &trialPVec, vector< vector<int> > &trialPMat, double &trialMcEnergy, vector<double> &trialAwParams); //enact dual switch and update trial->current variables

    //Checking
    void checkFidelity(); //check to ensure consistency

    //###### Write functions ######
    void writeDual(); //write out dual network
    void writePeriodicNetwork(); //calculate and write out network for periodic visualisation

public:
    void setIO(string in, string out); //set input/output properties
    void setProperties(bool per, bool readIn, vector<int> latDim, vector<int> ringLim, double alpha, vector<double> p); //set initial and target network properties
    void setPotential(double sep); //set parameters for potential model
    void setMonteCarlo(int seed, double t, int moves, double conv, double asf); //set monte carlo limits
    void setAnalysis(bool perVis); //set analysis tools
    void construct(ofstream &logfile); //aim to make network with supplied properties
    void write(); //write out control
};


#endif //DUAL_SWITCH_NETWORK_H
