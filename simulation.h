//read in parameters for simulation and control building and analysis of networks with different seeds

#ifndef DUAL_SWITCH_SIMULATION_H
#define DUAL_SWITCH_SIMULATION_H

#include "customIO.h"
#include "network.h"


class Simulation{
private:

    //IO
    string inPrefix, outPrefix; //for IO filenames

    //Network Properties
    bool periodic, load; //(a)periodic network, using existing network
    double alpha; //target aboav-weaire alpha parameter
    vector<int> latticeDimensions, ringSizeLimits; //initial lattice dimensions (nunmber of rings in x/y), min/max ring sizes
    vector<double> pVector; //target ring statistics

    //Simulation Parameters
    vector<int> randomSeedLimits; //random seeds
    int maxMoves; //for monte carlo
    double temperature, propConvergence, alphaEnergyScaling; //for monte carlo, convergence for network properties, factor for functional

    //Potential Model and Geometry Optimisation
    bool localGeomOpt, globalGeomOpt, resolveOverlaps; //local geometry optimisation every step, global at end, resolve edge overlaps in dual
    int localGeomOptMaxIt, globalGeomOptMaxIt; //iterations for local/global optimisation
    double atomDistance, harmonicK, localGeomOptConv, globalGeomOptConv, lineSearchStep; //between atoms in network, force constant, local/global force convergence criteria

    //Analysis
    bool periodicVis, spatialRdf, topoRdf; //periodic visualisation, spatial partial rdfs, topological rdfs
    int maxTopoShells; //shells for topo rdf
    double spatialRdfBinWidth; //histogram bin width

    //Logfile
    string logfileName; //name of log file

    //Calculation status
    vector<bool> networkConsistent, networkTargetsReached;

    void readInputFile(); //read in simulation parameters
    void initialiseNetwork(Network &network, int seed); //set up network with input properties and given random seed


public:
    void run(); //called by main to initialise simulation
};



#endif //DUAL_SWITCH_SIMULATION_H
