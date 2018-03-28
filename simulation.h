//read in parameters for simulation and control building and analysis of networks with different seeds

#ifndef DUAL_SWITCH_SIMULATION_H
#define DUAL_SWITCH_SIMULATION_H

#include <chrono>
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
    bool localGeomOpt, globalGeomOpt; //local geometry optimisation every step, global at end
    int localGeomOptMaxIt, globalGeomOptMaxIt; //iterations for local/global optimisation
    double atomDistance, harmonicK, localGeomOptConv, globalGeomOptConv, lineSearchStep; //between atoms in network, force constant, local/global force convergence criteria

    //Analysis
    bool convertDual, periodicVis, spatialRdf, topoRdf, assortative; //convert to atomic network, periodic visualisation, spatial partial rdfs, topological rdfs
    double spatialRdfBinWidth, spatialRdfExtent, topoRdfExtent; //histogram bin width and distance

    //Atomic Potential Model and Geometry Optimisation
    bool atomicGeomOpt; //geometry optimise atomic network
    double keatingA, keatingAlpha, keatingBeta; //keating potential model parameters, bond length, force constant and angle force constant
    int atomGeomOptMaxIt; //maximum iterations for geometry optimisation
    double atomGeomOptConv, atomLineSearchStep; //convergence condition and line search step size for atomic optimisation

    //Logfile
    string logfileName; //name of log file

    //Functions
    void readInputFile(); //read in simulation parameters
    void initialiseNetwork(Network &network, int seed); //set up network with input properties and given random seed


public:
    void run(); //called by main to initialise simulation
};



#endif //DUAL_SWITCH_SIMULATION_H
