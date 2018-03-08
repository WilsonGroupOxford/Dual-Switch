//network of single random seed

#ifndef DUAL_SWITCH_NETWORK_H
#define DUAL_SWITCH_NETWORK_H

#include "customIO.h"
#include "node.h"

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

    //Monte carlo dual switching properties
    int mcSeed, mcMaxMoves; //random seed, maximum proposed moves
    double mcTemperature, mcConvergence, mcAlphaScaleFactor; //temperature, properties convergence criteria, scale factor for energy functional

    //###### Construction variables ######

    //Further network properties
    int minRingSize, maxRingSize, nRingSizes, index6; //ring size limits, index of 6mr


    //Further monte carlo properties




    //###### Construction functions ######

    //Initialisation
    void initialiseNetworkProperties(); //additional network properites


    void initialiseMonteCarlo(); //set up variables and initialise random number generators


public:
    void setIO(string in, string out); //set input/output properties
    void setProperties(bool per, bool readIn, vector<int> latDim, vector<int> ringLim, double alpha, vector<double> p); //set initial and target network properties
    void setMonteCarlo(int seed, double t, int moves, double conv, double asf); //set monte carlo limits
    void construct(ofstream &logfile); //aim to make network with supplied properties
};


#endif //DUAL_SWITCH_NETWORK_H
