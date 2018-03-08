#include "network.h"

//###### SETTERS ######

void Network::setIO(string in, string out) {
    //set prefix for read in and write out to files
    inPrefix=in+"_";
    outPrefix=out+"_";
    return;
}

void Network::setProperties(bool per, bool readIn, vector<int> latDim, vector<int> ringLim, double alpha, vector<double> p) {
    //set initial lattice properties and target lattice properties
    periodic=per;
    load=readIn;
    initialLatticeDimensions=latDim;
    ringSizeLimits=ringLim;
    targetAlpha=alpha;
    targetPVector=p;
    return;
}

void Network::setMonteCarlo(int seed, double t, int moves, double conv, double asf) {
    //set up parameters for monte carlo dual switching
    mcSeed=seed;
    mcTemperature=t;
    mcMaxMoves=moves;
    mcConvergence=conv;
    mcAlphaScaleFactor=asf;
    return;
}

//###### Construction Main ######

void Network::construct(ofstream &logfile) {
    //attempt to build network with the specified properties

    writeFileLine(logfile,"constructing network seed "+to_string(mcSeed));
    initialiseNetworkProperties();


    return;
}

//###### Construction Initialisation ######

void Network::initialiseNetworkProperties() {
    //set up key variables for calculation

    minRingSize=ringSizeLimits[0];
    maxRingSize=ringSizeLimits[1];
    nRingSizes=maxRingSize-minRingSize+1;
    index6=6-minRingSize; //index of 6mr

    Node node();

    return;
}