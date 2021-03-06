#include "network.h"

//###### SETTERS ######

void Network::setIO(string in, string out) {
    //set prefix for read in and write out to files
    inPrefix=in+"_";
    outPrefix=out+"_";
    return;
}

void Network::setProperties(bool per, bool readIn, string latType, vector<int> latDim, bool ovrPbc, vector<double> ovrLatPbc, vector<int> ringLim, double alpha, vector<double> p) {
    //set initial lattice properties and target lattice properties
    periodic=per;
    load=readIn;
    initialLatticeType=latType;
    initialLatticeDimensions=latDim;
    overridePbc=ovrPbc;
    overridelatticePbc=ovrLatPbc;
    ringSizeLimits=ringLim;
    targetAlpha=alpha;
    targetPVector=p;
    return;
}

void Network::setPotential(double sep, double k, bool local, int localMaxIt, double localCC, bool global, int globalMaxIt,
                           double globalCC, double lineInc) {
    //set potential model parameters
    atomicSeparation=sep;
    harmonicK=k;
    localGeomOpt=local;
    localGeomOptMaxIt=localMaxIt;
    localGeomOptCC=localCC;
    globalGeomOpt=global;
    globalGeomOptMaxIt=globalMaxIt;
    globalGeomOptCC=globalCC;
    geomOptLineSearchInc=lineInc;
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

void Network::setAnalysis(bool convert, bool perVis, bool rdf, double rdfBw, double rdfExt, bool tRdf, double tRdfExt, bool rArea, bool clst, bool aMix){
    //specify analysis to perform
    convertToAtomic=convert;
    if(periodic) periodicVisualisation=perVis;
    else periodicVisualisation=false;
    spatialRdf=rdf;
    spatialRdfBinwidth=rdfBw;
    spatialRdfExtent=rdfExt;
    topoRdf=tRdf;
    topoRdfExtent=tRdfExt;
    areaLaw=rArea;
    clustering=clst;
    assortativeMix=aMix;
    return;
}

void Network::setAtomicPotential(bool opt, double kA, double kAlpha, double kBeta, int maxIt, double cc,
                                 double lineInc) {
    //set keating potential and geometry optimisation parameters
    atomicGeomOpt=opt;
    keatingA=kA;
    keatingAlpha=kAlpha;
    keatingBeta=kBeta;
    atomicGeomOptMaxIt=maxIt;
    atomicGeomOptCC=cc;
    atomicLineSearchInc=lineInc;
    return;
}
//###### GETTERS ######

bool Network::getConsistency() {
    //whether generated network is self consistent
    return consistent;
}

bool Network::getTargetStatus() {
    //whether network achieved targets
    return mcTargetReached;
}

bool Network::getIntersectionStatus() {
    //whether there are no intersections in the dual
    return noIntersections;
}

//###### Construction Main ######

void Network::construct(ofstream &logfile) {
    //attempt to build network with the specified properties

    name="network "+to_string(mcSeed);
    writeFileLine(logfile,"constructing "+name);

    initialiseNetworkProperties();
    initialisePotentialModel();
    if(load){
        loadPeriodicLattice();
        if(consistent) writeFileLine(logfile,name+" loaded successfully from "+inPrefix);
        else writeFileLine(logfile,name+" loaded unsuccesfully from "+inPrefix);
    }
    else if(periodic) initialisePeriodicLattice();
    else initialiseAperiodicLattice();

    initialiseMonteCarlo();

    if(initialLatticeType!="hexagonal"){
//        if(!periodic) initialiseAperiodicCrystalLattice();
        if(periodic) initialisePeriodicCrystalLattice(initialLatticeType);
    }

    monteCarlo();
    checkFidelity();

    writeFileLine(logfile,name+" constructed with "+to_string(mcProposedMoves)+" proposed mc moves with "+to_string(mcAcceptedMoves)+" accepted");
    if(localGeomOpt) writeFileLine(logfile,name+" had "+to_string(geomOptRejectCount)+" moves rejected due to unresolved geometry optimisation");
    if(consistent) writeFileLine(logfile,name+" checked for consistency and passed");
    else writeFileLine(logfile,name+" failed consistency test");
    if(mcTargetReached) writeFileLine(logfile, name+" targets met in "+to_string(mcProposedMoves)+" monte carlo moves");
    else writeFileLine(logfile, name+" targets not met");

    if(consistent){
        findNodeRings();
        checkGeometry();
        if(noIntersections) writeFileLine(logfile,name+" dual has no intersections");
        else writeFileLine(logfile,name+" dual has intersections");

        if(noIntersections && globalGeomOpt){
            if(periodic) globalMinimisationPeriodic();
            else globalMinimisationAperiodic();
            writeFileLine(logfile,name+" globally geometry optimised with "+to_string(geomOptIterations)+" iterations of steepest descent");
        }
    }

    return;
}

//###### Construction Initialisation ######

void Network::initialiseNetworkProperties() {
    //set up key variables for calculation

    minRingSize=ringSizeLimits[0];
    maxRingSize=ringSizeLimits[1];
    nRingSizes=maxRingSize-minRingSize+1;
    index6=6-minRingSize; //index of 6mr
    nodes.clear();
    nodeRings.clear();
    aboavWeaireParams.resize(3,0.0);

    return;
}

void Network::initialisePotentialModel() {
    //set up parameters for potential model

    //ideal separations between nodes of different sizes
    harmonicR0Matrix.resize(nRingSizes, (vector<double> (nRingSizes)));
    for (int i=0; i<nRingSizes; ++i){
        int n=i+minRingSize;
        for (int j=0; j<nRingSizes; ++j){
            int m=j+minRingSize;
            harmonicR0Matrix[i][j]=0.5*atomicSeparation*(1.0/tan(pi/n)+1.0/tan(pi/m));
        }
    }

    return;
}

void Network::initialisePeriodicLattice() {
    //make periodic hexagonal lattice of input dimensions

    nNodes=initialLatticeDimensions[0]*initialLatticeDimensions[1]; //set number of nodes as x*y

    //set up nodes of size 6, with ring size limits and not on edge
    for(int i=0; i<nNodes; ++i) nodes.push_back(Node(6,minRingSize,maxRingSize));

    //make coordinates of hexagonal lattice of alternating staggered layers (y layers each containing x nodes)
    int xNodes=initialLatticeDimensions[0], yNodes=initialLatticeDimensions[1]; //number of nodes in layer and number of layers
    bool stagger=false; //whether to stagger - will alternate
    double r66=harmonicR0Matrix[6-minRingSize][6-minRingSize]; //ideal 6-6 distance
    double xStagger=0.5*r66, yOffset=0.5*sqrt(3)*r66; //shift in x (intralayer), shift in y (interlayer)
    int nodeIndex; //index of node coordinates being set
    double xCrd; //successive x coordinates
    for(int layer=0; layer<yNodes; ++layer){
        if(stagger) xCrd=xStagger;
        else xCrd=0.0;
        for(int n=0; n<xNodes; ++n){
            nodeIndex=layer*xNodes+n;
            nodes[nodeIndex].coordinate.y=layer*yOffset;
            nodes[nodeIndex].coordinate.x=xCrd;
            xCrd=xCrd+r66;
        }
        stagger=!stagger;
    }

    //make connections between nodes - aperiodic first
    nodeIndex=0;
    for(int j=0; j<yNodes; ++j){//intralayer
        for(int i=0; i<xNodes-1; ++i){
            nodes[nodeIndex].addConnection(nodeIndex+1);
            nodes[nodeIndex+1].addConnection(nodeIndex);
            nodeIndex=++nodeIndex;
        }
        nodeIndex=++nodeIndex; //skip to next layer
    }
    nodeIndex=0;
    stagger=false;
    for(int j=0; j<yNodes-1;++j){//interlayer
        if(!stagger){
            nodes[nodeIndex].addConnection(nodeIndex+xNodes);
            nodes[nodeIndex+xNodes].addConnection(nodeIndex);
            nodeIndex=++nodeIndex;
            for(int i=1; i<xNodes; ++i){
                nodes[nodeIndex].addConnection(nodeIndex+xNodes-1);
                nodes[nodeIndex].addConnection(nodeIndex+xNodes);
                nodes[nodeIndex+xNodes-1].addConnection(nodeIndex);
                nodes[nodeIndex+xNodes].addConnection(nodeIndex);
                nodeIndex=++nodeIndex;
            }
        }
        else{
            for(int i=0; i<xNodes-1; ++i){//bottom layer <-> middle layer
                nodes[nodeIndex].addConnection(nodeIndex+xNodes);
                nodes[nodeIndex].addConnection(nodeIndex+xNodes+1);
                nodes[nodeIndex+xNodes].addConnection(nodeIndex);
                nodes[nodeIndex+xNodes+1].addConnection(nodeIndex);
                nodeIndex=++nodeIndex;
            }
            nodes[nodeIndex].addConnection(nodeIndex+xNodes);
            nodes[nodeIndex+xNodes].addConnection(nodeIndex);
            nodeIndex=++nodeIndex;
        }
        stagger=!stagger;
    }
    //add periodic connections
    nodeIndex=0;
    for(int j=0; j<yNodes; ++j){//intralayer
        nodes[nodeIndex].addConnection(nodeIndex+xNodes-1);
        nodes[nodeIndex+xNodes-1].addConnection(nodeIndex);
        nodeIndex=nodeIndex+xNodes; //skip to next layer
    }
    stagger=false;
    nodeIndex=0;
    for(int j=0; j<yNodes-1; ++j){//interlayer connections excluding bottom<->top
        if(!stagger){
            nodes[nodeIndex].addConnection(nodeIndex+2*xNodes-1);
            nodes[nodeIndex+2*xNodes-1].addConnection(nodeIndex);
        }
        else{
            nodes[nodeIndex+xNodes-1].addConnection(nodeIndex+xNodes);
            nodes[nodeIndex+xNodes].addConnection(nodeIndex+xNodes-1);
        }
        nodeIndex=nodeIndex+xNodes;
        stagger=!stagger;
    }
    for(int i=0; i<xNodes; ++i){//interlayer bottom<->top
        nodes[i].addConnection(nodeIndex+i);
        nodes[nodeIndex+i].addConnection(i);
    }
    for(int i=1; i<xNodes; ++i){//interlayer bottom<->top
        nodes[i].addConnection(nodeIndex+i-1);
        nodes[nodeIndex+i-1].addConnection(i);
    }
    nodes[0].addConnection(xNodes*yNodes-1); //bottom left to top right
    nodes[xNodes*yNodes-1].addConnection(0);

    //set up periodic dimensions
    periodicBoxX=r66*xNodes;
    periodicBoxY=yOffset*yNodes;
    rPeriodicBoxX=1.0/periodicBoxX;
    rPeriodicBoxY=1.0/periodicBoxY;

    //set up p vector and matrix
    pVector.resize(nRingSizes,0);
    pMatrix.resize(nRingSizes, (vector<int> (nRingSizes,0)));
    pVector[index6]=initialLatticeDimensions[0]*initialLatticeDimensions[1]; //set all rings as hexagons
    pMatrix[index6][index6]=nNodes*6;

    return;
}

void Network::initialiseAperiodicLattice() {
    //make aperiodic lattice of input dimensions

    //need extra nodes on edges and extra layer on top and bottom to ensure have inputed number of non-edge nodes
    nNodes=(initialLatticeDimensions[0]+2)*(initialLatticeDimensions[1]+2)-2; //set number of nodes as (x+2)*(y+2)-2, -2 as top and bottom layer need one few node

    //set up nodes of size 6, with ring size limits
    int xNodes=initialLatticeDimensions[0]+2, yNodes=initialLatticeDimensions[1]+2; //number of nodes in layer and number of layers
    for(int i=0; i<xNodes-1; ++i) nodes.push_back(Node(6,minRingSize,maxRingSize,true)); //bottom edge
    for(int i=0; i<yNodes-2; ++i){
        nodes.push_back(Node(6,minRingSize,maxRingSize,true)); //left edge
        for(int j=1; j<xNodes-1; ++j) nodes.push_back(Node(6,minRingSize,maxRingSize)); //bulk
        nodes.push_back(Node(6,minRingSize,maxRingSize,true)); //right edge
    }
    for(int i=0; i<xNodes-1; ++i) nodes.push_back(Node(6,minRingSize,maxRingSize,true)); //top edge

    //make coordinates of hexagonal lattice of alternating staggered layers
    bool stagger=false; //whether to stagger - will alternate
    double r66=harmonicR0Matrix[6-minRingSize][6-minRingSize]; //ideal 6-6 distance
    double xStagger=0.5*r66, yOffset=0.5*sqrt(3)*r66; //shift in x (intralayer), shift in y (interlayer)
    int nodeIndex=0; //index of node coordinates being set
    double xCrd; //successive x coordinates
    for(int j=0; j<xNodes-1; ++j){//bottom layer
        nodes[nodeIndex].coordinate.x=xStagger+r66*j;
        nodeIndex=++nodeIndex;
    }
    for(int i=1; i<yNodes-1; ++i){
        if(stagger) xCrd=xStagger;
        else xCrd=0.0;
        for(int j=0; j<xNodes; ++j){
            nodes[nodeIndex].coordinate.x=xCrd;
            nodes[nodeIndex].coordinate.y=i*yOffset;
            xCrd=xCrd+r66;
            nodeIndex=++nodeIndex;
        }
        stagger=!stagger;
    }
    for(int j=0; j<xNodes-1; ++j) {//top layer
        nodes[nodeIndex].coordinate.x=xStagger*2+r66*j;
        nodes[nodeIndex].coordinate.y=(yNodes-1)*yOffset;
        nodeIndex=++nodeIndex;
    }

    //make connections between nodes
    //intralayer
    nodeIndex=0;
    for(int i=0; i<xNodes-2; ++i){//bottom layer
        nodes[nodeIndex].addConnection(nodeIndex+1);
        nodes[nodeIndex+1].addConnection(nodeIndex);
        nodeIndex=++nodeIndex;
    }
    nodeIndex=++nodeIndex;
    for(int j=1; j<yNodes-1; ++j){//middle layers
        for(int i=0; i<xNodes-1; ++i){
            nodes[nodeIndex].addConnection(nodeIndex+1);
            nodes[nodeIndex+1].addConnection(nodeIndex);
            nodeIndex=++nodeIndex;
        }
        nodeIndex=++nodeIndex; //skip to next layer
    }
    for(int i=0; i<xNodes-2; ++i){//top layer
        nodes[nodeIndex].addConnection(nodeIndex+1);
        nodes[nodeIndex+1].addConnection(nodeIndex);
        nodeIndex=++nodeIndex;
    }
    //interlayer
    nodeIndex=0;
    for(int i=0; i<xNodes-1; ++i){//bottom layer <-> middle layer
        nodes[nodeIndex].addConnection(nodeIndex+xNodes-1);
        nodes[nodeIndex].addConnection(nodeIndex+xNodes);
        nodes[nodeIndex+xNodes-1].addConnection(nodeIndex);
        nodes[nodeIndex+xNodes].addConnection(nodeIndex);
        nodeIndex=++nodeIndex;
    }
    stagger=false;
    for(int j=1; j<yNodes-2;++j){//middle layer <-> middle layer
        if(!stagger){
            nodes[nodeIndex].addConnection(nodeIndex+xNodes);
            nodes[nodeIndex+xNodes].addConnection(nodeIndex);
            nodeIndex=++nodeIndex;
            for(int i=1; i<xNodes; ++i){
                nodes[nodeIndex].addConnection(nodeIndex+xNodes-1);
                nodes[nodeIndex].addConnection(nodeIndex+xNodes);
                nodes[nodeIndex+xNodes-1].addConnection(nodeIndex);
                nodes[nodeIndex+xNodes].addConnection(nodeIndex);
                nodeIndex=++nodeIndex;
            }
        }
        else{
            for(int i=0; i<xNodes-1; ++i){//bottom layer <-> middle layer
                nodes[nodeIndex].addConnection(nodeIndex+xNodes);
                nodes[nodeIndex].addConnection(nodeIndex+xNodes+1);
                nodes[nodeIndex+xNodes].addConnection(nodeIndex);
                nodes[nodeIndex+xNodes+1].addConnection(nodeIndex);
                nodeIndex=++nodeIndex;
            }
            nodes[nodeIndex].addConnection(nodeIndex+xNodes);
            nodes[nodeIndex+xNodes].addConnection(nodeIndex);
            nodeIndex=++nodeIndex;
        }
        stagger=!stagger;
    }
    //middle layer <-> top layer
    nodes[nodeIndex].addConnection(nodeIndex+xNodes);
    nodes[nodeIndex+xNodes].addConnection(nodeIndex);
    nodeIndex=++nodeIndex;
    for(int i=1; i<xNodes-1; ++i){
        nodes[nodeIndex].addConnection(nodeIndex+xNodes-1);
        nodes[nodeIndex].addConnection(nodeIndex+xNodes);
        nodes[nodeIndex+xNodes-1].addConnection(nodeIndex);
        nodes[nodeIndex+xNodes].addConnection(nodeIndex);
        nodeIndex=++nodeIndex;
    }
    nodes[nodeIndex].addConnection(nodeIndex+xNodes-1);
    nodes[nodeIndex+xNodes-1].addConnection(nodeIndex);

    //correct connectivities for edge nodes
    for(int i=0; i<nNodes; ++i){
        if(nodes[i].edge){
            nodes[i].size=nodes[i].connections.size();
            nodes[i].sizeIndex=nodes[i].size-minRingSize;
        }
    }

    //set up p vector and matrix
    int nodeRef;
    pVector.resize(nRingSizes,0);
    pMatrix.resize(nRingSizes, (vector<int> (nRingSizes,0)));
    pVector[index6]=initialLatticeDimensions[0]*initialLatticeDimensions[1]; //set all rings as hexagons
    for(int i=0; i<nNodes; ++i){
        if(!nodes[i].edge){
            for(int j=0; j<nodes[i].size; ++j){
                nodeRef=nodes[i].connections[j];
                if(!nodes[nodeRef].edge) pMatrix[nodes[i].sizeIndex][nodes[nodeRef].sizeIndex]=++pMatrix[nodes[i].sizeIndex][nodes[nodeRef].sizeIndex];
            }
        }
    }

    return;
}

void Network::initialisePeriodicCrystalLattice(string crystal){
    //set up a starting lattice which is not hexagonal by performing well-defined switch moves

    string latticeCode;
    if(crystal=="haeckelite") latticeCode="c0";

    if(latticeCode.substr(0,1)=="h"){
        int skip=stoi(latticeCode.substr(1,2));
        for(int i=0; i<initialLatticeDimensions[1]; i+=2) hexPeriodicCrystal(i,0,1,skip,false);
    }
    else if(latticeCode.substr(0,1)=="c"){
        int skip=stoi(latticeCode.substr(1,2));
        for(int i=0, j=0; i<initialLatticeDimensions[1]; i+=2, ++j) cubicPeriodicCrystal(i,j%2,1,skip);
    }
    else if(latticeCode=="xx") {
        for (int i = 0; i < initialLatticeDimensions[1]; i+=3) cubicPeriodicCrystal(i,0,1,0);
        for(int i=1; i<initialLatticeDimensions[1]; i+=6) hexPeriodicCrystal(i,0,1,0,false);
        for(int i=4; i<initialLatticeDimensions[1]; i+=6) hexPeriodicCrystal(i,1,1,0,true);
    }
//    else if(latticeCode=="xxx") {
//        for (int i = 0; i < initialLatticeDimensions[1]; i+=3) cubicPeriodicCrystal(i,0,2,2);
//        hexPeriodicCrystal(1,0,2,2,false);
//        hexPeriodicCrystal(4,1,2,2,true);
//        hexPeriodicCrystal(3,1,2,2,true);
////        hexPeriodicCrystal(9,0,1,0,true);
//        hexPeriodicCrystal(6,0,1,0,false);
//        hexPeriodicCrystal(9,0,1,0,true);
//        cubicPeriodicCrystal(2,1,1,10);
//        hexPeriodicCrystal(2,0,1,10,false);
//        for(int i=1; i<initialLatticeDimensions[1]; i+=6) hexPeriodicCrystal(i,0,1,0,false);
//        for(int i=4; i<initialLatticeDimensions[1]; i+=6) hexPeriodicCrystal(i,1,1,0,true);
//    }

}

void Network::hexPeriodicCrystal(int layer, int offset, int blockLength, int blockSkip, bool reverse) {
    //apply dual switch moves to crystal in hexagonal symmetry

    //block variables
    int repeatUnitLength=blockLength*2+blockSkip; //length of repeating unit (defect = 2 rings in length)
    int xNodes=initialLatticeDimensions[0], yNodes=initialLatticeDimensions[1]; //number of nodes in layer and number of layer

    //make list of nodes which need switching
    vector<bool> nodeSwitch(xNodes, false);
    int pos;
    for (int i=0, j=layer*xNodes; i<xNodes; ++i, ++j){
        pos=i%repeatUnitLength;
        if(pos<blockLength*2 && pos%2==offset){
            nodeSwitch[i]=true;
        }
    }

    //perform switches
    vector<int> switchTriangles(4);
    if(!reverse){
        if (layer % 2 == 0) {
            for (int i = 0, j = xNodes * layer; i < xNodes; ++i, ++j) {
                if (nodeSwitch[i]) {
                    switchTriangles[0] = j;
                    switchTriangles[1] = j + xNodes;
                    switchTriangles[2] = j + xNodes - 1;
                    switchTriangles[3] = j + 1;
                    if (i % xNodes == 0) switchTriangles[2] += xNodes;
                    else if ((i + 1) % xNodes == 0) switchTriangles[3] -= xNodes;
                    consoleVector(switchTriangles);
                    definedMove(switchTriangles);
                }
            }
        } else {
            for (int i = 0, j = xNodes * layer; i < xNodes; ++i, ++j) {
                if (nodeSwitch[i]) {
                    switchTriangles[0] = j;
                    switchTriangles[1] = j + xNodes + 1;
                    switchTriangles[2] = j + xNodes;
                    switchTriangles[3] = j + 1;
                    if ((i + 1) % xNodes == 0) {
                        switchTriangles[1] -= xNodes;
                        switchTriangles[3] -= xNodes;
                    }
                    if (layer + 1 == yNodes) {
                        switchTriangles[1] -= nNodes;
                        switchTriangles[2] -= nNodes;
                    }
                    definedMove(switchTriangles);
                }
            }
        }
    }
    else{
        if (layer % 2 == 0) {
            for (int i = 0, j = xNodes * layer; i < xNodes; ++i, ++j) {
                if (nodeSwitch[i]) {
                    switchTriangles[0] = j;
                    switchTriangles[1] = j + xNodes -1;
                    switchTriangles[2] = j + xNodes;
                    switchTriangles[3] = j - 1;
                    if (i % xNodes == 0){
                        switchTriangles[1] += xNodes;
                        switchTriangles[3] += xNodes;
                    }
                    definedMove(switchTriangles);
                }
            }
        } else {
            for (int i = 0, j = xNodes * layer; i < xNodes; ++i, ++j) {
                if (nodeSwitch[i]) {
                    switchTriangles[0] = j;
                    switchTriangles[1] = j + xNodes;
                    switchTriangles[2] = j + xNodes + 1;
                    switchTriangles[3] = j - 1;
                    if (i % xNodes == 0) switchTriangles[3] += xNodes;
                    if ((i + 1) % xNodes == 0) switchTriangles[2] -= xNodes;
                    if (layer + 1 == yNodes) {
                        switchTriangles[1] -= nNodes;
                        switchTriangles[2] -= nNodes;
                    }
                    definedMove(switchTriangles);
                }
            }
        }
    }


}

void Network::cubicPeriodicCrystal(int layer, int offset, int blockLength, int blockSkip) {
    //apply dual switch moves to crystal in cubic symmetry

    //block variables
    int repeatUnitLength=blockLength*2+blockSkip; //length of repeating unit (defect = 2 rings in length)
    int xNodes=initialLatticeDimensions[0], yNodes=initialLatticeDimensions[1]; //number of nodes in layer and number of layer

    //make list of nodes which need switching
    vector<bool> nodeSwitch(xNodes, false);
    int pos;
    for (int i=0, j=layer*xNodes; i<xNodes; ++i, ++j){
        pos=i%repeatUnitLength;
        if(pos<blockLength*2 && pos%2==offset){
            nodeSwitch[i]=true;
        }
    }

    //perform switches
    vector<int> switchTriangles(4);
    if(layer%2==0){
        for(int i=0,j=xNodes*layer; i<xNodes; ++i, ++j) {
            if(nodeSwitch[i]) {
                switchTriangles[0] = j;
                switchTriangles[1] = j + 1;
                switchTriangles[2] = j + xNodes;
                switchTriangles[3] = j - xNodes;
                if ((i + 1) % xNodes == 0) switchTriangles[1] -= xNodes;
                if(layer==0) switchTriangles[3] += nNodes;
                definedMove(switchTriangles);
            }
        }
    }
    else{
        for(int i=0,j=xNodes*layer; i<xNodes; ++i, ++j) {
            if(nodeSwitch[i]) {
                switchTriangles[0] = j;
                switchTriangles[1] = j + 1;
                switchTriangles[2] = j + xNodes + 1;
                switchTriangles[3] = j - xNodes + 1;
                if ((i + 1) % xNodes == 0){
                    switchTriangles[1] -= xNodes;
                    switchTriangles[2] -= xNodes;
                    switchTriangles[3] -= xNodes;
                }
                if(layer+1==yNodes) switchTriangles[2] -= nNodes;
                definedMove(switchTriangles);
            }
        }
    }


}


void Network::initialiseAperiodicCrystalLattice(){
    //set up a lattice which is not hexagonal by performing well defined dual switch moves
    //TESTER FOR PERIODIC CASE

    //vars
    int xSwitchFreq=2;
    int ySwitchFreq=1;

    //scan crystal from bottom to top and perform switches at required frequencies
    int xNodes=initialLatticeDimensions[0]+2, yNodes=initialLatticeDimensions[1]+2; //number of nodes in layer and number of layers
    int nodeIndex=0;
    bool stagger=false;
    bool switchMove;
    vector<int> switchTriangles(4);
    //bottom edge
    for(int j=0; j<xNodes-1; ++j){
        if(j % xSwitchFreq==1){
            switchTriangles[0] = nodeIndex;
            switchTriangles[1] = nodeIndex + xNodes - 1;
            switchTriangles[2] = nodeIndex - 1;
            switchTriangles[3] = nodeIndex + xNodes;
            definedMove(switchTriangles);
        }
        ++nodeIndex;
    }
    //middle section
    for(int i=1; i<yNodes-2; ++i){
        if(i % ySwitchFreq==0){
            for(int j=0; j<xNodes; ++j){
                if(stagger && j%xSwitchFreq==0){
                    switchTriangles[0] = nodeIndex;
                    switchTriangles[1] = nodeIndex + xNodes + 1;
                    switchTriangles[2] = nodeIndex + 1;
                    switchTriangles[3] = nodeIndex + xNodes;
                    definedMove(switchTriangles);
                };
                if(!stagger && j%xSwitchFreq==1){
                    switchTriangles[0] = nodeIndex;
                    switchTriangles[1] = nodeIndex + xNodes - 1;
                    switchTriangles[2] = nodeIndex - 1;
                    switchTriangles[3] = nodeIndex + xNodes;
                    definedMove(switchTriangles);
                };
                ++nodeIndex;
            }
        }
        else nodeIndex+=xNodes;
        stagger=!stagger;
    }
    //top edge
    for(int j=0; j<xNodes-1; ++j){
        if(j % xSwitchFreq==1){
            switchTriangles[0] = nodeIndex;
            switchTriangles[1] = nodeIndex + xNodes - 1;
            switchTriangles[2] = nodeIndex - 1;
            switchTriangles[3] = nodeIndex + xNodes;
            definedMove(switchTriangles);
        }
        ++nodeIndex;
    }


//    for(int i=0; i<yNodes-1; ++i){//scan inter-layer
//        if(i % ySwitchFreq==0){
//            ++nodeIndex; //skip edge node
//            for(int j=0; j<xNodes-2; ++j){//scan intra-layer
//                if(stagger && j % xSwitchFreq==0) switchMove=true;
//                else if(!stagger && j % xSwitchFreq==1) switchMove=true;
//                else switchMove=false;
//                if(switchMove) {
//                    switchTriangles[0] = nodeIndex;
//                    switchTriangles[1] = nodeIndex + xNodes + 1;
//                    switchTriangles[2] = nodeIndex + 1;
//                    switchTriangles[3] = nodeIndex + xNodes;
//                    consoleVector(switchTriangles);
//                }
//                ++nodeIndex;
//            }
//            ++nodeIndex; //skip edge node
//        }
//        else nodeIndex+=xNodes;
//        stagger=!stagger;
//    }
}

void Network::initialiseMonteCarlo() {
    //set up number generators and variables for mc

    //energy infinite so first move accepted
    mcEnergy=numeric_limits<double>::infinity();

    //reciprocal values for faster calculation
    if(mcTemperature==0.0) rMcTemperature=std::numeric_limits<double>::infinity();
    else rMcTemperature=1.0/mcTemperature;
    rTargetPVector.clear();
    rTargetPVector.resize(nRingSizes);
    for(int i=0; i<nRingSizes; ++i) {
        if(targetPVector[i]>0) rTargetPVector[i]=1.0/targetPVector[i];
        else rTargetPVector[i]=1.0;
    }
    targetMu=0.0;
    for(int i=0; i<nRingSizes;++i) targetMu=targetMu+(minRingSize+i)*(minRingSize+i)*targetPVector[i];
    targetMu=targetMu-36.0;
    rTargetMu=1.0/targetMu;

    //random generators
    randomGenerator.seed(mcSeed);
    zeroOneDistribution=uniform_real_distribution<double>(0.0,1.0);
    nodeDistribution=uniform_int_distribution<int>(0,nNodes-1); //-1 as inclusive

    return;
}

void Network::loadPeriodicLattice() {
    //read in existing periodic lattice to start simulation

    //number of nodes from main input file
    nNodes=initialLatticeDimensions[0]*initialLatticeDimensions[1]; //set number of nodes as x*y

    //read in dual data
    string dualCnxFilename=inPrefix+"dual_connectivity.out";
    string dualCrdFilename=inPrefix+"dual_coordinates.out";
    ifstream dualCnxFile(dualCnxFilename,ios::in);
    ifstream dualCrdFile(dualCrdFilename,ios::in);
    vector< vector<int> > dualCnxs;
    vector< vector<double> > dualCrds;
    readFileAdaptiveMatrix(dualCnxFile,dualCnxs,nNodes);
    readFileMatrix(dualCrdFile,dualCrds,nNodes,2);
    dualCnxFile.close();
    dualCrdFile.close();

    //set up nodes of read in size, with ring size limits and not on edge
    for(int i=0; i<nNodes; ++i){
        nodes.push_back(Node(dualCnxs[i].size()-1,minRingSize,maxRingSize));
    }

    //set up coordinates
    for(int i=0; i<nNodes; ++i){
        nodes[i].coordinate.x=dualCrds[i][0];
        nodes[i].coordinate.y=dualCrds[i][1];
    }

    //set up connections
    for(int i=0; i<nNodes; ++i){
        for(int j=1; j<dualCnxs[i].size(); ++j){
            nodes[i].addConnection(dualCnxs[i][j]);
        }
    }

    //set up periodic dimensions
    string pbcFilename=inPrefix+"periodic_lattice_dim.out";
    ifstream pbcFile(pbcFilename,ios::in);
    double originalPbcX, originalPbcY;
    readFileValue(pbcFile,originalPbcX);
    readFileValue(pbcFile,originalPbcY);
    pbcFile.close();
    if(!overridePbc){
        periodicBoxX=originalPbcX;
        periodicBoxY=originalPbcY;
    }
    else{
        //scale dual coordinates to fit in new box
        double sfX=overridelatticePbc[0]/originalPbcX;
        double sfY=overridelatticePbc[1]/originalPbcY;
        Crd2d originalBoxCentre(periodicBoxX*0.5,periodicBoxY*0.5);
        for(int i=0; i<nNodes; ++i){
            nodes[i].coordinate.x=sfX*(nodes[i].coordinate.x-originalBoxCentre.x)+originalBoxCentre.x;
            nodes[i].coordinate.y=sfY*(nodes[i].coordinate.y-originalBoxCentre.y)+originalBoxCentre.y;
        }
        periodicBoxX=overridelatticePbc[0];
        periodicBoxY=overridelatticePbc[1];
    }
    rPeriodicBoxX=1.0/periodicBoxX;
    rPeriodicBoxY=1.0/periodicBoxY;

    //set up p vector and matrix
    pVector.resize(nRingSizes,0);
    pMatrix.resize(nRingSizes, (vector<int> (nRingSizes,0)));
    int sizeIndexI, sizeIndexJ;
    for(int i=0; i<nNodes; ++i){
        sizeIndexI=nodes[i].sizeIndex;
        pVector[sizeIndexI]=++pVector[sizeIndexI];
        for(int j=0; j<nodes[i].size; ++j){
            sizeIndexJ=nodes[nodes[i].connections[j]].sizeIndex;
            pMatrix[sizeIndexI][sizeIndexJ]=++pMatrix[sizeIndexI][sizeIndexJ];
        }
    }

    //calculate aboav-weaire fit
    aboavWeaireParams=calculateAboavWeaireFit(pVector,pMatrix);
    checkFidelity();

    return;
}

//###### Construction monte carlo ######

void Network::definedMove(vector<int> &switchTriangles) {
    //perform single defined dual-switch move

    //move variables
    double moveMcEnergy;
    vector<int> movePVector;
    vector< vector<int> > movePMatrix;
    vector<double> moveAwParameters;

    if(periodic){
        movePVector=pVector;
        movePMatrix=pMatrix;
        calculateTrialPPeriodic(switchTriangles, movePVector, movePMatrix);
        moveAwParameters=calculateAboavWeaireFit(movePVector,movePMatrix);
        moveMcEnergy=mcEnergyFunctional(moveAwParameters,movePVector);
        mcTargetReached=acceptDualSwitchPeriodic(switchTriangles,movePVector,movePMatrix,moveMcEnergy,moveAwParameters);
    }
    else{
        movePVector=pVector;
        movePMatrix=pMatrix;
        calculateTrialPAperiodic(switchTriangles, movePVector, movePMatrix);
        moveAwParameters=calculateAboavWeaireFit(movePVector,movePMatrix);
        moveMcEnergy=mcEnergyFunctional(moveAwParameters,movePVector);
        mcTargetReached=acceptDualSwitchAperiodic(switchTriangles,movePVector,movePMatrix,moveMcEnergy,moveAwParameters);
    }
}

void Network::monteCarlo() {
    //main dual switch monte carlo process

    //trial variables
    bool acceptTrialMove;
    double trialMcEnergy;
    vector<int> trialPVector;
    vector< vector<int> > trialPMatrix;
    vector<double> trialAwParameters;

    //mc variables
    vector<int> switchTriangles(4); //nodes making triangles to switch
    mcAcceptedMoves=0; //counter for accepted moves
    if(localGeomOpt) geomOptRejectCount=0; //counter if rejection due to edge intersection
    mcTargetReached=false;
    mcProposedMoves=mcMaxMoves; //if target is not met will have max mc moves
    if(periodic){
        for(int move=0; move<mcMaxMoves; ++move){
            trialPVector=pVector;
            trialPMatrix=pMatrix;
            switchTriangles=pickRandomTrianglePairPeriodic();
            calculateTrialPPeriodic(switchTriangles, trialPVector, trialPMatrix);
            trialAwParameters=calculateAboavWeaireFit(trialPVector,trialPMatrix);
            trialMcEnergy=mcEnergyFunctional(trialAwParameters,trialPVector);
            testCounter=move+1;
            acceptTrialMove=evaluateMetropolisCondition(trialMcEnergy,mcEnergy);
            if(acceptTrialMove) mcTargetReached=acceptDualSwitchPeriodic(switchTriangles,trialPVector,trialPMatrix,trialMcEnergy,trialAwParameters);
            if(mcTargetReached){
                mcProposedMoves=move+1;
                break;
            }
            if((move+1)%1000000==0){
                checkFidelity();
                if (!consistent) break;
            }
        }
    }
    else{
        for(int move=0; move<mcMaxMoves; ++move){
            trialPVector=pVector;
            trialPMatrix=pMatrix;
            switchTriangles=pickRandomTrianglePairAperiodic();
//            //***** HACKED FOR DEVELOPMENT ******
//            switchTriangles=pickDefinedTrianglePair(move);
//            //****** HACK END ******
            calculateTrialPAperiodic(switchTriangles, trialPVector, trialPMatrix);
            trialAwParameters=calculateAboavWeaireFit(trialPVector,trialPMatrix);
            trialMcEnergy=mcEnergyFunctional(trialAwParameters,trialPVector);
            testCounter=move;
            acceptTrialMove=evaluateMetropolisCondition(trialMcEnergy,mcEnergy);
//            //****** HACK BEGIN ******
//            acceptTrialMove=true;
//            //******HACK END ******
            if(acceptTrialMove) mcTargetReached=acceptDualSwitchAperiodic(switchTriangles,trialPVector,trialPMatrix,trialMcEnergy,trialAwParameters);
            if(mcTargetReached){
                mcProposedMoves=move+1;
                break;
            }
            if((move+1)%1000000==0){
                checkFidelity();
                if (!consistent) break;
            }
        }
    }

    return;
}

vector<int> Network::pickRandomTrianglePairPeriodic() {
    //pick four nodes that make up two triangles in dual, making sure doesn't violate ring size limits

    //indices 0,1 as connected pair, 2,3 as bridge.
    //0,1 can't both be on edge
    //0,1 when decremented can't be less than min ring size. 2,3 when incremented can't exceed max ring size
    int ref0, ref1, ref2, ref3;
    vector<int> trianglePair(4);

    bool picked=false;
    vector<int> commonNodes; //shared nodes, should be two
    do{//loop until find valid pair
        ref0=pickRandomNode();
        if(nodes[ref0].size>minRingSize){//select node within size limit
            ref1=nodes[ref0].connections[pickRandomConnection(nodes[ref0].size)];
            if(nodes[ref1].size>minRingSize){//select node within size limit
                commonNodes=getCommonValuesBetweenVectors(nodes[ref0].connections, nodes[ref1].connections);
                if(commonNodes.size()==2){//should always be two
                    ref2=commonNodes[0];
                    ref3=commonNodes[1];
                    if(nodes[ref2].size<maxRingSize && nodes[ref3].size<maxRingSize) picked=true;
                }
            }
        }
    }while(!picked);

    trianglePair[0]=ref0;
    trianglePair[1]=ref1;
    trianglePair[2]=ref2;
    trianglePair[3]=ref3;

    return trianglePair;
}

vector<int> Network::pickRandomTrianglePairAperiodic() {
    //pick four nodes that make up two triangles in dual, making sure connected pair not on edge and doesn't violate ring size limits

    //indices 0,1 as connected pair, 2,3 as bridge.
    //0,1 can't both be on edge
    //0,1 when decremented can't be less than min ring size. 2,3 when incremented can't exceed max ring size
    int ref0, ref1, ref2, ref3;
    vector<int> trianglePair(4);

    bool picked=false;
    vector<int> commonNodes; //shared nodes, should be two
    do{//loop until find valid pair
        ref0=pickRandomNode();
        if(!nodes[ref0].edge && nodes[ref0].size>minRingSize){//select node not on edge and within size limit
            ref1=nodes[ref0].connections[pickRandomConnection(nodes[ref0].size)];
            if(nodes[ref1].size>minRingSize){//select node within size limit
                commonNodes=getCommonValuesBetweenVectors(nodes[ref0].connections, nodes[ref1].connections);
                if(commonNodes.size()==2){//should always be two
                    ref2=commonNodes[0];
                    ref3=commonNodes[1];
                    if(nodes[ref2].size<maxRingSize && nodes[ref3].size<maxRingSize) picked=true;
                }
            }
        }
    }while(!picked);

    trianglePair[0]=ref0;
    trianglePair[1]=ref1;
    trianglePair[2]=ref2;
    trianglePair[3]=ref3;

    return trianglePair;
}

vector<int> Network::pickDefinedTrianglePair(int m) {
    //**** FOR DEVELOPMENT, PICK TRIANGLE PAIR BASED ON MOVE ****
    //**** NO CONSISTENCY CHECKING ****

    int ref0, ref1, ref2, ref3;
    vector<int> trianglePair(4);

    bool picked=false;
    vector<int> commonNodes; //shared nodes, should be two

    int moveI=m/4, moveJ=m%4, moveK=m/16;
    int moveStart;

    if(moveI==0) moveStart=250;
    else if(moveI==1) moveStart=251;
    else if(moveI==2) moveStart=272;
    else if(moveI==3) moveStart=184;
    else if(moveI==4) moveStart=251;
    else if(moveI==5) moveStart=253;
    else if(moveI==6) moveStart=226;

    if(moveI==0){
        if(moveJ==0){
            ref0=moveStart;
            ref1=ref0-22;
        }
        else if(moveJ==1){
            ref0=moveStart+1;
            ref1=ref0-22;
        }
        else if(moveJ==2){
            ref0=moveStart-23;
            ref1=ref0-21;
        }
        else if(moveJ==3){
            ref0=moveStart-22;
            ref1=ref0-21;
        }
    }
    if(moveI==1){
        if(moveJ==0){
            ref0=moveStart;
            ref1=ref0-23;
        }
        else if(moveJ==1){
            ref0=moveStart+1;
            ref1=ref0-23;
        }
        else if(moveJ==2){
            ref0=moveStart-24;
            ref1=ref0-22;
        }
        else if(moveJ==3){
            ref0=moveStart-23;
            ref1=ref0-22;
        }
    }
    if(moveI==2){
        if(moveJ==0){
            ref0=moveStart;
            ref1=ref0-21;
        }
        else if(moveJ==1){
            ref0=moveStart+1;
            ref1=ref0-21;
        }
        else if(moveJ==2){
            ref0=moveStart-23;
            ref1=ref0-23;
        }
        else if(moveJ==3){
            ref0=moveStart-22;
            ref1=ref0-23;
        }
    }
    if(moveI==3){
        if(moveJ==0){
            ref0=moveStart;
            ref1=ref0+22;
        }
        else if(moveJ==1){
            ref0=moveStart-1;
            ref1=ref0+22;
        }
        else if(moveJ==2){
            ref0=moveStart+24;
            ref1=ref0+22;
        }
        else if(moveJ==3){
            ref0=moveStart+23;
            ref1=ref0+22;
        }
    }

//    if(moveI==0) moveStart=250;
//    else if(moveI==1) moveStart=252;
//    else if(moveI==2) moveStart=205;
//    else if(moveI==3) moveStart=207;
//    else if(moveI==4) moveStart=251;
//    else if(moveI==5) moveStart=253;
//    else if(moveI==6) moveStart=226;
//
//    if(moveK==0){
//        if(moveJ==0){
//            ref0=moveStart;
//            ref1=ref0-22;
//        }
//        else if(moveJ==1){
//            ref0=moveStart+1;
//            ref1=ref0-22;
//        }
//        else if(moveJ==2){
//            ref0=moveStart-23;
//            ref1=ref0-21;
//        }
//        else if(moveJ==3){
//            ref0=moveStart-22;
//            ref1=ref0-21;
//        }
//    }
//    if(moveI==4 or moveI==5){
//        if(moveJ==0){
//            ref0=moveStart;
//            ref1=ref0-23;
//        }
//        if(moveJ==1){
//            ref0=moveStart+1;
//            ref1=ref0-23;
//        }
//        if(moveJ==2){
//            ref0=moveStart-24;
//            ref1=ref0-22;
//        }
//        if(moveJ==3){
//            ref0=moveStart-23;
//            ref1=ref0-22;
//        }
//    }
//    if(moveI==6){
//        if(moveJ==0){
//            ref0=moveStart;
//            ref1=ref0-21;
//        }
//        if(moveJ==1){
//            ref0=moveStart+2;
//            ref1=ref0+1;
//        }
//        if(moveJ==2){
//            ref0=moveStart+27;
//            ref1=ref0+1;
//        }
//    }
//    if(m==27){
//        ref0=187;
//        ref1=186;
//    }
//    if(m==28){
//        ref0=185;
//        ref1=209;
//    }

    ref0=20;
    ref1=13;

    commonNodes=getCommonValuesBetweenVectors(nodes[ref0].connections, nodes[ref1].connections);
    ref2=commonNodes[0];
    ref3=commonNodes[1];

    trianglePair[0]=ref0;
    trianglePair[1]=ref1;
    trianglePair[2]=ref2;
    trianglePair[3]=ref3;

    consoleVector(trianglePair);
    return trianglePair;
}

void Network::calculateTrialPPeriodic(vector<int> &triangles, vector<int> &trialPVector, vector<vector<int> > &trialPMatrix) {
    //calculate trial p vector and matrix for proposed dual switch

    //current vs trial sizes and whether edge
    vector<int> currSizeIndex(4), trialSizeIndex(4);
    for(int i=0; i<2; ++i){
        currSizeIndex[i]=nodes[triangles[i]].sizeIndex;
        trialSizeIndex[i]=currSizeIndex[i]-1;
    }
    for(int i=2; i<4; ++i){
        currSizeIndex[i]=nodes[triangles[i]].sizeIndex;
        trialSizeIndex[i]=currSizeIndex[i]+1;
    }

    //p vector
    for(int i=0; i<4; ++i){
        trialPVector[currSizeIndex[i]]=--trialPVector[currSizeIndex[i]];
        trialPVector[trialSizeIndex[i]]=++trialPVector[trialSizeIndex[i]];
    }

    //p matrix
    //consider only connections to nodes not in triangles first
    vector<int> cnxs;
    int nbSizeIndex;
    for(int i=0; i<4; ++i){
        cnxs=nodes[triangles[i]].connections;
        removeValuesFromVectorByRef(cnxs,triangles);
        for(int j=0; j<cnxs.size(); ++j){//loop over remaining connections
            if(!nodes[cnxs[j]].edge){
                nbSizeIndex=nodes[cnxs[j]].sizeIndex;
                trialPMatrix[currSizeIndex[i]][nbSizeIndex]=--trialPMatrix[currSizeIndex[i]][nbSizeIndex];
                trialPMatrix[nbSizeIndex][currSizeIndex[i]]=--trialPMatrix[nbSizeIndex][currSizeIndex[i]];
                trialPMatrix[trialSizeIndex[i]][nbSizeIndex]=++trialPMatrix[trialSizeIndex[i]][nbSizeIndex];
                trialPMatrix[nbSizeIndex][trialSizeIndex[i]]=++trialPMatrix[nbSizeIndex][trialSizeIndex[i]];
            }
        }
    }
    //consider connections for nodes in triangles
    for(int i=0; i<2; ++i){
        for(int j=i+1; j<4; ++j){//break connections
            trialPMatrix[currSizeIndex[i]][currSizeIndex[j]]=--trialPMatrix[currSizeIndex[i]][currSizeIndex[j]];
            trialPMatrix[currSizeIndex[j]][currSizeIndex[i]]=--trialPMatrix[currSizeIndex[j]][currSizeIndex[i]];
        }
        for(int j=2; j<4; ++j){//make connections except 2->3
            trialPMatrix[trialSizeIndex[i]][trialSizeIndex[j]]=++trialPMatrix[trialSizeIndex[i]][trialSizeIndex[j]];
            trialPMatrix[trialSizeIndex[j]][trialSizeIndex[i]]=++trialPMatrix[trialSizeIndex[j]][trialSizeIndex[i]];
        }
    }
    trialPMatrix[trialSizeIndex[2]][trialSizeIndex[3]]=++trialPMatrix[trialSizeIndex[2]][trialSizeIndex[3]];
    trialPMatrix[trialSizeIndex[3]][trialSizeIndex[2]]=++trialPMatrix[trialSizeIndex[3]][trialSizeIndex[2]];

    return;
}

void Network::calculateTrialPAperiodic(vector<int> &triangles, vector<int> &trialPVector, vector<vector<int> > &trialPMatrix) {
    //calculate trial p vector and matrix for proposed dual switch

    //current vs trial sizes and whether edge
    vector<int> currSizeIndex(4), trialSizeIndex(4);
    vector<bool> edge(4);
    for(int i=0; i<2; ++i){
        currSizeIndex[i]=nodes[triangles[i]].sizeIndex;
        trialSizeIndex[i]=currSizeIndex[i]-1;
        edge[i]=nodes[triangles[i]].edge;
    }
    for(int i=2; i<4; ++i){
        currSizeIndex[i]=nodes[triangles[i]].sizeIndex;
        trialSizeIndex[i]=currSizeIndex[i]+1;
        edge[i]=nodes[triangles[i]].edge;
    }

    //p vector
    for(int i=0; i<4; ++i){
        if(!edge[i]){
            trialPVector[currSizeIndex[i]]=--trialPVector[currSizeIndex[i]];
            trialPVector[trialSizeIndex[i]]=++trialPVector[trialSizeIndex[i]];
        }
    }

    //p matrix
    //consider only connections to nodes not in triangles first
    vector<int> cnxs;
    int nbSizeIndex;
    for(int i=0; i<4; ++i){
        if(!edge[i]){//if not on edge
            cnxs=nodes[triangles[i]].connections;
            removeValuesFromVectorByRef(cnxs,triangles);
            for(int j=0; j<cnxs.size(); ++j){//loop over remaining connections
                if(!nodes[cnxs[j]].edge){
                    nbSizeIndex=nodes[cnxs[j]].sizeIndex;
                    trialPMatrix[currSizeIndex[i]][nbSizeIndex]=--trialPMatrix[currSizeIndex[i]][nbSizeIndex];
                    trialPMatrix[nbSizeIndex][currSizeIndex[i]]=--trialPMatrix[nbSizeIndex][currSizeIndex[i]];
                    trialPMatrix[trialSizeIndex[i]][nbSizeIndex]=++trialPMatrix[trialSizeIndex[i]][nbSizeIndex];
                    trialPMatrix[nbSizeIndex][trialSizeIndex[i]]=++trialPMatrix[nbSizeIndex][trialSizeIndex[i]];
                }
            }
        }
    }
    //consider connections for nodes in triangles
    for(int i=0; i<2; ++i){
        if(!edge[i]){
            for(int j=i+1; j<4; ++j){//break connections
                if(!edge[j]){
                    trialPMatrix[currSizeIndex[i]][currSizeIndex[j]]=--trialPMatrix[currSizeIndex[i]][currSizeIndex[j]];
                    trialPMatrix[currSizeIndex[j]][currSizeIndex[i]]=--trialPMatrix[currSizeIndex[j]][currSizeIndex[i]];
                }
            }
            for(int j=2; j<4; ++j){//make connections except 2->3
                if(!edge[j]){
                    trialPMatrix[trialSizeIndex[i]][trialSizeIndex[j]]=++trialPMatrix[trialSizeIndex[i]][trialSizeIndex[j]];
                    trialPMatrix[trialSizeIndex[j]][trialSizeIndex[i]]=++trialPMatrix[trialSizeIndex[j]][trialSizeIndex[i]];
                }
            }
        }
    }
    if(!edge[2] && !edge[3]){//2->3 connection
        trialPMatrix[trialSizeIndex[2]][trialSizeIndex[3]]=++trialPMatrix[trialSizeIndex[2]][trialSizeIndex[3]];
        trialPMatrix[trialSizeIndex[3]][trialSizeIndex[2]]=++trialPMatrix[trialSizeIndex[3]][trialSizeIndex[2]];
    }

    return;
}

vector<double> Network::calculateAboavWeaireFit(vector<int> &pVec, vector<vector<int> > &pMat) {
    //calculate 6(n-6) and nmn and perform linear regression
    vector<Crd2d> fitData;
    fitData.clear();

    double mn,norm;///aw alpha and mean ring size about each ring,normalisation factor
    //find x,y for rings present in system
    for(int i=0; i<nRingSizes; ++i){
        if(pVec[i]>0){
            mn=0.0;
            norm=1.0/accumulate(pMat[i].begin(), pMat[i].end(), 0.0);
            for (int j=0; j<nRingSizes; ++j) mn=mn+(minRingSize+j)*norm*pMat[i][j];
            fitData.push_back(Crd2d(6.0*(minRingSize+i-6.0),mn*(minRingSize+i)));
        }
    }

    //perform linear regression if not just one ring type
    vector<double> fitParameters(3,0.0); //alpha, mu and rsq
    if(fitData.size()>1){
        fitParameters=leastSquaresLinearRegression(fitData);
        fitParameters[0]=1.0-fitParameters[0]; //alpha from gradient
        fitParameters[1]=fitParameters[1]-36.0; //mu from intercept // rsq unchanged
    }
    return fitParameters;
}

double Network::mcEnergyFunctional(vector<double> &awParams, vector<int> &pVec) {
    //evaluate energy for monte carlo condition from difference between current and target network properties

    //difference of alpha to target
    double e=fabs(awParams[0]-targetAlpha)*mcAlphaScaleFactor;

    //difference of fitted mu to mu
    e=e+fabs(awParams[1]-targetMu)*rTargetMu;

    //difference of ring statistics to target
    double norm=1.0/accumulate(pVec.begin(), pVec.end(), 0.0);
    for(int i=0; i<nRingSizes; ++i){
        e=e+fabs(norm*pVec[i]-targetPVector[i])*rTargetPVector[i];
    }

    return e;
}

bool Network::evaluateMetropolisCondition(double &trialEnergy, double &currEnergy) {
    //condition to accept/reject mc move
    bool accept=true;
    double deltaE=trialEnergy-currEnergy;
    if(deltaE<0.0) accept=true;
    else{
        double condition=exp(-deltaE*rMcTemperature);
        if(condition>metropolisRandomNum()) accept=true;
        else accept=false;
    }

    return accept;
}

bool Network::acceptDualSwitchPeriodic(vector<int> &switchTriangles, vector<int> &trialPVec, vector<vector<int> > &trialPMat, double &trialMcEnergy, vector<double> &trialAwParams) {
    //enact dual switch on nodes and update trial-> current variables

    //make/break connections
    nodes[switchTriangles[0]].breakConnection(switchTriangles[1]);
    nodes[switchTriangles[1]].breakConnection(switchTriangles[0]);
    nodes[switchTriangles[2]].makeConnection(switchTriangles[3]);
    nodes[switchTriangles[3]].makeConnection(switchTriangles[2]);

    //locally minimise
    if(localGeomOpt){
        int optimisationStatus=localMinimisationPeriodic(switchTriangles);
        //if failed reverse make/break and don't update vectors
        if(optimisationStatus==0){
            nodes[switchTriangles[0]].makeConnection(switchTriangles[1]);
            nodes[switchTriangles[1]].makeConnection(switchTriangles[0]);
            nodes[switchTriangles[2]].breakConnection(switchTriangles[3]);
            nodes[switchTriangles[3]].breakConnection(switchTriangles[2]);
            geomOptRejectCount=++geomOptRejectCount;
            return false;
        }
    }

    //update vectors
    pVector=trialPVec;
    pMatrix=trialPMat;
    mcEnergy=trialMcEnergy;
    aboavWeaireParams=trialAwParams;
    mcAcceptedMoves=++mcAcceptedMoves;

    cout<<testCounter<<" "<<aboavWeaireParams[0]<<" "<<aboavWeaireParams[1]<<" "<<aboavWeaireParams[2]<<" "<<pVector[1]<<" "<<mcEnergy<<endl;

    if(mcEnergy<=mcConvergence) return true;
    else return false;
}

bool Network::acceptDualSwitchAperiodic(vector<int> &switchTriangles, vector<int> &trialPVec, vector<vector<int> > &trialPMat, double &trialMcEnergy, vector<double> &trialAwParams) {
    //enact dual switch on nodes and update trial-> current variables

    //make/break connections
    nodes[switchTriangles[0]].breakConnection(switchTriangles[1]);
    nodes[switchTriangles[1]].breakConnection(switchTriangles[0]);
    nodes[switchTriangles[2]].makeConnection(switchTriangles[3]);
    nodes[switchTriangles[3]].makeConnection(switchTriangles[2]);

    //locally minimise
    if(localGeomOpt){
        int optimisationStatus=localMinimisationAperiodic(switchTriangles);
        //if failed reverse make/break and don't update vectors
        if(optimisationStatus==0){
            nodes[switchTriangles[0]].makeConnection(switchTriangles[1]);
            nodes[switchTriangles[1]].makeConnection(switchTriangles[0]);
            nodes[switchTriangles[2]].breakConnection(switchTriangles[3]);
            nodes[switchTriangles[3]].breakConnection(switchTriangles[2]);
            geomOptRejectCount=++geomOptRejectCount;
            return false;
        }
    }

    //update vectors
    pVector=trialPVec;
    pMatrix=trialPMat;
    mcEnergy=trialMcEnergy;
    aboavWeaireParams=trialAwParams;
    mcAcceptedMoves=++mcAcceptedMoves;

//    cout<<testCounter<<" "<<aboavWeaireParams[0]<<" "<<aboavWeaireParams[1]<<" "<<aboavWeaireParams[2]<<" "<<mcEnergy<<endl;
    if(mcEnergy<=mcConvergence) return true;
    else return false;
}

double Network::metropolisRandomNum() {
    //random between 0->1
    return zeroOneDistribution(randomGenerator);
}

int Network::pickRandomNode() {
    //random between 0->nNodes-1
    return nodeDistribution(randomGenerator);
}

int Network::pickRandomConnection(int nCnxs) {
    //random between 0->nCnxs
    double rand=zeroOneDistribution(randomGenerator)*(nCnxs-1);
    return round(rand);
}

void Network::findNodeRings() {
    //find rings of nodes in dual graph around each central node

    vector<int> ring; //ring around central node
    for(int i=0; i<nNodes; ++i){//loop over all nodes in graph not on edge
        if(!nodes[i].edge){
            ring=findNodeRing(i);
            nodeRings.push_back(Ring(nodes[i].size,ring,i));
        }
    }
    nRings=nodeRings.size(); //number of rings in network

    return;
}

int Network::localMinimisationPeriodic(vector<int> &switchTriangles) {
    //minimse within 1 connection of triangle pair

    //get next two topological shells
    vector<int> minimisationRegion=switchTriangles;
    vector<int> nextShell=nextTopologicalShell(switchTriangles);
    vector<int> fixedRegion=nextTopologicalShell(nextShell,switchTriangles);

    //minimisation region contains all shells, fixed region just outermost shell
    addValuesToVectorByRef(minimisationRegion,nextShell);
    addValuesToVectorByRef(minimisationRegion,fixedRegion);

    //set up global maps
    int nMinNodes=minimisationRegion.size();
    map<int,int> globalToLocalMap;
    vector<Crd2d> localCoordinates(nMinNodes);
    for(int i=0; i<nMinNodes; ++i){
        globalToLocalMap[minimisationRegion[i]]=i;
        localCoordinates[i]=nodes[minimisationRegion[i]].coordinate;
    }

    //make unique list of harmonic pairs
    vector<Pair> localHarmonicPairs;
    vector<double> localHarmonicR0;
    localHarmonicPairs.clear();
    localHarmonicR0.clear();
    int globalA, globalB, indexA, indexB;
    for(int i=0; i<nMinNodes; ++i){
        globalA=minimisationRegion[i];
        for(int j=0; j<nodes[globalA].size; ++j){
            globalB=nodes[globalA].connections[j];
            if(globalToLocalMap.count(globalB)) {//only include if in minimisation region
                if(globalA<globalB){//prevent double counting
                    localHarmonicPairs.push_back(Pair(globalToLocalMap[globalA],globalToLocalMap[globalB]));
                    indexA=nodes[globalA].sizeIndex;
                    indexB=nodes[globalB].sizeIndex;
                    localHarmonicR0.push_back(harmonicR0Matrix[indexA][indexB]);
                }
            }
        }
    }

    //get list of line intersections to prevent overlapping
    int ringSize;
    vector<int> ring;
    vector<DoublePair> localLinePairs;
    DoublePair linePair;
    int lineA, lineB, lineC, lineD;
    map<string,int> pairIDMap;
    string pairID;
    for(int i=0; i<nMinNodes-fixedRegion.size(); ++i){//find rings around each node not in fixed region
        ringSize=nodes[minimisationRegion[i]].size;
        ring=findNodeRing(minimisationRegion[i]);
        ring.push_back(ring[0]); //make ring periodic
        lineA=globalToLocalMap[minimisationRegion[i]];
        for(int j=0; j<ringSize; ++j){
            lineB=globalToLocalMap[ring[j]];
            for(int k=0; k<ringSize; ++k){
                lineC=globalToLocalMap[ring[k]];
                lineD=globalToLocalMap[ring[k+1]];
                if(lineC!=lineB && lineD!=lineB){//neglect lines containing node B
                    linePair=DoublePair(lineA,lineB,lineC,lineD);
                    //ensure each pair is only used once
                    linePair.sort();
                    pairID=linePair.getID();
                    if(!pairIDMap.count(pairID)){
                        localLinePairs.push_back(linePair);
                        pairIDMap[pairID]=1;
                    }
                }
            }
        }
    }

    //convert fixed region to local
    for(int i=0; i<fixedRegion.size(); ++i) fixedRegion[i]=globalToLocalMap[fixedRegion[i]];

    //minimise
    vector<Trio> emptyAngles; //placeholder
    emptyAngles.clear();
    HarmonicPeriodicGO optimiser;
    optimiser.setCoordinates(localCoordinates);
    optimiser.setSystemParameters(localHarmonicPairs,emptyAngles,fixedRegion,localLinePairs);
    optimiser.setPotentialParameters(harmonicK,localHarmonicR0);
    optimiser.setOptisationParameters(localGeomOptCC,geomOptLineSearchInc,localGeomOptMaxIt,true);
    optimiser.setPeriodicBoundary(periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
    int status=optimiser.steepestDescent();
    localCoordinates=optimiser.getMinimisedCoordinates();

    //update global coordinates
    for(int i=0; i<nMinNodes; ++i) nodes[minimisationRegion[i]].coordinate=localCoordinates[i];

    return status;
}

int Network::localMinimisationAperiodic(vector<int> &switchTriangles) {
    //minimse within 1 connection of triangle pair

    //get next two topological shells
    vector<int> minimisationRegion=switchTriangles;
    vector<int> nextShell=nextTopologicalShell(switchTriangles);
    vector<int> fixedRegion=nextTopologicalShell(nextShell,switchTriangles);

    //minimisation region contains all shells, fixed region just outermost shell
    addValuesToVectorByRef(minimisationRegion,nextShell);
    addValuesToVectorByRef(minimisationRegion,fixedRegion);

    //set up global maps
    int nMinNodes=minimisationRegion.size();
    map<int,int> globalToLocalMap;
    vector<Crd2d> localCoordinates(nMinNodes);
    for(int i=0; i<nMinNodes; ++i){
        globalToLocalMap[minimisationRegion[i]]=i;
        localCoordinates[i]=nodes[minimisationRegion[i]].coordinate;
    }

    //make unique list of harmonic pairs
    vector<Pair> localHarmonicPairs;
    vector<double> localHarmonicR0;
    localHarmonicPairs.clear();
    localHarmonicR0.clear();
    int globalA, globalB, indexA, indexB;
    for(int i=0; i<nMinNodes; ++i){
        globalA=minimisationRegion[i];
        for(int j=0; j<nodes[globalA].size; ++j){
            globalB=nodes[globalA].connections[j];
            if(globalToLocalMap.count(globalB)) {//only include if in minimisation region
                if(globalA<globalB){//prevent double counting
                    localHarmonicPairs.push_back(Pair(globalToLocalMap[globalA],globalToLocalMap[globalB]));
                    if(nodes[globalA].edge) indexA=index6;
                    else indexA=nodes[globalA].sizeIndex;
                    if(nodes[globalB].edge) indexB=index6;
                    else indexB=nodes[globalB].sizeIndex;
                    localHarmonicR0.push_back(harmonicR0Matrix[indexA][indexB]);
                }
            }
        }
    }

    //get list of line intersections to prevent overlapping
    int ringSize;
    vector<int> ring;
    vector<DoublePair> localLinePairs;
    DoublePair linePair;
    int lineA, lineB, lineC, lineD;
    map<string,int> pairIDMap;
    string pairID;
    for(int i=0; i<nMinNodes-fixedRegion.size(); ++i){//find rings around each node not in fixed region
        if(!nodes[minimisationRegion[i]].edge){//not on edge
            ringSize=nodes[minimisationRegion[i]].size;
            ring=findNodeRing(minimisationRegion[i]);
            ring.push_back(ring[0]); //make ring periodic
            lineA=globalToLocalMap[minimisationRegion[i]];
            for(int j=0; j<ringSize; ++j){
                lineB=globalToLocalMap[ring[j]];
                for(int k=0; k<ringSize; ++k){
                    lineC=globalToLocalMap[ring[k]];
                    lineD=globalToLocalMap[ring[k+1]];
                    if(lineC!=lineB && lineD!=lineB){//neglect lines containing node B
                        linePair=DoublePair(lineA,lineB,lineC,lineD);
                        //ensure each pair is only used once
                        linePair.sort();
                        pairID=linePair.getID();
                        if(!pairIDMap.count(pairID)){
                            localLinePairs.push_back(linePair);
                            pairIDMap[pairID]=1;
                        }
                    }
                }

            }
        }
    }

    //convert fixed region to local
    for(int i=0; i<fixedRegion.size(); ++i) fixedRegion[i]=globalToLocalMap[fixedRegion[i]];

    //minimise
    vector<Trio> emptyAngles; //placeholder
    emptyAngles.clear();
    HarmonicAperiodicGO optimiser;
//    KeatingAperiodicBondOnlyGO optimiser;
    optimiser.setCoordinates(localCoordinates);
    optimiser.setSystemParameters(localHarmonicPairs,emptyAngles,fixedRegion,localLinePairs);
    optimiser.setPotentialParameters(harmonicK,localHarmonicR0);
    optimiser.setOptisationParameters(localGeomOptCC,geomOptLineSearchInc,localGeomOptMaxIt,true);
    int status=optimiser.steepestDescent();
    localCoordinates=optimiser.getMinimisedCoordinates();


    //update global coordinates
    for(int i=0; i<nMinNodes; ++i) nodes[minimisationRegion[i]].coordinate=localCoordinates[i];

    return status;
}

void Network::globalMinimisationPeriodic() {
    //minimse entire network

    //make list of coordinates
    vector<Crd2d> globalCoordinates(nNodes);
    for(int i=0; i<nNodes; ++i) globalCoordinates[i]=nodes[i].coordinate;

    //make unique list of harmonic pairs
    vector<Pair> globalHarmonicPairs;
    vector<double> globalHarmonicR0;
    globalHarmonicPairs.clear();
    globalHarmonicR0.clear();
    int globalA, globalB, indexA, indexB;
    for(int i=0; i<nNodes; ++i){
        globalA=i;
        for(int j=0; j<nodes[globalA].size; ++j){
            globalB=nodes[globalA].connections[j];
            if(globalA<globalB){//prevent double counting
                globalHarmonicPairs.push_back(Pair(globalA,globalB));
                indexA=nodes[globalA].sizeIndex;
                indexB=nodes[globalB].sizeIndex;
                globalHarmonicR0.push_back(harmonicR0Matrix[indexA][indexB]);
            }
        }
    }

    //get list of unique line intersections
    vector<DoublePair> linePairs;
    DoublePair linePair;
    int lineA, lineB, lineC, lineD;
    map<string,int> pairIDMap;
    string pairID;
    linePairs.clear();
    for(int i=0; i<nRings; ++i){//loop over ring centres
        lineA=nodeRings[i].id;
        for(int j=0; j<nodeRings[i].size; ++j){//loop over central->chain lines
            lineB=nodeRings[i].chain[j];
            for(int k=0; k<nodeRings[i].size;++k){//loop over chain->chain lines
                lineC=nodeRings[i].chain[k];
                lineD=nodeRings[i].chain[k+1];
                if(lineC!=lineB && lineD!=lineB){//neglect chain->chain lines containing B
                    linePair=DoublePair(lineA,lineB,lineC,lineD);
                    //ensure each pair is only used once
                    linePair.sort();
                    pairID=linePair.getID();
                    if(!pairIDMap.count(pairID)){
                        linePairs.push_back(linePair);
                        pairIDMap[pairID]=1;
                    }
                }
            }
        }
    }

    //no fixed region
    vector<int> fixedRegion;
    fixedRegion.clear();

    //minimise
    vector<Trio> emptyAngles; //placeholder
    emptyAngles.clear();
    HarmonicPeriodicGO optimiser;
    optimiser.setCoordinates(globalCoordinates);
    optimiser.setSystemParameters(globalHarmonicPairs,emptyAngles,fixedRegion,linePairs);
    optimiser.setPotentialParameters(harmonicK,globalHarmonicR0);
    optimiser.setOptisationParameters(globalGeomOptCC,geomOptLineSearchInc,globalGeomOptMaxIt,true);
    optimiser.setPeriodicBoundary(periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
    int status=optimiser.steepestDescent();
    globalCoordinates=optimiser.getMinimisedCoordinates();


    //update coordinates
    for(int i=0; i<nNodes; ++i) nodes[i].coordinate=globalCoordinates[i];

    //get global energy and iterations
    geomOptEnergy=optimiser.getEnergy();
    geomOptIterations=optimiser.getIterations();

    return;
}

void Network::globalMinimisationAperiodic() {
    //minimse entire network

    //make list of coordinates
    vector<Crd2d> globalCoordinates(nNodes);
    for(int i=0; i<nNodes; ++i) globalCoordinates[i]=nodes[i].coordinate;

    //make unique list of harmonic pairs
    vector<Pair> globalHarmonicPairs;
    vector<double> globalHarmonicR0;
    globalHarmonicPairs.clear();
    globalHarmonicR0.clear();
    int globalA, globalB, indexA, indexB;
    for(int i=0; i<nNodes; ++i){
        globalA=i;
        for(int j=0; j<nodes[globalA].size; ++j){
            globalB=nodes[globalA].connections[j];
            if(globalA<globalB){//prevent double counting
                globalHarmonicPairs.push_back(Pair(globalA,globalB));
                if(nodes[globalA].edge) indexA=index6;
                else indexA=nodes[globalA].sizeIndex;
                if(nodes[globalB].edge) indexB=index6;
                else indexB=nodes[globalB].sizeIndex;
                globalHarmonicR0.push_back(harmonicR0Matrix[indexA][indexB]);
            }
        }
    }

    //get list of unique line intersections
    vector<DoublePair> linePairs;
    DoublePair linePair;
    int lineA, lineB, lineC, lineD;
    map<string,int> pairIDMap;
    string pairID;
    linePairs.clear();
    for(int i=0; i<nRings; ++i){//loop over ring centres
        lineA=nodeRings[i].id;
        for(int j=0; j<nodeRings[i].size; ++j){//loop over central->chain lines
            lineB=nodeRings[i].chain[j];
            for(int k=0; k<nodeRings[i].size;++k){//loop over chain->chain lines
                lineC=nodeRings[i].chain[k];
                lineD=nodeRings[i].chain[k+1];
                if(lineC!=lineB && lineD!=lineB){//neglect chain->chain lines containing B
                    linePair=DoublePair(lineA,lineB,lineC,lineD);
                    //ensure each pair is only used once
                    linePair.sort();
                    pairID=linePair.getID();
                    if(!pairIDMap.count(pairID)){
                        linePairs.push_back(linePair);
                        pairIDMap[pairID]=1;
                    }
                }
            }
        }
    }

    //no fixed region
    vector<int> fixedRegion;
    fixedRegion.clear();

    //minimise
    vector<Trio> emptyAngles; //placeholder
    emptyAngles.clear();
    HarmonicAperiodicGO optimiser;
//    KeatingAperiodicBondOnlyGO optimiser;
    optimiser.setCoordinates(globalCoordinates);
    optimiser.setSystemParameters(globalHarmonicPairs,emptyAngles,fixedRegion,linePairs);
    optimiser.setPotentialParameters(harmonicK,globalHarmonicR0);
    optimiser.setOptisationParameters(globalGeomOptCC,geomOptLineSearchInc,globalGeomOptMaxIt,true);
    int status=optimiser.steepestDescent();
    globalCoordinates=optimiser.getMinimisedCoordinates();

    //update coordinates
    for(int i=0; i<nNodes; ++i) nodes[i].coordinate=globalCoordinates[i];

    //get global energy and iterations
    geomOptEnergy=optimiser.getEnergy();
    geomOptIterations=optimiser.getIterations();

    return;
}

vector<int> Network::findNodeRing(int n){
    //find ring around a central node

    int nSize=nodes[n].size;
    vector<int> ring(nSize); //ring around central node
    vector<int> centralCnxs=nodes[n].connections;
    vector<int> outerCnxs;
    ring[0]=centralCnxs[0]; //pick first connection as starting point
    outerCnxs=nodes[ring[0]].connections;
    for(int j=0; j<outerCnxs.size(); ++j){//find second connection in ring
        if(checkVectorForValue(centralCnxs,outerCnxs[j])){
            ring[1]=outerCnxs[j];
            break;
        };
    }
    for(int j=2; j<nSize; ++j){//find remaining connections
        outerCnxs=nodes[ring[j-1]].connections;
        removeValueFromVectorByRef(outerCnxs,ring[j-2]); //remove previous node in chain
        for(int k=0; k<outerCnxs.size(); ++k){
            if(checkVectorForValue(centralCnxs,outerCnxs[k])) {
                ring[j] = outerCnxs[k];
                break;
            }
        }
    }

    return ring;
}

vector<int> Network::nextTopologicalShell(vector<int> &currentShell) {
    //find next topological shell of nodes around central nodes
    vector<int> nextShell;
    nextShell.clear();
    for(int i=0; i<currentShell.size(); ++i) addValuesToVectorByRef(nextShell, nodes[currentShell[i]].connections);
    nextShell=sortVectorRemoveDuplicates(nextShell);
    removeValuesFromVectorByRef(nextShell,currentShell);
    return nextShell;
}

vector<int> Network::nextTopologicalShell(vector<int> &currentShell, vector<int> &prevShell) {
    //find next topological shell of nodes around central nodes
    vector<int> nextShell;
    nextShell.clear();
    for(int i=0; i<currentShell.size(); ++i) addValuesToVectorByRef(nextShell, nodes[currentShell[i]].connections);
    nextShell=sortVectorRemoveDuplicates(nextShell);
    removeValuesFromVectorByRef(nextShell,currentShell);
    removeValuesFromVectorByRef(nextShell,prevShell);
    return nextShell;
}

//###### Construction checking ######

void Network::checkFidelity() {
    //check self consistency of dual and p vector/matrix

    bool selfConsistent=true;
    int nodeRef;
    for(int i=0; i<nNodes; ++i){
        if(nodes[i].size!=nodes[i].connections.size()) selfConsistent=false;
        if(nodes[i].size!=nodes[i].sizeIndex+minRingSize) selfConsistent=false;
        if(selfConsistent){
            for(int j=0; j<nodes[i].size; ++j) if(!checkVectorForValue(nodes[nodes[i].connections[j]].connections,i)) selfConsistent=false;
        }
        else break;
    }
//    cout<<"consistent: "<<selfConsistent<<endl;

    vector<int> checkPVector(nRingSizes,0);
    vector< vector<int> > checkPMatrix(nRingSizes,(vector<int> (nRingSizes,0)));
    //reconstruct p vector/matrix from the dual and check with current p vector/matrix
    int indexI, indexJ;
    for(int i=0; i<nNodes; ++i){
        if(!nodes[i].edge){
            indexI=nodes[i].sizeIndex;
            checkPVector[indexI]=++checkPVector[indexI];
            for(int j=0; j<nodes[i].size; ++j){
                if(!nodes[nodes[i].connections[j]].edge){
                    indexJ=nodes[nodes[i].connections[j]].sizeIndex;
                    checkPMatrix[indexI][indexJ]=++checkPMatrix[indexI][indexJ];
                }
            }
        }
    }
    bool pVectorPass=true, pMatrixPass=true;
    for(int i=0; i<nRingSizes; ++i){
        if(pVector[i]!=checkPVector[i]) pVectorPass=false;
        for(int j=0; j<nRingSizes; ++j) if(pMatrix[i][j]!=checkPMatrix[i][j]) pMatrixPass=false;
    }
//    cout<<"p vector: "<<pVectorPass<<" p matrix: "<<pMatrixPass<<endl;

    if(selfConsistent && pVectorPass && pMatrixPass) consistent=true;
    else consistent=false;

    return;
}

void Network::checkGeometry() {
    //check for node edge overlap which would lead to ring overlap

    //get list of unique line intersections
    vector<DoublePair> linePairs;
    DoublePair linePair;
    int lineA, lineB, lineC, lineD;
    map<string,int> pairIDMap;
    string pairID;
    linePairs.clear();
    for(int i=0; i<nRings; ++i){//loop over ring centres
        lineA=nodeRings[i].id;
        for(int j=0; j<nodeRings[i].size; ++j){//loop over central->chain lines
            lineB=nodeRings[i].chain[j];
            for(int k=0; k<nodeRings[i].size;++k){//loop over chain->chain lines
                lineC=nodeRings[i].chain[k];
                lineD=nodeRings[i].chain[k+1];
                if(lineC!=lineB && lineD!=lineB){//neglect chain->chain lines containing B
                    linePair=DoublePair(lineA,lineB,lineC,lineD);
                    //ensure each pair is only used once
                    linePair.sort();
                    pairID=linePair.getID();
                    if(!pairIDMap.count(pairID)){
                        linePairs.push_back(linePair);
                        pairIDMap[pairID]=1;
                    }
                }
            }
        }
    }

    //check if any intersections overlap
    noIntersections=true;
    if(periodic){
        Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
        for(int i=0; i<linePairs.size();++i){
            line1a=nodes[linePairs[i].a].coordinate;
            line1b=minimumImageCrd(line1a,nodes[linePairs[i].b].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            line2a=minimumImageCrd(line1a,nodes[linePairs[i].c].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            line2b=minimumImageCrd(line2a,nodes[linePairs[i].d].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            noIntersections=!properIntersectionLines(line1a,line1b,line2a,line2b);
            if(!noIntersections) break;
        }
    }
    else{
        Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
        for(int i=0; i<linePairs.size();++i){
            line1a=nodes[linePairs[i].a].coordinate;
            line1b=nodes[linePairs[i].b].coordinate;
            line2a=nodes[linePairs[i].c].coordinate;
            line2b=nodes[linePairs[i].d].coordinate;
            noIntersections=!properIntersectionLines(line1a,line1b,line2a,line2b);
            if(!noIntersections) break;
        }
    }
    return;
}

//###### Analysis main ######

void Network::analyse(ofstream &logfile) {
    //control analysis of generated network

    //unoptional analysis
    analyseRingStatistics();
    analyseAboavWeaire();
    //optional analysis
    if(convertToAtomic){
        convertDualToAtomicNetwork();
        if(atomicGeomOpt && periodic) geometryOptimiseAtomicNetworkPeriodic();
        if(atomicGeomOpt && !periodic) geometryOptimiseAtomicNetworkAperiodic();
        analyseAtomicGeometry();
        if(areaLaw) analyseRingAreas();
    }
    if(periodic && spatialRdf) analysePartialSpatialRdfs();
    if(periodic && topoRdf) analysePartialTopologicalRdfs();
    if(periodic && clustering) analyseClusters();
    if(assortativeMix) analyseAssortativity();

    writeFileLine(logfile,name+" analysed");
    return;
}

//###### Analysis individual functions ######

void Network::convertDualToAtomicNetwork() {
    //triangulate nodes to generate atomic configuration

    //make unique vertex list and vertex rings from node rings
    nVertices=0;
    vertices.clear();
    vertexRings.clear();
    Triangle nodeTriangle;
    vector<int> ringPath;
    string triangleID;
    map<string,int> triangleVertexMap; //triangle id to vertex index map
    for(int i=0; i<nRings; ++i){//loop over node rings
        ringPath.clear();
        for(int j=0; j<nodeRings[i].size; ++j){//loop over all node triangles (=vertices) in ring
            nodeTriangle=Triangle(nodeRings[i].id,nodeRings[i].chain[j],nodeRings[i].chain[j+1]);
            nodeTriangle.sort();
            triangleID=nodeTriangle.getID();
            if(triangleVertexMap.count(triangleID)==0){//if vertex not already generated, make vertex
                vertices.push_back(Vertex(nodeTriangle.a, nodeTriangle.b, nodeTriangle.c));
                triangleVertexMap[triangleID]=nVertices;
                nVertices=++nVertices;
            }
            ringPath.push_back(triangleVertexMap[triangleID]);
        }
        vertexRings.push_back(Ring(ringPath.size(),ringPath,i)); //make vertex ring
    }

    //make vertex connections
    for(int i=0; i<vertexRings.size(); ++i){//triply degenerate but vertex only accepts addition if unique
        for(int j=0; j<vertexRings[i].size; ++j){
            vertices[vertexRings[i].chain[j]].addConnection(vertexRings[i].chain[j+1]);
            vertices[vertexRings[i].chain[j+1]].addConnection(vertexRings[i].chain[j]);
        }
    }

    //make vertex coordinates
    vector<Crd2d> nodeTriCrds(3);
    if(!periodic){//aperiodic
        for(int i=0; i<nVertices; ++i){
            for(int j=0; j<3; ++j) nodeTriCrds[j]=nodes[vertices[i].dualNodes[j]].coordinate;
            vertices[i].coordinate=crdCentreOfMass(nodeTriCrds);
        }
    }
    else{//periodic
        for(int i=0; i<nVertices; ++i){
            for(int j=0; j<3; ++j) nodeTriCrds[j]=nodes[vertices[i].dualNodes[j]].coordinate;
            for(int j=1; j<3; ++j) nodeTriCrds[j]=minimumImageCrd(nodeTriCrds[0],nodeTriCrds[j],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            vertices[i].coordinate=crdCentreOfMass(nodeTriCrds);
            applyPeriodicBoundary(vertices[i].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
        }
    }

    return;
}

void Network::analyseRingStatistics() {
    //find ring statistics from p vector and normalise p matrix

    //ring statistics of entire network
    ringStatistics.clear();
    ringStatistics.resize(nRingSizes);
    double normalisation=1.0/accumulate(pVector.begin(), pVector.end(), 0.0);
    for(int i=0; i<nRingSizes; ++i) ringStatistics[i]=normalisation*pVector[i];

    //normalised p matrix (pi)
    piMatrix.clear();
    piMatrix.resize(nRingSizes,(vector<double> (nRingSizes)));
    for(int i=0; i<nRingSizes; ++i){
        normalisation=accumulate(pMatrix[i].begin(), pMatrix[i].end(), 0.0);
        if(normalisation>0.0) normalisation=1.0/normalisation;
        for(int j=0; j<nRingSizes; ++j) piMatrix[i][j]=pMatrix[i][j]*normalisation;
    }

    return;
}

void Network::analyseAboavWeaire() {
    //get aboav weaire parameters - will be the same as calculated in the monte carlo
    aboavWeaireParams=calculateAboavWeaireFit(pVector,pMatrix);
    return;
}

void Network::analysePartialSpatialRdfs() {
    //calculate partial rdfs for each ring size node in dual real space ~com of ring

    //get densities of each node
    spatialRdfDensities.resize(nRingSizes);
    for(int i=0; i<nRingSizes; ++i) spatialRdfDensities[i]=pVector[i]/(periodicBoxX*periodicBoxY);

    //make sure maximum extent is within bounds of network (L/2)
    double maxExtent;
    if(periodicBoxX<periodicBoxY) maxExtent=0.5*periodicBoxX;
    else maxExtent=0.5*periodicBoxY;
    if(spatialRdfExtent<maxExtent) maxExtent=spatialRdfExtent;

    //set up rdfs;
    spatialPartialRdfs.clear(); //set up rdfs for each ring size pair
    for(int i=minRingSize; i<minRingSize+nRingSizes; ++i){
        for(int j=minRingSize; j<minRingSize+nRingSizes; ++j){
            spatialPartialRdfs.push_back(Rdf(spatialRdfBinwidth,spatialRdfExtent,to_string(i)+"-"+to_string(j)));
        }
    }

    //loop over all nodes pairs and add to rdf
    int sizeIndexI, sizeIndexJ, rdfIndexIJ, rdfIndexJI; //size indices and rdfs to add to
    Vec2d vec;
    double d; //distance between points
    for(int i=0; i<nNodes-1; ++i){
        sizeIndexI=nodes[i].sizeIndex;
        for(int j=i+1; j<nNodes; ++j){
            sizeIndexJ=nodes[j].sizeIndex;
            vec=Vec2d(nodes[i].coordinate,nodes[j].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            d=vec.length();
            rdfIndexIJ=sizeIndexI*nRingSizes+sizeIndexJ;
            rdfIndexJI=sizeIndexJ*nRingSizes+sizeIndexI;
            spatialPartialRdfs[rdfIndexIJ].addValue(d);
            spatialPartialRdfs[rdfIndexJI].addValue(d);
        }
    }

    return;
}

void Network::analysePartialTopologicalRdfs() {
    //calculate partial topological rdfs for each ring size node in dual real space ~com of ring

    //set up rdfs;
    topoRdfExtent=++topoRdfExtent; //increment to make inclusive of size in input file
    topoRdfShellSizes.resize(nRingSizes, (vector<int> (int(topoRdfExtent), 0)));
    topoPartialRdfs.clear();
    for(int i=minRingSize; i<minRingSize+nRingSizes; ++i){
        for(int j=minRingSize; j<minRingSize+nRingSizes; ++j){
            topoPartialRdfs.push_back(Rdf(1.0,topoRdfExtent,to_string(i)+"-"+to_string(j))); //bin size 1 as incremental in discrete shells
        }
    }

    //construct rdfs
    int sizeIndexI, sizeIndexJ;
    vector<int> shell0, shell1, shell2; //prev, current, next shell
    for(int i=0; i<nNodes; ++i){//loop over all nodes
        shell0.clear();
        shell1.clear();
        shell1.push_back(i);
        sizeIndexI=nodes[i].sizeIndex;
        //add self terms for "shell 0"
        topoRdfShellSizes[sizeIndexI][0]=++topoRdfShellSizes[sizeIndexI][0];
        topoPartialRdfs[sizeIndexI*nRingSizes+sizeIndexI].addValue(0);
        shell2=nextTopologicalShell(shell1,shell0);
        for(int shell=1; shell<topoRdfExtent; ++shell){//loop over all shells (excluding incremented as will not be included in rdf)
            shell2=nextTopologicalShell(shell1,shell0);
            topoRdfShellSizes[sizeIndexI][shell]=topoRdfShellSizes[sizeIndexI][shell]+shell2.size(); //increase number of members in shell by size of shell
            for(int j=0; j<shell2.size(); ++j){//add to histogram
                sizeIndexJ=nodes[shell2[j]].sizeIndex;
                topoPartialRdfs[sizeIndexI*nRingSizes+sizeIndexJ].addValue(shell);
            }
            shell0=shell1;
            shell1=shell2;
        }
    }

    return;
}

void Network::analyseClusters() {
    //find clusters of each ring size and calculate sizes

    //find clusters of the same ring size and analyse
    vector< vector<int> > clusterSizes(nRingSizes);
    int n0,n1;
    vector<bool> checkNode(nNodes,true); //marks if node already part of a cluster
    for(int i=0, j=minRingSize; i<nRingSizes; ++i, ++j){//loop over ring sizes
        vector< vector<int> > clusters;
        vector<int> cluster;
        //find clusters
        for(int k=0; k<nNodes; ++k){//loop over nodes
            if(checkNode[k] && nodes[k].size==j){//select starting node of correct size
                cluster.clear();
                cluster.push_back(k);
                checkNode[k]=false;
                bool searchComplete=false;
                vector<int> prevSearch,search;
                prevSearch.clear();
                prevSearch.push_back(k);
                do {//find members of same cluster
                    search.clear();
                    for(int l=0; l<prevSearch.size(); ++l){
                        n0=prevSearch[l];
                        for(int m=0; m<nodes[n0].size; ++m){
                            n1=nodes[n0].connections[m];
                            if(checkNode[n1] && nodes[n1].size==j){
                                search.push_back(n1);
                                checkNode[n1]=false;
                            }
                        }
                    }
                    for(int m=0; m<search.size(); ++m) cluster.push_back(search[m]);
                    prevSearch=search;
                    if(search.size()==0) searchComplete=true;
                }while(!searchComplete);
                clusters.push_back(cluster);
            }
        }
        //get cluster sizes
        for(int k=0; k<clusters.size(); ++k) clusterSizes[i].push_back(clusters[k].size());
    }

    //calculate cluster size distributions
    clusterValues.resize(nRingSizes);
    clusterDistributions.resize(nRingSizes);
    for(int i=0;i<nRingSizes;++i){
        //get unique cluster size values
        vector<int> rawClusterSizes=clusterSizes[i];
        sort(rawClusterSizes.begin(),rawClusterSizes.end());
        vector<int> values;
        values.clear();
        if(rawClusterSizes.size()>0) {
            values.push_back(rawClusterSizes[0]);
            for (int j = 1; j < rawClusterSizes.size(); ++j) if (rawClusterSizes[j] != values.rbegin()[0]) values.push_back(rawClusterSizes[j]);
        }
        //get count of each cluster size
        map<int,int> refMap;
        refMap.clear();
        for(int j=0; j<values.size(); ++j) refMap[values[j]]=j;
        vector<int> dist(values.size(),0);
        for(int j=0; j<rawClusterSizes.size(); ++j) ++dist[refMap[rawClusterSizes[j]]];
        clusterValues[i]=values;
        clusterDistributions[i]=dist;
    }
}


void Network::analyseAssortativity() {
    //calculate assortative mixing Pearson correlation coefficient by two methods, and theoretically from alpha

    assortativeMixing.resize(3,0.0);
    //calculate useful moments
    double n, n1=0.0, n2=0.0, n3=0.0;
    int sizeI, sizeJ;
    for(int i=0; i<nRingSizes; ++i){
        sizeI=i+minRingSize;
        n=sizeI*ringStatistics[i];
        n1=n1+n;
        n2=n2+n*sizeI;
        n3=n3+n*sizeI*sizeI;
    }

    //method 1, using pi matrix
    double connectedCorrFunc=0.0;
    for(int i=0; i<nRingSizes; ++i){
        sizeI=i+minRingSize;
        for(int j=0; j<nRingSizes; ++j){
            sizeJ=j+minRingSize;
            connectedCorrFunc=connectedCorrFunc+(sizeI-1.0)*(sizeJ-1.0)*sizeI*ringStatistics[i]*(piMatrix[i][j]-sizeJ*ringStatistics[j]/6.0);
        }
    }
    assortativeMixing[0]=connectedCorrFunc/(n3-n2*n2/n1);

    //method 2, using edges
    int indexI, indexJ;
    double sum1=0.0, sum2=0.0, sum3=0.0, nEdges=0.0;
    for(int i=0; i<nNodes; ++i){
        indexI=i;
        for(int j=0; j<nodes[i].size; ++j){
            indexJ=nodes[i].connections[j];
            if(indexI<indexJ){//prevent double counting
                sizeI=nodes[indexI].size;
                sizeJ=nodes[indexJ].size;
                sum1=sum1+sizeI*sizeJ;
                sum2=sum2+sizeI+sizeJ;
                sum3=sum3+sizeI*sizeI+sizeJ*sizeJ;
                nEdges=nEdges+1.0;
            }
        }
    }
    sum1=sum1/nEdges;
    sum2=pow(0.5*sum2/nEdges,2);
    sum3=0.5*sum3/nEdges;
    assortativeMixing[1]=(sum1-sum2)/(sum3-sum2);

    //theoretical calculation from alpha and aboav-weaire law
    connectedCorrFunc=0.0;
    double mi;
    for(int i=0; i<nRingSizes; ++i){
        sizeI=i+minRingSize;
        mi=0.0;
        for(int j=0; j<nRingSizes; ++j){
            sizeJ=j+minRingSize;
            mi=mi+sizeJ*piMatrix[i][j];
        }
        mi=n2+n1*(1.0-aboavWeaireParams[0])*(sizeI-n1);
        mi=mi/sizeI;
        connectedCorrFunc=connectedCorrFunc+sizeI*sizeI*ringStatistics[i]*mi;
//        connectedCorrFunc=connectedCorrFunc+sizeI*ringStatistics[i]*(n2+n1*(1.0-aboavWeaireParams[0])*(sizeI-n1));
//        connectedCorrFunc=connectedCorrFunc+sizeI*ringStatistics[i]*n2;
//        connectedCorrFunc=connectedCorrFunc-sizeI*ringStatistics[i]*n1*n1;
//        connectedCorrFunc=connectedCorrFunc+sizeI*sizeI*ringStatistics[i]*n1;
//        connectedCorrFunc=connectedCorrFunc+sizeI*ringStatistics[i]*n1*aboavWeaireParams[0]*n1;
//        connectedCorrFunc=connectedCorrFunc-sizeI*sizeI*ringStatistics[i]*n1*aboavWeaireParams[0];
    }
//    connectedCorrFunc=2*n2*n1-n1*n1*n1+n1*n1*n1*aboavWeaireParams[0]-n1*n2*aboavWeaireParams[0];
//    connectedCorrFunc=2*n2*n1-n1*n1*n1+aboavWeaireParams[0]*n1*(n1*n1-n2);
    connectedCorrFunc=connectedCorrFunc-n2*n2/n1;
    assortativeMixing[2]=connectedCorrFunc/(n3-n2*n2/n1);

//    cout<<aboavWeaireParams[0]<<" "<<assortativeMixing[0]<<" "<<assortativeMixing[1]<<" "<<assortativeMixing[2]<<endl;
    return;
}

void Network::analyseRingAreas() {
    //calculate area of each ring in the network using shoelace method

    //loop over each ring and calculate dimensionless area
    ringAreas.resize(nRings);
    int ringSize;
    double area;
    Crd2d ringOrigin, p0, p1;
    for(int i=0; i<vertexRings.size(); ++i){
        //recentre coordinates on first point with pbcs
        ringSize=vertexRings[i].size;
        ringOrigin=nodes[vertexRings[i].id].coordinate;
        vector<Crd2d> ringCrds(ringSize+1);
        for(int j=0; j<ringSize+1;++j){
            ringCrds[j]=recentreCrdByCrd(ringOrigin,vertices[vertexRings[i].chain[j]].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
        }
        //calculate area
        area=0.0;
        for(int j=0; j<ringSize; ++j){
            p0=ringCrds[j];
            p1=ringCrds[j+1];
            area+=p0.y*p1.x-p0.x*p1.y;
//            cout<<p0.x<<" "<<p0.y<<", ";
        }
        area*=0.5;
        ringAreas[i]=fabs(area);
    }
}

void Network::analyseAtomicGeometry() {
    //calculate mean and standard deviation of bond lengths and angles in different ring sizes
    //ONLY ACCURATE FOR PERIODIC AS DOUBLE COUNTS NON-EDGE BOND LENGTHS

    atomicGeomBondMean.resize(nRingSizes+1);
    atomicGeomBondStd.resize(nRingSizes+1);
    atomicGeomAngleMean.resize(nRingSizes+1);
    atomicGeomAngleStd.resize(nRingSizes+1);
    int ringSize;
    vector<double> ringBonds, ringAngles;
    vector<Crd2d> ringCrds;
    Crd2d ringCentre;
    Vec2d v0, v1;
    double d, mean, stdev;
    atomicBondDistribution.clear();
    atomicAngleDistribution.clear();
    //loop over all ring sizes
    for(int i=0, j=minRingSize; i<nRingSizes; ++i, ++j){
        //reset variables
        ringBonds.clear();
        ringAngles.clear();
        //loop over all rings and find those which have search ring size
        for(int k=0; k<vertexRings.size(); ++k){
            ringSize=vertexRings[k].size;
            if(ringSize==j){//ring size matches
                ringCrds.resize(ringSize+2); //end repeated at start and vice versa
                ringCentre=nodes[vertexRings[k].id].coordinate;
                ringCrds[0]=recentreCrdByCrd(ringCentre,vertices[vertexRings[k].chain[ringSize-1]].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
                for(int l=0; l<ringSize+1; ++l){
                    ringCrds[l+1]=recentreCrdByCrd(ringCentre,vertices[vertexRings[k].chain[l]].coordinate,periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
                }
                //calculate bond lengths
                for(int l=1; l<=ringSize; ++l){
                    v0=Vec2d(ringCrds[l],ringCrds[l+1],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
                    d=v0.length();
                    ringBonds.push_back(d);
                    atomicBondDistribution.push_back(d);
                }
                //calculate bond angles
                for(int l=1; l<=ringSize; ++l){
                    v0=Vec2d(ringCrds[l],ringCrds[l+1],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
                    v1=Vec2d(ringCrds[l],ringCrds[l-1],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
                    v0.normalise();
                    v1.normalise();
                    d=vectorDotProduct(v0,v1);
                    d=acos(d)*180.0/M_PI;
                    ringAngles.push_back(d);
                    atomicAngleDistribution.push_back(d);
                }
            }
        }
        meanAndStdDeviation(ringBonds,mean,stdev);
        atomicGeomBondMean[i]=mean;
        atomicGeomBondStd[i]=stdev;
        meanAndStdDeviation(ringAngles,mean,stdev);
        atomicGeomAngleMean[i]=mean;
        atomicGeomAngleStd[i]=stdev;
    }
    meanAndStdDeviation(atomicBondDistribution,mean,stdev);
    atomicGeomBondMean[nRingSizes]=mean;
    atomicGeomBondStd[nRingSizes]=stdev;
    meanAndStdDeviation(atomicAngleDistribution,mean,stdev);
    atomicGeomAngleMean[nRingSizes]=mean;
    atomicGeomAngleStd[nRingSizes]=stdev;
}

void Network::meanAndStdDeviation(vector<double> &values, double &mean, double &stdev) {
    //calculate mean and standard deviation of provied values
    mean=0.0;
    stdev=0.0;
    int n=values.size();
    if(n>0) {
        for (int i = 0; i < n; ++i) {
            mean += values[i];
            stdev += values[i] * values[i];
        }
        mean /= n;
        stdev = sqrt(stdev / n - mean * mean);
    }
}

void Network::geometryOptimiseAtomicNetworkPeriodic() {
    //optimise atomic network with keating potential

    //make list of coordinates
    vector<Crd2d> vertexCoordinates(nVertices);
    for(int i=0; i<nVertices; ++i) vertexCoordinates[i]=vertices[i].coordinate;

    //make list of bonds
    vector<Pair> vertexBonds;
    Pair bond;
    for(int i=0; i<nVertices; ++i){//all unique pairs
        bond.a=i;
        for(int j=0; j<vertices[i].coordination; ++j){
            bond.b=vertices[i].connections[j];
            if(bond.a<bond.b) vertexBonds.push_back(bond); //prevent double counting
        }
    }

    //make list of angles
    vector<Trio> vertexAngles;
    Trio angle;
    for(int i=0; i<nVertices; ++i){//all unique angles about each central atom
        angle.b=i; //central atom as b
        angle.a=vertices[i].connections[0];
        angle.c=vertices[i].connections[1];
        vertexAngles.push_back(angle);
        angle.a=vertices[i].connections[1];
        angle.c=vertices[i].connections[2];
        vertexAngles.push_back(angle);
        angle.a=vertices[i].connections[0];
        angle.c=vertices[i].connections[2];
        vertexAngles.push_back(angle);
    }

    //make list of line intersections - each edge of ring with other edges of ring
    vector<DoublePair> edges;
    vector<Pair> ringEdges;
    edges.clear();
    for(int i=0; i<vertexRings.size(); ++i){//loop over rings
        ringEdges.clear();
        //make list of edges
        for(int j=0; j<vertexRings[i].size; ++j) ringEdges.push_back(Pair(vertexRings[i].chain[j],vertexRings[i].chain[j+1]));
        //prevent overlap of first edge with rest except neighbouring edges
        for(int j=2; j<vertexRings[i].size-1; ++j){
            edges.push_back(DoublePair(ringEdges[0].a,ringEdges[0].b,ringEdges[j].a,ringEdges[j].b));
        }
        //prevent overlap of remaining edges with rest
        for(int j=1; j<ringEdges.size()-1; ++j){
            for(int k=j+2; k<ringEdges.size(); ++k){
                edges.push_back(DoublePair(ringEdges[j].a,ringEdges[j].b,ringEdges[k].a,ringEdges[k].b));
            }
        }
    }

    //minimise
    vector<int> emptyFixed; //placeholder
    emptyFixed.clear();
    KeatingPeriodicGO optimiser;
    optimiser.setCoordinates(vertexCoordinates);
    optimiser.setSystemParameters(vertexBonds,vertexAngles,emptyFixed,edges);
    optimiser.setPotentialParameters(keatingA,keatingAlpha,keatingBeta);
    optimiser.setOptisationParameters(atomicGeomOptCC,atomicLineSearchInc,atomicGeomOptMaxIt,false);
    optimiser.setPeriodicBoundary(periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
    atomicGeomOptStatus=optimiser.steepestDescent();
    vertexCoordinates=optimiser.getMinimisedCoordinates();

    //check entire atomic network for overlaps
    vector<Pair> checkEdges;
    checkEdges.clear();
    for(int i=0; i<vertexRings.size(); ++i){//make list of edges, not unique but ok as will double check
        for(int j=0; j<vertexRings[i].size; ++j) checkEdges.push_back(Pair(vertexRings[i].chain[j],vertexRings[i].chain[j+1]));
    }
    bool intersection=false;
    Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
    for(int i=0; i<checkEdges.size()-1; ++i){//check all edges against each other
        for(int j=i+1; j<checkEdges.size(); ++j){
            line1a=vertexCoordinates[checkEdges[i].a];
            line1b=minimumImageCrd(line1a,vertexCoordinates[checkEdges[i].b],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            line2a=minimumImageCrd(line1a,vertexCoordinates[checkEdges[j].a],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            line2b=minimumImageCrd(line2a,vertexCoordinates[checkEdges[j].b],periodicBoxX,periodicBoxY,rPeriodicBoxX,rPeriodicBoxY);
            intersection=properIntersectionLines(line1a,line1b,line2a,line2b);
            if (intersection) break;
        }
    }
    if(intersection) atomicGeomOptStatus=0;

    //if no intersections get energy and iteration information and update coordinates
    if(atomicGeomOptStatus==1){
        atomicGeomOptEnergy=optimiser.getEnergy();
        atomicGeomOptIterations=optimiser.getIterations();
        for(int i=0; i<nVertices; ++i) vertices[i].coordinate=vertexCoordinates[i];
    }

    return;
}

void Network::geometryOptimiseAtomicNetworkAperiodic() {
    //optimise atomic network with keating potential

    //make list of coordinates
    vector<Crd2d> vertexCoordinates(nVertices);
    for(int i=0; i<nVertices; ++i) vertexCoordinates[i]=vertices[i].coordinate;

    //make list of bonds
    vector<Pair> vertexBonds;
    Pair bond;
    for(int i=0; i<nVertices; ++i){//all unique pairs
        bond.a=i;
        for(int j=0; j<vertices[i].coordination; ++j){
            bond.b=vertices[i].connections[j];
            if(bond.a<bond.b) vertexBonds.push_back(bond); //prevent double counting
        }
    }

    //make list of angles
    vector<Trio> vertexAngles;
    Trio angle;
    for(int i=0; i<nVertices; ++i){//all unique angles about each central atom
        angle.b=i; //central atom as b
        //first angle for both two and three coordinate vertices
        angle.a=vertices[i].connections[0];
        angle.c=vertices[i].connections[1];
        vertexAngles.push_back(angle);
        //if three coordinate add two more angles
        if(vertices[i].coordination==3){
            angle.a=vertices[i].connections[1];
            angle.c=vertices[i].connections[2];
            vertexAngles.push_back(angle);
            angle.a=vertices[i].connections[0];
            angle.c=vertices[i].connections[2];
            vertexAngles.push_back(angle);
        }
    }

    //make list of line intersections - each edge of ring with other edges of ring
    vector<DoublePair> edges;
    vector<Pair> ringEdges;
    edges.clear();
    for(int i=0; i<vertexRings.size(); ++i){//loop over rings
        ringEdges.clear();
        //make list of edges
        for(int j=0; j<vertexRings[i].size; ++j) ringEdges.push_back(Pair(vertexRings[i].chain[j],vertexRings[i].chain[j+1]));
        //prevent overlap of first edge with rest except neighbouring edges
        for(int j=2; j<vertexRings[i].size-1; ++j){
            edges.push_back(DoublePair(ringEdges[0].a,ringEdges[0].b,ringEdges[j].a,ringEdges[j].b));
        }
        //prevent overlap of remaining edges with rest
        for(int j=1; j<ringEdges.size()-1; ++j){
            for(int k=j+2; k<ringEdges.size(); ++k){
                edges.push_back(DoublePair(ringEdges[j].a,ringEdges[j].b,ringEdges[k].a,ringEdges[k].b));
            }
        }
    }

    //minimise
    vector<int> emptyFixed; //placeholder
    emptyFixed.clear();
    KeatingAperiodicGO optimiser;
    optimiser.setCoordinates(vertexCoordinates);
    optimiser.setSystemParameters(vertexBonds,vertexAngles,emptyFixed,edges);
    optimiser.setPotentialParameters(keatingA,keatingAlpha,keatingBeta);
    optimiser.setOptisationParameters(atomicGeomOptCC,atomicLineSearchInc,atomicGeomOptMaxIt,false);
    atomicGeomOptStatus=optimiser.steepestDescent();
    vertexCoordinates=optimiser.getMinimisedCoordinates();

    //check entire atomic network for overlaps
    vector<Pair> checkEdges;
    checkEdges.clear();
    for(int i=0; i<vertexRings.size(); ++i){//make list of edges, not unique but ok as will double check
        for(int j=0; j<vertexRings[i].size; ++j) checkEdges.push_back(Pair(vertexRings[i].chain[j],vertexRings[i].chain[j+1]));
    }
    bool intersection=false;
    for(int i=0; i<checkEdges.size()-1; ++i){//check all edges against each other
        for(int j=i+1; j<checkEdges.size(); ++j){
            intersection=properIntersectionLines(vertexCoordinates[checkEdges[i].a],vertexCoordinates[checkEdges[i].b],vertexCoordinates[checkEdges[j].a],vertexCoordinates[checkEdges[j].b]);
            if (intersection) break;
        }
    }
    if(intersection) atomicGeomOptStatus=0;

    //if no intersections get energy and iteration information and update coordinates
    if(atomicGeomOptStatus==1){
        atomicGeomOptEnergy=optimiser.getEnergy();
        atomicGeomOptIterations=optimiser.getIterations();
        for(int i=0; i<nVertices; ++i) vertices[i].coordinate=vertexCoordinates[i];
    }

    return;
}

//###### Write out ######

void Network::write() {
    //write out to files

    writeDual();
    if(periodicVisualisation) writePeriodicDualNetwork();
    if(globalGeomOpt) writeGeometryOptimisationEnergy();
    writeRingStatistics();
    writeAboavWeaire();
    if(periodic && spatialRdf) writeSpatialPartialRdfs();
    if(periodic && topoRdf) writeTopoPartialRdfs();
    if(periodic && clustering) writeClusters();
    if(assortativeMix) writeAssortativeMixing();
    if(convertToAtomic){
        writeAtomicNetwork();
        if(atomicGeomOpt) writeAtomicGeometryOptimisation();
        writeAtomicGeometrySummary();
        if(areaLaw) writeRingAreas();
    }

    return;
}

void Network::writeDual() {
    //write out dual graph information

    string crdOutputFileName=outPrefix+"dual_coordinates.out";
    ofstream crdOutputFile(crdOutputFileName, ios::in|ios::trunc);
    crdOutputFile<<fixed<<showpoint<<setprecision(6);
    for(int i=0; i<nNodes; ++i) writeFileCrd(crdOutputFile,nodes[i].coordinate);
    crdOutputFile.close();

    string cnxOutputFileName=outPrefix+"dual_connectivity.out";
    ofstream cnxOutputFile(cnxOutputFileName, ios::in|ios::trunc);
    vector<int> cnxs;
    for(int i=0; i<nNodes; ++i){
        cnxs=nodes[i].connections;
        cnxs.insert(cnxs.begin(),i); //add node id to front of vector
        writeFileRowVector(cnxOutputFile,cnxs);
    }
    cnxOutputFile.close();

    string sizeOutputFileName=outPrefix+"dual_size.out";
    ofstream sizeOutputFile(sizeOutputFileName, ios::in|ios::trunc);
    vector<int> line(2); //size and whether an edge
    for(int i=0; i<nNodes; ++i){
        line[0]=nodes[i].size;
        line[1]=nodes[i].edge;
        writeFileRowVector(sizeOutputFile,line);
    }
    sizeOutputFile.close();

    return;
}

void Network::writeAtomicNetwork() {
    //write out atomic coordinates, connectivities, rings and ring sizes

    string crdOutputFileName=outPrefix+"graph_coordinates.out";
    ofstream crdOutputFile(crdOutputFileName, ios::in|ios::trunc);
    crdOutputFile<<fixed<<showpoint<<setprecision(6);
    for(int i=0; i<nVertices; ++i) writeFileCrd(crdOutputFile,vertices[i].coordinate);
    crdOutputFile.close();

    string cnxOutputFileName=outPrefix+"graph_connectivity.out";
    ofstream cnxOutputFile(cnxOutputFileName, ios::in|ios::trunc);
    vector<int> cnxs;
    for(int i=0; i<nVertices; ++i){
        cnxs.clear();
        cnxs.push_back(i);
        for(int j=0; j<vertices[i].coordination; ++j) cnxs.push_back(vertices[i].connections[j]);
        writeFileRowVector(cnxOutputFile,cnxs);
    }
    cnxOutputFile.close();

    string sizeOutputFileName=outPrefix+"graph_size.out";
    ofstream sizeOutputFile(sizeOutputFileName, ios::in|ios::trunc);
    for(int i=0; i<vertexRings.size(); ++i) writeFileLine(sizeOutputFile,vertexRings[i].size);
    sizeOutputFile.close();

    string ringOutputFileName=outPrefix+"graph_rings.out";
    ofstream ringOutputFile(ringOutputFileName, ios::in|ios::trunc);
    for(int i=0; i<vertexRings.size(); ++i) writeFileRowVector(ringOutputFile,vertexRings[i].path());
    ringOutputFile.close();

    return;
}

void Network::writePeriodicDualNetwork() {
    //only for visualisation, calculate periodic images of network for visualisation

    //set limits of visualisation region
    double imageProportion=1.0;
    double leftLimit=-imageProportion*periodicBoxX, rightLimit=(1.0+imageProportion)*periodicBoxX;
    double bottomLimit=-imageProportion*periodicBoxY, topLimit=(1.0+imageProportion)*periodicBoxY;

    //find nodes which would periodically fit in vis region, make new nodes
    Crd2d imageCrd; //coordinate of image
    int nPeriodicNodes=nNodes;
    vector<Node> periodicNodes=nodes;
    map <int, vector<int> > nodeImageMap; //node references for all images of node
    for(int i=0; i<nNodes; ++i) nodeImageMap[i]=vector<int> (1,i); //add image in original box
    for(int n=0; n<nNodes; ++n){//loop over all nodes
        for(int i=-1; i<=1; ++i){//loop over nearest images
            for(int j=-1; j<=1; ++j){
                if(abs(i)+abs(j)!=0){
                    imageCrd=periodicNodes[n].coordinate;
                    imageCrd.x=imageCrd.x+i*periodicBoxX;
                    imageCrd.y=imageCrd.y+j*periodicBoxY;
                    if(imageCrd.x>leftLimit && imageCrd.x<rightLimit && imageCrd.y>bottomLimit && imageCrd.y<topLimit){//if in region make new coordinate
                        periodicNodes.push_back(periodicNodes[n]);
                        periodicNodes[nPeriodicNodes].coordinate=imageCrd;
                        nodeImageMap[n].push_back(nPeriodicNodes);
                        nPeriodicNodes=++nPeriodicNodes;
                    }
                }
            }
        }
    }

    //remake connections to nearest image
    int nodeRef;
    Vec2d imageVec;
    vector<int> images, cnxCopy;
    vector<double> distances;
    double micX=periodicBoxX*0.5, micY=periodicBoxY*0.5; //minimum image convention
    for(int i=0; i<nPeriodicNodes; ++i){
        cnxCopy=periodicNodes[i].connections;
        periodicNodes[i].connections.clear();
        for(int j=0; j<periodicNodes[i].size; ++j){
            nodeRef=cnxCopy[j];
            images=nodeImageMap[nodeRef];
            for(int k=0; k<images.size(); ++k){
                imageVec=Vec2d(periodicNodes[i].coordinate,periodicNodes[images[k]].coordinate);
                if(fabs(imageVec.x)<=micX && fabs(imageVec.y)<=micY){
                    periodicNodes[i].addConnection(images[k]);
                    break;
                }
            }
        }
    }
    //set undercoordinated nodes as edges
    for(int i=0; i<nPeriodicNodes; ++i){
        if(periodicNodes[i].size!=periodicNodes[i].connections.size()) periodicNodes[i].edge=true;
    }

    //write dual coordinates and connectivities
    string crdOutputFileName=outPrefix+"dual_periodic_coordinates.out";
    ofstream crdOutputFile(crdOutputFileName, ios::in|ios::trunc);
    crdOutputFile<<fixed<<showpoint<<setprecision(6);
    for(int i=0; i<nPeriodicNodes; ++i) writeFileCrd(crdOutputFile,periodicNodes[i].coordinate);
    crdOutputFile.close();

    string sizeOutputFileName=outPrefix+"dual_periodic_size.out";
    ofstream sizeOutputFile(sizeOutputFileName, ios::in|ios::trunc);
    for(int i=0; i<nPeriodicNodes; ++i) writeFileLine(sizeOutputFile,periodicNodes[i].size,int(periodicNodes[i].edge));
    sizeOutputFile.close();

    string cnxOutputFileName=outPrefix+"dual_periodic_connectivity.out";
    ofstream cnxOutputFile(cnxOutputFileName, ios::in|ios::trunc);
    vector<int> cnxs;
    for(int i=0; i<nPeriodicNodes; ++i){
        cnxs=periodicNodes[i].connections;
        cnxs.insert(cnxs.begin(),i); //add node id to front of vector
        writeFileRowVector(cnxOutputFile,cnxs);
    }
    cnxOutputFile.close();

    string latticeOutputFileName=outPrefix+"periodic_lattice_dim.out";
    ofstream latticeOutputFile(latticeOutputFileName, ios::in|ios::trunc);
    writeFileValue(latticeOutputFile,periodicBoxX);
    writeFileValue(latticeOutputFile,periodicBoxY);
    latticeOutputFile.close();

    if(convertToAtomic) writePeriodicAtomicNetwork(periodicNodes);

    return;
}

void Network::writePeriodicAtomicNetwork(vector<Node> &periodicNodes) {
    //only for visualisation, calculate periodic images of network for visualisation

    int nPeriodicNodes=periodicNodes.size();

    //find periodic node rings, created from find node ring(s) functions
    vector<Ring> periodicNodeRings; //periodic rings
    for(int i=0; i<nPeriodicNodes; ++i){//find node rings
        if(!periodicNodes[i].edge){

            //find single node ring
            int nSize=periodicNodes[i].size;
            vector<int> ring(nSize); //ring around central node
            vector<int> centralCnxs=periodicNodes[i].connections;
            vector<int> outerCnxs;
            ring[0]=centralCnxs[0]; //pick first connection as starting point
            outerCnxs=periodicNodes[ring[0]].connections;
            for(int j=0; j<outerCnxs.size(); ++j){//find second connection in ring
                if(checkVectorForValue(centralCnxs,outerCnxs[j])){
                    ring[1]=outerCnxs[j];
                    break;
                };
            }
            for(int j=2; j<nSize; ++j){//find remaining connections
                outerCnxs=periodicNodes[ring[j-1]].connections;
                removeValueFromVectorByRef(outerCnxs,ring[j-2]); //remove previous node in chain
                for(int k=0; k<outerCnxs.size(); ++k){
                    if(checkVectorForValue(centralCnxs,outerCnxs[k])) {
                        ring[j] = outerCnxs[k];
                        break;
                    }
                }
            }

            //add to vector
            periodicNodeRings.push_back(Ring(periodicNodes[i].size,ring,i));
        }
    }

    //make unique vertex list and vertex rings from node rings, taken from convert dual to atomic network
    int nPeriodicRings=periodicNodeRings.size();
    int nPeriodicVertices=0;
    vector<Vertex> periodicVertices;
    vector<Ring> periodicVertexRings;
    periodicVertices.clear();
    periodicVertexRings.clear();
    Triangle nodeTriangle;
    vector<int> ringPath;
    string triangleID;
    map<string,int> triangleVertexMap; //triangle id to vertex index map
    for(int i=0; i<nPeriodicRings; ++i){//loop over node rings
        ringPath.clear();
        for(int j=0; j<periodicNodeRings[i].size; ++j){//loop over all node triangles (=vertices) in ring
            nodeTriangle=Triangle(periodicNodeRings[i].id,periodicNodeRings[i].chain[j],periodicNodeRings[i].chain[j+1]);
            nodeTriangle.sort();
            triangleID=nodeTriangle.getID();
            if(triangleVertexMap.count(triangleID)==0){//if vertex not already generated, make vertex
                periodicVertices.push_back(Vertex(nodeTriangle.a, nodeTriangle.b, nodeTriangle.c));
                triangleVertexMap[triangleID]=nPeriodicVertices;
                nPeriodicVertices=++nPeriodicVertices;
            }
            ringPath.push_back(triangleVertexMap[triangleID]);
        }
        periodicVertexRings.push_back(Ring(ringPath.size(),ringPath,i)); //make vertex ring
    }

    //make vertex connections
    for(int i=0; i<periodicVertexRings.size(); ++i){//triply degenerate but vertex only accepts addition if unique
        for(int j=0; j<periodicVertexRings[i].size; ++j){
            periodicVertices[periodicVertexRings[i].chain[j]].addConnection(periodicVertexRings[i].chain[j+1]);
            periodicVertices[periodicVertexRings[i].chain[j+1]].addConnection(periodicVertexRings[i].chain[j]);
        }
    }

    //make vertex coordinates
    vector<Crd2d> nodeTriCrds(3);
    for(int i=0; i<nPeriodicVertices; ++i){
        for(int j=0; j<3; ++j) nodeTriCrds[j]=periodicNodes[periodicVertices[i].dualNodes[j]].coordinate;
        periodicVertices[i].coordinate=crdCentreOfMass(nodeTriCrds);
    }

    //generate list of whether original or image
    vector<int> imageList(nPeriodicRings,0);
    for(int i=0; i<nRings; ++i) imageList[i]=1;

    //write atomic coordinates
    string crdOutputFileName=outPrefix+"graph_periodic_coordinates.out";
    ofstream crdOutputFile(crdOutputFileName, ios::in|ios::trunc);
    crdOutputFile<<fixed<<showpoint<<setprecision(6);
    for(int i=0; i<nPeriodicVertices; ++i) writeFileCrd(crdOutputFile,periodicVertices[i].coordinate);
    crdOutputFile.close();

    string sizeOutputFileName=outPrefix+"graph_periodic_size.out";
    ofstream sizeOutputFile(sizeOutputFileName, ios::in|ios::trunc);
    for(int i=0; i<periodicVertexRings.size(); ++i){
        writeFileLine(sizeOutputFile,periodicVertexRings[i].size,false);
        writeFileLine(sizeOutputFile,"  ",false);
        writeFileLine(sizeOutputFile,imageList[i]);
    }
    sizeOutputFile.close();

    string ringOutputFileName=outPrefix+"graph_periodic_rings.out";
    ofstream ringOutputFile(ringOutputFileName, ios::in|ios::trunc);
    for(int i=0; i<periodicVertexRings.size(); ++i) writeFileRowVector(ringOutputFile,periodicVertexRings[i].path());
    ringOutputFile.close();

    return;
}

void Network::writeGeometryOptimisationEnergy() {
    //write out final energy of geometry optimisation

    string geomOutputFileName=outPrefix+"analysis_geometry_opt.out";
    ofstream geomOutputFile(geomOutputFileName, ios::in|ios::trunc);
    geomOutputFile<<fixed<<showpoint<<setprecision(6);
    writeFileValue(geomOutputFile, geomOptEnergy);
    writeFileValue(geomOutputFile, geomOptIterations);
    geomOutputFile.close();

    return;
}

void Network::writeRingStatistics() {
    //write ring statistics and pi matrix

    string rsOutputFileName=outPrefix+"analysis_ring_statistics.out";
    ofstream rsOutputFile(rsOutputFileName, ios::in|ios::trunc);
    rsOutputFile<<fixed<<showpoint<<setprecision(6);
    writeFileRowVector(rsOutputFile,ringStatistics);
    rsOutputFile.close();

    string piOutputFileName=outPrefix+"analysis_pi_matrix.out";
    ofstream piOutputFile(piOutputFileName, ios::in|ios::trunc);
    piOutputFile<<fixed<<showpoint<<setprecision(6);
    writeFileMatrix(piOutputFile,piMatrix);
    piOutputFile.close();

    return;
}

void Network::writeAboavWeaire() {
    //write alpha, mu and rsq
    string awOutputFileName=outPrefix+"analysis_aw.out";
    ofstream awOutputFile(awOutputFileName, ios::in|ios::trunc);
    awOutputFile<<fixed<<showpoint<<setprecision(6);
    for(int i=0; i<3; ++i) writeFileValue(awOutputFile,aboavWeaireParams[i]);
    awOutputFile.close();
    return;
}

void Network::writeSpatialPartialRdfs() {
    //write out partial rdfs to single file

    string rdfOutputFileName=outPrefix+"analysis_spatial_rdfs.out";
    ofstream rdfOutputFile(rdfOutputFileName, ios::in|ios::trunc);
    rdfOutputFile<<setprecision(8);
    writeFileLine(rdfOutputFile,"Partial RDF Analysis, bin width, nBins, densities, N and ideal distances");
    writeFileValue(rdfOutputFile,spatialRdfBinwidth);
    writeFileLine(rdfOutputFile,int(floor(spatialRdfExtent/spatialRdfBinwidth)));
    writeFileRowVector(rdfOutputFile,spatialRdfDensities);
    writeFileRowVector(rdfOutputFile,pVector);
    writeFileMatrix(rdfOutputFile,harmonicR0Matrix);
    for(int i=0; i<nRingSizes*nRingSizes; ++i){
        writeFileLine(rdfOutputFile,"Non-normalised partial rdf "+spatialPartialRdfs[i].id);
        writeFileColVector(rdfOutputFile,spatialPartialRdfs[i].vectorHistogram());
    }
    rdfOutputFile.close();
    return;
}

void Network::writeTopoPartialRdfs() {
    //write out partial rdfs to single file

    string rdfOutputFileName=outPrefix+"analysis_topo_rdfs.out";
    ofstream rdfOutputFile(rdfOutputFileName, ios::in|ios::trunc);
    rdfOutputFile<<setprecision(8);
    for(int i=0; i<nRingSizes; ++i){
            writeFileLine(rdfOutputFile,"Topological shell sizes for central node size: "+to_string(minRingSize+i));
            writeFileColVector(rdfOutputFile,topoRdfShellSizes[i]);
    }
    for(int i=0; i<nRingSizes*nRingSizes; ++i){
        writeFileLine(rdfOutputFile,"Non-normalised partial topo rdf "+topoPartialRdfs[i].id);
        writeFileColVector(rdfOutputFile,topoPartialRdfs[i].vectorHistogram());
    }
    rdfOutputFile.close();
    return;
}

void Network::writeClusters() {
    //write out cluster size distribution to single file

    string clusterOutputFileName=outPrefix+"analysis_clusters.out";
    ofstream clusterOutputFile(clusterOutputFileName, ios::in|ios::trunc);
    
    //check clusters include all nodes
    int sum=0;
    for(int i=0; i<nRingSizes; ++i){
        for(int j=0; j<clusterValues[i].size(); ++j) sum+=clusterValues[i][j]*clusterDistributions[i][j];
    }
    if(sum!=nNodes) writeFileLine(clusterOutputFile,"ERROR");

    //write cluster distribution
    for(int i=0, j=minRingSize; i<nRingSizes; ++i,++j){
        writeFileLine(clusterOutputFile,"Cluster size distributions for rings of size: "+to_string(j));
        writeFileLine(clusterOutputFile,"Number of entries: "+to_string(clusterValues[i].size()));
        for(int j=0; j<clusterValues[i].size(); ++j){
            writeFileLine(clusterOutputFile,clusterValues[i][j],clusterDistributions[i][j]);
        }
    }
    clusterOutputFile.close();
    return;
}

void Network::writeAssortativeMixing() {
    //write pearson correlation coefficient calculated in three different ways
    string amOutputFileName=outPrefix+"analysis_assortative_mix.out";
    ofstream amOutputFile(amOutputFileName, ios::in|ios::trunc);
    amOutputFile<<fixed<<showpoint<<setprecision(6);
    writeFileRowVector(amOutputFile,assortativeMixing);
    amOutputFile.close();
    return;
}

void Network::writeAtomicGeometryOptimisation() {
    //write out success and final energy/iterations of atomic geometry optimisation
    string geomOutputFileName=outPrefix+"analysis_atomic_geometry_opt.out";
    ofstream geomOutputFile(geomOutputFileName, ios::in|ios::trunc);
    geomOutputFile<<fixed<<showpoint<<setprecision(6);
    writeFileValue(geomOutputFile, atomicGeomOptStatus);
    if(atomicGeomOptStatus==1){
        writeFileValue(geomOutputFile, atomicGeomOptEnergy);
        writeFileValue(geomOutputFile, atomicGeomOptIterations);
    }
    geomOutputFile.close();
    return;
}

void Network::writeAtomicGeometrySummary() {
    //write out average bond length and angle and entire distribution
    string geomOutputFileName=outPrefix+"analysis_atomic_geometry_sum.out";
    ofstream geomOutputFile(geomOutputFileName, ios::in|ios::trunc);
    geomOutputFile<<fixed<<showpoint<<setprecision(6);
    writeFileRowVector(geomOutputFile, atomicGeomBondMean);
    writeFileRowVector(geomOutputFile, atomicGeomBondStd);
    writeFileRowVector(geomOutputFile, atomicGeomAngleMean);
    writeFileRowVector(geomOutputFile, atomicGeomAngleStd);
    geomOutputFile.close();
    geomOutputFileName=outPrefix+"analysis_atomic_geometry_full.out";
    ofstream geomFullOutputFile(geomOutputFileName, ios::in|ios::trunc);
    geomFullOutputFile<<fixed<<showpoint<<setprecision(6);
    int index=0;
    for(int i=0, ringSize=minRingSize; i<nRingSizes; ++i, ++ringSize){
        int n=ringStatistics[i]*nRings*ringSize;
        for(int j=0; j<n; ++j){
            geomFullOutputFile<<ringSize<<"  "<<atomicBondDistribution[index]<<"  "<<atomicAngleDistribution[index]<<endl;
            ++index;
        }
    }
    geomFullOutputFile.close();
    return;
}

void Network::writeRingAreas(){
    //write absolute area of rings
    //write dimensionless area averaged over: ideal bond, all bonds, bonds of rings of same size (if have information)
    string areaOutputFileName=outPrefix+"analysis_ring_area.out";
    ofstream areaOutputFile(areaOutputFileName, ios::in|ios::trunc);
    areaOutputFile<<fixed<<showpoint<<setprecision(6);
    int n,m;
    double d0, d1, d2, d3;
    for(int i=0; i<nRings; ++i){
        n=vertexRings[i].size;
        m=n-minRingSize;
        d0=ringAreas[i];
        d1=d0/(keatingA*keatingA);
        d2=d0/(atomicGeomBondMean[nRingSizes]*atomicGeomBondMean[nRingSizes]);
        d3=d0/(atomicGeomBondMean[m]*atomicGeomBondMean[m]);
        areaOutputFile<<n<<"  "<<d0<<"  "<<d1<<"  "<<d2<<"  "<<d3<<endl;
    }
    areaOutputFile.close();
    return;
}
