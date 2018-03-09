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

void Network::setPotential(double sep) {
    //set potential model parameters
    atomicSeparation=sep;
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
    initialisePotentialModel();
    if(periodic) initialisePeriodicLattice();
    else initialiseAperiodicLattice();
    checkFidelity();
    initialiseMonteCarlo();

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
            nodeIndex=layer*yNodes+n;
            nodes[nodeIndex].coordinate.y=layer*yOffset;
            nodes[nodeIndex].coordinate.x=xCrd;
            xCrd=xCrd+r66;
        }
        stagger=!stagger;
    }
    //make aperiodic first


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
    rTargetMu=0.0;
    for(int i=0; i<nRingSizes;++i) rTargetMu=rTargetMu+(minRingSize+i)*(minRingSize+i)*targetPVector[i];
    rTargetMu=rTargetMu-36.0;
    rTargetMu=1.0/rTargetMu;

    //random generators
    metropolisGenerator.seed(mcSeed);
    nodePickGenerator.seed(mcSeed+1);
    metropolisDistribution=uniform_real_distribution<double>(0.0,1.0);
    nodePickDistribution=uniform_int_distribution<int>(0,nNodes-1); //-1 as inclusive

    return;
}

//###### Construction monte carlo ######

double Network::metropolisRandomNum() {
    //random between 0->1
    return metropolisDistribution(metropolisGenerator);
}

int Network::pickRandomNode() {
    //random between 0->nNodes-1
    return nodePickDistribution(nodePickGenerator);
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
    cout<<"consistent: "<<selfConsistent<<endl;

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
    cout<<"p vector: "<<pVectorPass<<" p matrix: "<<pMatrixPass<<endl;

    return;
}

//###### Write out ######

void Network::write() {
    //write out to files

    writeDual();

    return;
}

void Network::writeDual() {
    //write out dual graph information

    string crdOutputFileName=outPrefix+"dual_coordinates.out";
    ofstream crdOutputFile(crdOutputFileName, ios::in|ios::trunc);
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