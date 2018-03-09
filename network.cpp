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

void Network::setAnalysis(bool perVis) {
    //specify analysis to perform
    if(periodic) periodicVisualisation=perVis;
    else periodicVisualisation=false;
    return;
}

//###### GETTERS ######

//###### Construction Main ######

void Network::construct(ofstream &logfile) {
    //attempt to build network with the specified properties

    name="network "+to_string(mcSeed);
    writeFileLine(logfile,"constructing "+name);

    initialiseNetworkProperties();
    initialisePotentialModel();
    if(periodic) initialisePeriodicLattice();
    else initialiseAperiodicLattice();
    initialiseMonteCarlo();
    monteCarlo();
    checkFidelity();

    writeFileLine(logfile,name+" constructed");
    if(consistent) writeFileLine(logfile,name+" checked for consistency and passed");
    else writeFileLine(logfile,name+" failed consistency test");
    if(mcTargetReached) writeFileLine(logfile, name+" targets met in "+to_string(mcProposedMoves)+" monte carlo moves");
    else writeFileLine(logfile, name+" targets not met");




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

//###### Construction monte carlo ######

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
            acceptTrialMove=evaluateMetropolisCondition(trialMcEnergy,mcEnergy);
            if(acceptTrialMove) mcTargetReached=acceptDualSwitch(switchTriangles,trialPVector,trialPMatrix,trialMcEnergy,trialAwParameters);
            if(mcTargetReached){
                mcProposedMoves=move+1;
                break;
            }
        }
    }
    else{
        for(int move=0; move<mcMaxMoves; ++move){
            trialPVector=pVector;
            trialPMatrix=pMatrix;
            switchTriangles=pickRandomTrianglePairAperiodic();
            calculateTrialPAperiodic(switchTriangles, trialPVector, trialPMatrix);
            trialAwParameters=calculateAboavWeaireFit(trialPVector,trialPMatrix);
            trialMcEnergy=mcEnergyFunctional(trialAwParameters,trialPVector);
            acceptTrialMove=evaluateMetropolisCondition(trialMcEnergy,mcEnergy);
            if(acceptTrialMove) mcTargetReached=acceptDualSwitch(switchTriangles,trialPVector,trialPMatrix,trialMcEnergy,trialAwParameters);
            if(mcTargetReached){
                mcProposedMoves=move+1;
                break;
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
    for(int i=0; i<nRingSizes; ++i) e=e+fabs(norm*pVec[i]-targetPVector[i])*rTargetPVector[i];

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

bool Network::acceptDualSwitch(vector<int> &switchTriangles, vector<int> &trialPVec, vector<vector<int> > &trialPMat, double &trialMcEnergy, vector<double> &trialAwParams) {
    //enact dual switch on nodes and update trial-> current variables

    //make/break connections
    nodes[switchTriangles[0]].breakConnection(switchTriangles[1]);
    nodes[switchTriangles[1]].breakConnection(switchTriangles[0]);
    nodes[switchTriangles[2]].makeConnection(switchTriangles[3]);
    nodes[switchTriangles[3]].makeConnection(switchTriangles[2]);

    pVector=trialPVec;
    pMatrix=trialPMat;
    mcEnergy=trialMcEnergy;
    aboavWeaireParams=trialAwParams;

    cout<<aboavWeaireParams[0]<<" "<<aboavWeaireParams[1]<<" "<<aboavWeaireParams[2]<<" "<<mcEnergy<<endl;
//    consoleVector(switchTriangles);
//    checkFidelity();

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

//###### Write out ######

void Network::write() {
    //write out to files

    writeDual();

    if(periodicVisualisation) writePeriodicNetwork();

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

void Network::writePeriodicNetwork() {
    //only for visualisation, calculate periodic images of network for visualisation

    //set limits of visualisation region
    double imageProportion=0.5;
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
    Crd2d imageVec;
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
                imageVec=vectorFromCrds(periodicNodes[i].coordinate,periodicNodes[images[k]].coordinate);
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


    return;
}