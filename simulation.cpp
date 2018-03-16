#include "simulation.h"

void Simulation::run() {
    //control all simulations

    readInputFile(); //get simulation parameters
    ofstream logfile(logfileName,ios::in|ios::app); //initialise logfile

    //perform network generation and analysis for the given random seeds
    writeFileLine(logfile,"Network Generation and Analysis");
    writeFileDashes(logfile);
    networkConsistent.clear();
    networkTargetsReached.clear();
    for(int seed=randomSeedLimits[0]; seed<=randomSeedLimits[1]; ++seed){
        cout<<seed<<endl;

        Network network;
        initialiseNetwork(network,seed);
        network.construct(logfile);

        networkConsistent.push_back(network.getConsistency());
        networkTargetsReached.push_back(network.getConsistency());
        networkIntersectionFree.push_back(network.getIntersectionStatus());

        if(networkConsistent.rbegin()[0]){
            network.analyse(logfile);
            network.write();
        }


    }
    logfile.close();
}

void Simulation::readInputFile() {
    //read parameters from input file

    //main input file
    ifstream mainInputFile("dual_switch.inpt", ios::in); //open input file
    readFileSkipLines(mainInputFile,1);
    readFileValue(mainInputFile,inPrefix);
    readFileValue(mainInputFile,outPrefix);
    readFileSkipLines(mainInputFile,2);
    readFileValue(mainInputFile,periodic);
    readFileValue(mainInputFile,load);
    readFileRowVector(mainInputFile,latticeDimensions,2);
    readFileRowVector(mainInputFile,ringSizeLimits,2);
    readFileValue(mainInputFile,alpha);
    readFileRowVector(mainInputFile,pVector,ringSizeLimits[1]-ringSizeLimits[0]+1);
    readFileSkipLines(mainInputFile,2);
    readFileRowVector(mainInputFile,randomSeedLimits,2);
    readFileValue(mainInputFile,temperature);
    readFileValue(mainInputFile,maxMoves);
    readFileValue(mainInputFile,propConvergence);
    readFileValue(mainInputFile,alphaEnergyScaling);
    readFileSkipLines(mainInputFile,2);
    readFileValue(mainInputFile,atomDistance);
    readFileValue(mainInputFile,harmonicK);
    readFileValue(mainInputFile,localGeomOpt);
    readFileValue(mainInputFile,localGeomOptMaxIt);
    readFileValue(mainInputFile,localGeomOptConv);
    readFileValue(mainInputFile,lineSearchStep);
    readFileValue(mainInputFile,globalGeomOpt);
    readFileValue(mainInputFile,globalGeomOptMaxIt);
    readFileValue(mainInputFile,globalGeomOptConv);
    readFileSkipLines(mainInputFile,2);
    readFileValue(mainInputFile,periodicVis);
    readFileValue(mainInputFile,spatialRdf);
    readFileValue(mainInputFile,spatialRdfBinWidth);
    readFileValue(mainInputFile,spatialRdfExtent);
    readFileValue(mainInputFile,topoRdf);
    readFileValue(mainInputFile,maxTopoShells);
    mainInputFile.close();

    logfileName=outPrefix+".log";
    ofstream logfile(logfileName,ios::in|ios::trunc);
    writeFileLine(logfile,"Dual Switch Calculation "+outPrefix);
    writeFileDashes(logfile);
    writeFileLine(logfile,"Input Parameters");
    writeFileDashes(logfile);
    string line;
    if(periodic) line="periodic ";
    else line="aperiodic ";
    line=line+"lattice of dimensions "+to_string(latticeDimensions[0])+"x"+to_string(latticeDimensions[1]);
    writeFileLine(logfile,line);
    line="ring size limits: "+to_string(ringSizeLimits[0])+"-"+to_string(ringSizeLimits[1]);
    writeFileLine(logfile,line);
    writeFileLine(logfile,"target alpha: "+to_string(alpha));
    writeFileLine(logfile,"target ring statistics:");
    writeFileRowVector(logfile,pVector);
    line="maximum "+to_string(maxMoves)+" moves of Monte Carlo dual switching at temperature "+to_string(temperature);
    writeFileLine(logfile,line);
    if(localGeomOpt) writeFileLine(logfile,"stepwise local geometry optimisation enabled");
    else writeFileLine(logfile,"stepwise local geometry optimisation disabled");
    if(globalGeomOpt) writeFileLine(logfile,"final global geometry optimisation enabled");
    else writeFileLine(logfile,"final global geometry optimisation disabled");
    writeFileDashes(logfile);

    logfile.close();

    return;
}

void Simulation::initialiseNetwork(Network &network, int seed) {
    //set up network with input properties and given random seed

    network.setIO(inPrefix,outPrefix+"_"+to_string(seed));
    network.setProperties(periodic,load,latticeDimensions,ringSizeLimits,alpha,pVector);
    network.setPotential(atomDistance, harmonicK, localGeomOpt, localGeomOptMaxIt, localGeomOptConv,
                         globalGeomOpt, globalGeomOptMaxIt, globalGeomOptConv, lineSearchStep);
    network.setMonteCarlo(seed,temperature,maxMoves,propConvergence,alphaEnergyScaling);
    network.setAnalysis(periodicVis,spatialRdf,spatialRdfBinWidth,spatialRdfExtent);

    return;
}