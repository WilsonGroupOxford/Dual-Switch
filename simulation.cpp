#include "simulation.h"

void Simulation::run() {
    //control all simulations

    readInputFile(); //get simulation parameters
    ofstream logfile(logfileName,ios::in|ios::app); //initialise logfile

    //perform network generation and analysis for the given random seeds
    writeFileLine(logfile,"Network Generation and Analysis");
    writeFileDashes(logfile);
    ofstream seedStatusFile(outPrefix+"_"+to_string(randomSeedLimits[0])+"_"+to_string(randomSeedLimits[1])+"_seed_status.out",ios::in|ios::trunc); //initialise status file
    writeFileLine(seedStatusFile,"seed           consistent     target         intersection free");
    vector<int> seedStatus(4); //record whether network consistent, met target and is intersection free
    ofstream timingFile(outPrefix+"_timing.out",ios::in|ios::trunc); //initialise timing file
    writeFileLine(timingFile,"seed    minutes    seconds");
    chrono::high_resolution_clock::time_point timeBegin, timeEnd; //time points for each seed
    vector<int> seedTime(3); //time to complete each seed in minutes and seconds
    for(int seed=randomSeedLimits[0]; seed<=randomSeedLimits[1]; ++seed){
        timeBegin=chrono::high_resolution_clock::now();
        cout<<"Seed "<<seed<<endl;
        Network network;
        initialiseNetwork(network,seed);
        network.construct(logfile);

        seedStatus[0]=seed;
        seedStatus[1]=network.getConsistency();
        seedStatus[2]=network.getTargetStatus();
        seedStatus[3]=network.getIntersectionStatus();

        if(seedStatus[1]==1){//continue if consistent
            network.analyse(logfile);
            network.write();
        }

        writeFileRowVector(seedStatusFile, seedStatus);

        timeEnd=chrono::high_resolution_clock::now();
        int time=chrono::duration_cast<chrono::seconds>(timeEnd-timeBegin).count();
        seedTime[0]=seed;
        seedTime[1]=time/60;
        seedTime[2]=time-seedTime[1]*60;
        writeFileRowVector(timingFile,seedTime);
    }
    seedStatusFile.close();
    timingFile.close();
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
    readFileValue(mainInputFile,convertDual);
    readFileValue(mainInputFile,periodicVis);
    readFileValue(mainInputFile,spatialRdf);
    readFileValue(mainInputFile,spatialRdfBinWidth);
    readFileValue(mainInputFile,spatialRdfExtent);
    readFileValue(mainInputFile,topoRdf);
    readFileValue(mainInputFile,topoRdfExtent);
    readFileValue(mainInputFile,assortative);
    readFileSkipLines(mainInputFile,2);
    readFileValue(mainInputFile,atomicGeomOpt);
    readFileValue(mainInputFile,keatingA);
    readFileValue(mainInputFile,keatingAlpha);
    readFileValue(mainInputFile,keatingBeta);
    readFileValue(mainInputFile,atomGeomOptMaxIt);
    readFileValue(mainInputFile,atomGeomOptConv);
    readFileValue(mainInputFile,atomLineSearchStep);
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
    network.setAnalysis(convertDual,periodicVis,spatialRdf,spatialRdfBinWidth,spatialRdfExtent,topoRdf,topoRdfExtent,assortative);
    network.setAtomicPotential(atomicGeomOpt,keatingA,keatingAlpha,keatingBeta,atomGeomOptMaxIt,atomGeomOptConv,atomLineSearchStep);
    return;
}