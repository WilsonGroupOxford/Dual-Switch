#include "geometryopt.h"

//###### BASE CLASS ######

GeometryOptimiser::GeometryOptimiser(vector<Crd2d> crds, double convCriteria, double lineInc, int maxIt){
    coordinates=crds;
    convergenceCriteria=convCriteria;
    lineSearchIncrement=lineInc;
    maxIterations=maxIt;
    nPnts=coordinates.size();
    zeroForce=Crd2d();
    zeroForces.resize(nPnts,Crd2d());
    return;
}

void GeometryOptimiser::setSystemParameters(vector<Pair> sysBonds, vector<Trio> sysAngles, vector<int> sysFixed,
                                            vector<DoublePair> sysInt) {
    //set system parameters, and get number of each interaction type
    bonds=sysBonds;
    angles=sysAngles;
    fixed=sysFixed;
    intersectionLines=sysInt;

    nBonds=bonds.size();
    nAngles=angles.size();
    nFixed=fixed.size();
    nIntersectionLines=intersectionLines.size();
    return;
}

int GeometryOptimiser::steepestDescent() {
    //steepest descent optimisation

    //check/resolve initial structure and if invalid do not perform optimisation
    if(checkInitialStructure()==0) return 0;

    iterations=0; //intialise iteration counter
    previousEnergy=numeric_limits<double>::infinity(); //set initial energy to infinite
    deltaEZeroCount=0; //initialise no energy change counter
    complete=false; //calculation complete
    converged=false; //calculation converged
    for(;;){
        calculateForces();
        lineSearch();
        checkConvergence();
        cout<<iterations<<" "<<currentEnergy<<endl;
        if(complete) break;
    }

    return 1;
}

void GeometryOptimiser::calculateForces() {
    //calculate force on each point due to each interaction type

    forces=zeroForces; //set forces to zero
    calculateBondForces(); //forces due to bonds
    calculateAngleForces(); //forces due to angles
    for(int i=0; i<nFixed; ++i) forces[fixed[i]]=zeroForce; //kill forces on fixed points

    return;
}

void GeometryOptimiser::lineSearch() {
    //move points along direction of force, recalculate energy and find minimum
    double e0, e1;
    e0=calculateEnergy();
    for(;;){
        for(int i=0; i<nPnts; ++i){//update positions
            coordinates[i].x=coordinates[i].x+lineSearchIncrement*forces[i].x;
            coordinates[i].y=coordinates[i].y+lineSearchIncrement*forces[i].y;
        }
        e1=calculateEnergy(); //calculate new energy
        if(e1-e0>1e-8){//if energy increases, passed through minimum so backtrack
            for(int i=0; i<nPnts; ++i){//update positions
                coordinates[i].x=coordinates[i].x-lineSearchIncrement*forces[i].x;
                coordinates[i].y=coordinates[i].y-lineSearchIncrement*forces[i].y;
            }
            currentEnergy=calculateEnergy();
            break;
        }
        else e0=e1;
    }
    return;
}

double GeometryOptimiser::calculateEnergy() {
    //calculate energy due to line intersections, bonds and angles
    double energy=0.0;
    if(checkIntersections()) return numeric_limits<double>::infinity(); //if intersections infinite energy
    energy=energy+calculateBondEnergies();
    energy=energy+calculateAngleEnergies();
    return energy;
}

void GeometryOptimiser::checkConvergence() {
    //check if energy converged, or hit maximum iterations

    if(iterations==maxIterations){//complete but not converged
        complete=true;
        return;
    }
    iterations=++iterations;
    if(currentEnergy<convergenceCriteria){//complete and converged
        complete=true;
        converged=true;
        return;
    }
    double deltaE=previousEnergy-currentEnergy;
    previousEnergy=currentEnergy;
    if(deltaE<1e-6) deltaEZeroCount=++deltaEZeroCount;
    else deltaEZeroCount=0; //reset counter
    if(deltaEZeroCount==10){//complete and converged and non-zero
        complete=true;
        converged=true;
        return;
    }
    return;
}

vector<Crd2d> GeometryOptimiser::getMinimisedCoordinates() {
    //return minimised coordinates
    return coordinates;
}

double GeometryOptimiser::getEnergy() {
    //return final minimised energy
    return currentEnergy;
}

int GeometryOptimiser::getIterations() {
    //return number of iterations of minimisation
    return iterations;
}