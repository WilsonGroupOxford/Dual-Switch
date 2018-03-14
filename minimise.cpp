#include "minimise.h"

HarmonicMinimiser::HarmonicMinimiser() {
    //default invokes tests
    runTests();
    return;
}

HarmonicMinimiser::HarmonicMinimiser(vector<Crd2d> crds, vector<Pair> harmPairs, vector<int> fixedPnts,
                                     vector<double> harmR0, double harmK, double cc, double inc, int maxIt, vector<DoublePair> lineInt) {
    coordinates=crds;
    harmonicPairs=harmPairs;
    fixedPoints=fixedPnts;
    harmonicR0=harmR0;
    harmonicK=harmK;
    convergenceCriteria=cc;
    lineSearchIncrement=inc;
    maxIterations=maxIt;
    intersectionPairs=lineInt;
    nPnts=coordinates.size();
    nHarmonicPairs=harmonicPairs.size();
    nFixedPnts=fixedPnts.size();
    nIntersectionPairs=intersectionPairs.size();
    zeroForce=Crd2d();
    zeroForces.resize(nPnts,Crd2d());
    return;
}

void HarmonicMinimiser::steepestDescent() {
    //perform steepest descent until energy converged or hits iteration limit

    iterations=0; //intialise iteration counter
    previousEnergy=numeric_limits<double>::infinity(); //set initial energy to infinite
    deltaEZeroCount=0; //initialise no energy change counter
    complete=false; //calculation complete
    converged=false; //calculation converged
    for(;;){
        calculateForces();
        lineSearch();
        checkConvergence();
//        cout<<iterations<<" "<<currentEnergy<<endl;
        if(complete) break;
    }
    return;
}

void HarmonicMinimiser::calculateForces() {
    //calculate force on each point due to harmonic pairs - set to zero if a fixed point
    forces=zeroForces;
    //calculate force from each harmonic pair
    Crd2d force;
    int a, b;
    for(int i=0; i<nHarmonicPairs; ++i){
        a=harmonicPairs[i].a;
        b=harmonicPairs[i].b;
        force=harmonicForce(coordinates[a], coordinates[b], harmonicR0[i], harmonicK);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    //kill forces on fixed points
    for(int i=0; i<nFixedPnts; ++i) forces[fixedPoints[i]]=zeroForce;
    return;
}

double HarmonicMinimiser::calculateEnergy() {
    //calculate energy of harmonic interactions
    double energy=0.0;
    bool intersection=false;
    for(int i=0; i<nIntersectionPairs; ++i){
        intersection=properIntersectionLines(coordinates[intersectionPairs[i].a],coordinates[intersectionPairs[i].b],
                                coordinates[intersectionPairs[i].c], coordinates[intersectionPairs[i].d]);
        if(intersection) return numeric_limits<double>::infinity();
    }
    for(int i=0; i<nHarmonicPairs; ++i){
        energy=energy+harmonicEnergy(coordinates[harmonicPairs[i].a], coordinates[harmonicPairs[i].b], harmonicR0[i], harmonicK);
    }
    return energy;
}

Crd2d HarmonicMinimiser::harmonicForce(Crd2d &c1, Crd2d &c2, double &r0, double &k) {
    //F=-k(r-r0)
    Crd2d f, fDir; //force and direction
    double fMag, r; //magnitude of force, distance between points
    fDir.x=c2.x-c1.x;
    fDir.y=c2.y-c1.y;
    r=sqrt(fDir.x*fDir.x+fDir.y*fDir.y);
    fMag=k*(r-r0)/r; //divide by r to cancel length built into fDir
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

double HarmonicMinimiser::harmonicEnergy(Crd2d &c1, Crd2d &c2, double &r0, double &k) {
    //U=0.5k(r-r0)^2
    double r=sqrt(pow((c2.x-c1.x),2)+pow((c2.y-c1.y),2));
    return 0.5*k*pow((r-r0),2);
}

void HarmonicMinimiser::lineSearch() {
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

void HarmonicMinimiser::checkConvergence() {
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

vector<Crd2d> HarmonicMinimiser::getMinimisedCoordinates() {
    //return minimised coordinates
    return coordinates;
}

void HarmonicMinimiser::runTests() {
    //test minimiser

}