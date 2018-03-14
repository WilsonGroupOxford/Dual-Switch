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

//###### Steepest descent without PBCs ######

void HarmonicMinimiser::steepestDescent() {
    //perform steepest descent until energy converged or hits iteration limit

    cout<<"Starting intersection "<<checkIntersections()<<endl;
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
    if(checkIntersections()) return numeric_limits<double>::infinity(); //check if line intersections
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

bool HarmonicMinimiser::checkIntersections() {
    bool intersection=false;
    for(int i=0; i<nIntersectionPairs; ++i){
        intersection=properIntersectionLines(coordinates[intersectionPairs[i].a],coordinates[intersectionPairs[i].b],
                                             coordinates[intersectionPairs[i].c], coordinates[intersectionPairs[i].d]);
        if(intersection) return true;
    }
    return intersection;
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

//###### Steepest Descent With PBCs ######

void HarmonicMinimiser::steepestDescent(double x, double y, double rx, double ry) {
    //perform steepest descent until energy converged or hits iteration limit, overloaded for periodic boundary conditions

    //setup periodic boundary variables
    pbcX=x;
    pbcY=y;
    pbcRX=rx;
    pbcRY=ry;
    cout<<"Starting intersection "<<checkIntersectionsPBC()<<endl;
    iterations=0; //intialise iteration counter
    previousEnergy=numeric_limits<double>::infinity(); //set initial energy to infinite
    deltaEZeroCount=0; //initialise no energy change counter
    complete=false; //calculation complete
    converged=false; //calculation converged
    for(;;){
        calculateForcesPBC();
        lineSearchPBC();
        checkConvergence();
//        cout<<iterations<<" "<<currentEnergy<<endl;
        if(complete) break;
    }
    wrapAroundCoordinates(); //apply periodic bcs to coordiates

    return;
}

void HarmonicMinimiser::calculateForcesPBC() {
    //calculate force on each point due to harmonic pairs - set to zero if a fixed point
    forces=zeroForces;
    //calculate force from each harmonic pair
    Crd2d force;
    int a, b;
    for(int i=0; i<nHarmonicPairs; ++i){
        a=harmonicPairs[i].a;
        b=harmonicPairs[i].b;
        force=harmonicForcePBC(coordinates[a], coordinates[b], harmonicR0[i], harmonicK);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    //kill forces on fixed points
    for(int i=0; i<nFixedPnts; ++i) forces[fixedPoints[i]]=zeroForce;
    return;
}

Crd2d HarmonicMinimiser::harmonicForcePBC(Crd2d &c1, Crd2d &c2, double &r0, double &k) {
    //F=-k(r-r0)
    Crd2d f, fDir; //force and direction
    double fMag, r; //magnitude of force, distance between points
    fDir.x=c2.x-c1.x;
    fDir.y=c2.y-c1.y;
    fDir.x=fDir.x-pbcX*round(fDir.x*pbcRX); //apply pbc
    fDir.y=fDir.y-pbcY*round(fDir.y*pbcRY); //apply pbc
    r=sqrt(fDir.x*fDir.x+fDir.y*fDir.y);
    fMag=k*(r-r0)/r; //divide by r to cancel length built into fDir
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void HarmonicMinimiser::lineSearchPBC() {
    //move points along direction of force, recalculate energy and find minimum
    double e0, e1;
    e0=calculateEnergyPBC();
    for(;;){
        for(int i=0; i<nPnts; ++i){//update positions
            coordinates[i].x=coordinates[i].x+lineSearchIncrement*forces[i].x;
            coordinates[i].y=coordinates[i].y+lineSearchIncrement*forces[i].y;
        }
        e1=calculateEnergyPBC(); //calculate new energy
        if(e1-e0>1e-8){//if energy increases, passed through minimum so backtrack
            for(int i=0; i<nPnts; ++i){//update positions
                coordinates[i].x=coordinates[i].x-lineSearchIncrement*forces[i].x;
                coordinates[i].y=coordinates[i].y-lineSearchIncrement*forces[i].y;
            }
            currentEnergy=calculateEnergyPBC();
            break;
        }
        else e0=e1;
    }
    return;
}

double HarmonicMinimiser::calculateEnergyPBC() {
    //calculate energy of harmonic interactions
    double energy=0.0;
    if(checkIntersectionsPBC()) return numeric_limits<double>::infinity(); //check if line intersections
    for(int i=0; i<nHarmonicPairs; ++i){
        energy=energy+harmonicEnergyPBC(coordinates[harmonicPairs[i].a], coordinates[harmonicPairs[i].b], harmonicR0[i], harmonicK);
    }
    return energy;
}

double HarmonicMinimiser::harmonicEnergyPBC(Crd2d &c1, Crd2d &c2, double &r0, double &k) {
    //U=0.5k(r-r0)^2
    Crd2d c3((c2.x-c1.x),(c2.y-c1.y));
    c3.x=c3.x-pbcX*round(c3.x*pbcRX);
    c3.y=c3.y-pbcY*round(c3.y*pbcRY);
    double r=sqrt(c3.x*c3.x+c3.y*c3.y);
    return 0.5*k*pow((r-r0),2);
}

bool HarmonicMinimiser::checkIntersectionsPBC() {
    bool intersection=false;
    Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
    for(int i=0; i<nIntersectionPairs; ++i){
        line1a=coordinates[intersectionPairs[i].a];
        line1b=minimumImageCrd(line1a,coordinates[intersectionPairs[i].b],pbcX,pbcY,pbcRX,pbcRY);
        line2a=minimumImageCrd(line1a,coordinates[intersectionPairs[i].c],pbcX,pbcY,pbcRX,pbcRY);
        line2b=minimumImageCrd(line1a,coordinates[intersectionPairs[i].d],pbcX,pbcY,pbcRX,pbcRY);
        intersection=properIntersectionLines(line1a,line1b,line2a,line2b);
        if(intersection) return true;
    }
    return intersection;
}

void HarmonicMinimiser::wrapAroundCoordinates() {
    //apply periodic boundary conditions to coordinates
    for(int i=0; i<nPnts; ++i){
        coordinates[i].x=coordinates[i].x-pbcX*round(coordinates[i].x*pbcRX);
        coordinates[i].y=coordinates[i].y-pbcY*round(coordinates[i].y*pbcRY);
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