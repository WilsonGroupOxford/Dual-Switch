#include "geometryopt.h"

//###### OPTIMISATION BASE CLASS ######

GeometryOptimiser::GeometryOptimiser() {
    return;
}

void GeometryOptimiser::setCoordinates(vector<Crd2d> crds) {
    //set coordinates for minimisation
    coordinates=crds;
    return;
}

void GeometryOptimiser::setOptisationParameters(double convCriteria, double lineInc, int maxIt, bool res) {
    //set parameters for optimisation functions
    convergenceCriteria=convCriteria;
    lineSearchIncrement=lineInc;
    maxIterations=maxIt;
    nPnts=coordinates.size();
    resolveStructure=res;
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

int GeometryOptimiser::checkInitialStructure() {
    //check no initial line overlaps
    bool resolved=!checkIntersections();
    if(resolved) return 1;
    if(resolveStructure) {
        resolved = resolveIntersections();
        if (resolved) return 1;
    }
    return 0;
}

bool GeometryOptimiser::resolveIntersections() {
    //attempt to resolve intersections, return success or failure

    //get intersecting lines
    vector<DoublePair> intersectingLines = getIntersections();
    int nIntersectingLines = intersectingLines.size();

    //find main line common to all
    int majorIntersectionCount;
    Pair majorIntersection = findMajorIntersection(intersectingLines, majorIntersectionCount);

    //abort if not one major line
    if (majorIntersectionCount != nIntersectingLines) return false;

    //get list of unique points
    vector<int> uniquePoints = findUniquePoints(intersectingLines, majorIntersection);

    //resolve intersections
    return moveIntersectingPoints(uniquePoints, majorIntersection, nIntersectingLines);
}

Pair GeometryOptimiser::findMajorIntersection(vector<DoublePair> &intersectingLines, int &majorCount) {
    //find intersecting line common to all
    map<string,int> lineMap;
    Pair line1, line2;
    string code1, code2;
    for(int i=0; i<intersectingLines.size();++i){//count number of instances of each line
        line1=Pair(intersectingLines[i].a,intersectingLines[i].b);
        line2=Pair(intersectingLines[i].c,intersectingLines[i].d);
        code1="#"+to_string(line1.a)+"#"+to_string(line1.b);
        code2="#"+to_string(line2.a)+"#"+to_string(line2.b);
        if(!lineMap.count(code1)) lineMap[code1]=1;
        else lineMap[code1]=++lineMap[code1];
        if(!lineMap.count(code2)) lineMap[code2]=1;
        else lineMap[code2]=++lineMap[code2];
    }
    int maxCount=0;
    Pair majorLine;
    for(int i=0; i<intersectingLines.size();++i){//find most common line
        line1=Pair(intersectingLines[i].a,intersectingLines[i].b);
        line2=Pair(intersectingLines[i].c,intersectingLines[i].d);
        code1="#"+to_string(line1.a)+"#"+to_string(line1.b);
        code2="#"+to_string(line2.a)+"#"+to_string(line2.b);
        if(lineMap[code1]>maxCount){
            majorLine=line1;
            maxCount=lineMap[code1];
        }
        if(lineMap[code2]>maxCount){
            majorLine=line2;
            maxCount=lineMap[code2];
        }
    }
    majorCount=maxCount;
    return majorLine;
}

vector<int> GeometryOptimiser::findUniquePoints(vector<DoublePair> &intersectingLines, Pair &majorLine) {
    //find unique list of points which might need to move
    vector<int> points;
    points.clear();
    int p1, p2;
    for(int i=0; i<intersectingLines.size(); ++i){//loop over lines
        if(intersectingLines[i].a==majorLine.a){//pick points from line which is not the major line
            p1=intersectingLines[i].c;
            p2=intersectingLines[i].d;
        }
        else{
            p1=intersectingLines[i].a;
            p2=intersectingLines[i].b;
        }
        if(!checkVectorForValue(points,p1)) points.push_back(p1);
        if(!checkVectorForValue(points,p2)) points.push_back(p2);
    }
    return points;
}

//###### APERIODIC OPTIMISATION BASE CLASS ######

AperiodicGeometryOptimiser::AperiodicGeometryOptimiser() {
    return;
}

int AperiodicGeometryOptimiser::steepestDescent() {
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
//        cout<<iterations<<" "<<currentEnergy<<endl;
        if(complete) break;
    }

    return 1;
}

bool AperiodicGeometryOptimiser::checkIntersections() {
    //check for overlap of lines
    bool intersection=false;
    for(int i=0; i<nIntersectionLines; ++i){
        intersection=properIntersectionLines(coordinates[intersectionLines[i].a],coordinates[intersectionLines[i].b],
                                             coordinates[intersectionLines[i].c], coordinates[intersectionLines[i].d]);
        if(intersection) return true;
    }
    return intersection;
}

vector<DoublePair> AperiodicGeometryOptimiser::getIntersections() {
    vector<DoublePair> intersectingLines; //stored intersecting pairs
    intersectingLines.clear();
    for(int i=0; i<nIntersectionLines; ++i){
        if(properIntersectionLines(coordinates[intersectionLines[i].a],coordinates[intersectionLines[i].b],
                                   coordinates[intersectionLines[i].c], coordinates[intersectionLines[i].d])){
            intersectingLines.push_back(intersectionLines[i]);
        };
    }
    return intersectingLines;
}

bool AperiodicGeometryOptimiser::moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection, int &nIntersections) {
    //move points above and below major intersection and check if resolved intersections

    //set up axes relative to major intersection
    Vec2d xAxis(coordinates[majorIntersection.a],coordinates[majorIntersection.b]);
    xAxis.normalise();
    Vec2d yAxis=xAxis;
    yAxis.rotate90();

    //find points above and below x-axis (major intersection)
    vector<int> pointsAboveX, pointsBelowX;
    vector<Crd2d> initialCrdsAboveX, initialCrdsBelowX; //save intial coordinates if need to revert
    Vec2d posVec;
    pointsAboveX.clear();
    pointsBelowX.clear();
    initialCrdsAboveX.clear();
    initialCrdsBelowX.clear();
    double dotProd;
    for(int i=0; i<uniquePoints.size(); ++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[uniquePoints[i]]);
        dotProd=vectorDotProduct(yAxis,posVec);
        if(dotProd>0.0){
            pointsAboveX.push_back(uniquePoints[i]);
            initialCrdsAboveX.push_back(coordinates[uniquePoints[i]]);
        }
        else{
            pointsBelowX.push_back(uniquePoints[i]);
            initialCrdsBelowX.push_back(coordinates[uniquePoints[i]]);
        }
    }

    //try moving points above x to below x
    for(int i=0; i<pointsAboveX.size();++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsAboveX[i]]);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsAboveX[i]].x=coordinates[pointsAboveX[i]].x+direction.x;
        coordinates[pointsAboveX[i]].y=coordinates[pointsAboveX[i]].y+direction.y;
    }

    //check if successful and if not revert and move points below x to above x
    bool moveSuccessful=!checkIntersections();
    if(moveSuccessful) return true;

    for(int i=0; i<pointsAboveX.size(); ++i) coordinates[pointsAboveX[i]]=initialCrdsAboveX[i]; //revert
    for(int i=0; i<pointsBelowX.size();++i){//move points below x to above
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsBelowX[i]]);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsBelowX[i]].x=coordinates[pointsBelowX[i]].x+direction.x;
        coordinates[pointsBelowX[i]].y=coordinates[pointsBelowX[i]].y+direction.y;
    }

    //check if successful and if not revert and try moving slightly in x-direction
    moveSuccessful=!checkIntersections();
    if(moveSuccessful) return true;


    for (int i = 0; i < pointsBelowX.size(); ++i) coordinates[pointsBelowX[i]] = initialCrdsBelowX[i]; //revert
    //move above to below again, but shift in x
    Vec2d directionInc=xAxis;
    directionInc.scale(0.1);
    for(int i=0; i<pointsAboveX.size();++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsAboveX[i]]);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsAboveX[i]].x=coordinates[pointsAboveX[i]].x+direction.x-xAxis.x;
        coordinates[pointsAboveX[i]].y=coordinates[pointsAboveX[i]].y+direction.y-xAxis.y;
    }
    //increment along in x direction
    for(int j=0; j<=20; ++j){
        moveSuccessful=!checkIntersections();
        if(moveSuccessful) return true;
        for(int i=0; i<pointsAboveX.size();++i){
            coordinates[pointsAboveX[i]].x=coordinates[pointsAboveX[i]].x+directionInc.x;
            coordinates[pointsAboveX[i]].y=coordinates[pointsAboveX[i]].y+directionInc.y;
        }
    }

    //if not successful try same with points below
    for(int i=0; i<pointsAboveX.size(); ++i) coordinates[pointsAboveX[i]]=initialCrdsAboveX[i]; //revert
    //move above to below again, but shift in x
    for(int i=0; i<pointsBelowX.size();++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsBelowX[i]]);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsBelowX[i]].x=coordinates[pointsBelowX[i]].x+direction.x-xAxis.x;
        coordinates[pointsBelowX[i]].y=coordinates[pointsBelowX[i]].y+direction.y-xAxis.y;
    }
    //increment along in x direction
    for(int j=0; j<=20; ++j){
        moveSuccessful=!checkIntersections();
        if(moveSuccessful) return true;
        for(int i=0; i<pointsBelowX.size();++i){
            coordinates[pointsBelowX[i]].x=coordinates[pointsBelowX[i]].x+directionInc.x;
            coordinates[pointsBelowX[i]].y=coordinates[pointsBelowX[i]].y+directionInc.y;
        }
    }

    //if not successful revert
    for(int i=0; i<pointsBelowX.size(); ++i) coordinates[pointsBelowX[i]]=initialCrdsBelowX[i]; //revert

    return false;
}

//###### PERIODIC OPTIMISATION BASE CLASS ######

PeriodicGeometryOptimiser::PeriodicGeometryOptimiser() {
    return;
}

void PeriodicGeometryOptimiser::setPeriodicBoundary(double x, double y, double rx, double ry) {
    //periodic boundary conditions for minimum image convention and wrapping
    pbcX=x;
    pbcY=y;
    pbcRX=rx;
    pbcRY=ry;
    return;
}

int PeriodicGeometryOptimiser::steepestDescent() {
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
//        cout<<iterations<<" "<<currentEnergy<<endl;
        if(complete) break;
    }

    wrapCoordinates();

    return 1;
}

bool PeriodicGeometryOptimiser::checkIntersections() {
    bool intersection=false;
    Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
    for(int i=0; i<nIntersectionLines; ++i){
        line1a=coordinates[intersectionLines[i].a];
        line1b=minimumImageCrd(line1a,coordinates[intersectionLines[i].b],pbcX,pbcY,pbcRX,pbcRY);
        line2a=minimumImageCrd(line1a,coordinates[intersectionLines[i].c],pbcX,pbcY,pbcRX,pbcRY);
        line2b=minimumImageCrd(line2a,coordinates[intersectionLines[i].d],pbcX,pbcY,pbcRX,pbcRY);
        intersection=properIntersectionLines(line1a,line1b,line2a,line2b);
        if(intersection) return true;
    }
    return intersection;
}

vector<DoublePair> PeriodicGeometryOptimiser::getIntersections() {
    vector<DoublePair> intersectingLines; //stored intersecting pairs
    intersectingLines.clear();
    Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
    for(int i=0; i<nIntersectionLines; ++i){
        line1a=coordinates[intersectionLines[i].a];
        line1b=minimumImageCrd(line1a,coordinates[intersectionLines[i].b],pbcX,pbcY,pbcRX,pbcRY);
        line2a=minimumImageCrd(line1a,coordinates[intersectionLines[i].c],pbcX,pbcY,pbcRX,pbcRY);
        line2b=minimumImageCrd(line2a,coordinates[intersectionLines[i].d],pbcX,pbcY,pbcRX,pbcRY);
        if(properIntersectionLines(line1a,line1b,line2a,line2b)) intersectingLines.push_back(intersectionLines[i]);
    }
    return intersectingLines;
}

bool PeriodicGeometryOptimiser::moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection,
                                               int &nIntersections) {
    //move points above and below major intersection and check if resolved intersections
    //take minimum images relative to major point a


    //set up axes relative to major intersection
    Vec2d xAxis(coordinates[majorIntersection.a],coordinates[majorIntersection.b],pbcX,pbcY,pbcRX,pbcRY);
    xAxis.normalise();
    Vec2d yAxis=xAxis;
    yAxis.rotate90();

    //find points above and below x-axis (major intersection)
    vector<int> pointsAboveX, pointsBelowX;
    vector<Crd2d> initialCrdsAboveX, initialCrdsBelowX; //save intial coordinates if need to revert
    Vec2d posVec;
    pointsAboveX.clear();
    pointsBelowX.clear();
    initialCrdsAboveX.clear();
    initialCrdsBelowX.clear();
    double dotProd;
    for(int i=0; i<uniquePoints.size(); ++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[uniquePoints[i]],pbcX,pbcY,pbcRX,pbcRY);
        dotProd=vectorDotProduct(yAxis,posVec);
        if(dotProd>0.0){
            pointsAboveX.push_back(uniquePoints[i]);
            initialCrdsAboveX.push_back(coordinates[uniquePoints[i]]);
        }
        else{
            pointsBelowX.push_back(uniquePoints[i]);
            initialCrdsBelowX.push_back(coordinates[uniquePoints[i]]);
        }
    }

    //try moving points above x to below x
    for(int i=0; i<pointsAboveX.size();++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsAboveX[i]],pbcX,pbcY,pbcRX,pbcRY);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsAboveX[i]].x=coordinates[pointsAboveX[i]].x+direction.x;
        coordinates[pointsAboveX[i]].y=coordinates[pointsAboveX[i]].y+direction.y;
    }

    //check if successful and if not revert and move points below x to above x
    bool moveSuccessful=!checkIntersections();
    if(moveSuccessful) return true;

    for(int i=0; i<pointsAboveX.size(); ++i) coordinates[pointsAboveX[i]]=initialCrdsAboveX[i]; //revert
    for(int i=0; i<pointsBelowX.size();++i){//move points below x to above
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsBelowX[i]],pbcX,pbcY,pbcRX,pbcRY);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsBelowX[i]].x=coordinates[pointsBelowX[i]].x+direction.x;
        coordinates[pointsBelowX[i]].y=coordinates[pointsBelowX[i]].y+direction.y;
    }

    //check if successful and if not revert and try moving slightly in x-direction
    moveSuccessful=!checkIntersections();
    if(moveSuccessful) return true;


    for (int i = 0; i < pointsBelowX.size(); ++i) coordinates[pointsBelowX[i]] = initialCrdsBelowX[i]; //revert
    //move above to below again, but shift in x
    Vec2d directionInc=xAxis;
    directionInc.scale(0.1);
    for(int i=0; i<pointsAboveX.size();++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsAboveX[i]],pbcX,pbcY,pbcRX,pbcRY);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsAboveX[i]].x=coordinates[pointsAboveX[i]].x+direction.x-xAxis.x;
        coordinates[pointsAboveX[i]].y=coordinates[pointsAboveX[i]].y+direction.y-xAxis.y;
    }
    //increment along in x direction
    for(int j=0; j<=20; ++j){
        moveSuccessful=!checkIntersections();
        if(moveSuccessful) return true;
        for(int i=0; i<pointsAboveX.size();++i){
            coordinates[pointsAboveX[i]].x=coordinates[pointsAboveX[i]].x+directionInc.x;
            coordinates[pointsAboveX[i]].y=coordinates[pointsAboveX[i]].y+directionInc.y;
        }
    }

    //if not successful try same with points below
    for(int i=0; i<pointsAboveX.size(); ++i) coordinates[pointsAboveX[i]]=initialCrdsAboveX[i]; //revert
    //move above to below again, but shift in x
    for(int i=0; i<pointsBelowX.size();++i){
        posVec=Vec2d(coordinates[majorIntersection.a],coordinates[pointsBelowX[i]],pbcX,pbcY,pbcRX,pbcRY);
        dotProd=vectorDotProduct(yAxis,posVec);
        Vec2d direction=yAxis;
        direction.scale(dotProd*1.01);
        direction.invert();
        coordinates[pointsBelowX[i]].x=coordinates[pointsBelowX[i]].x+direction.x-xAxis.x;
        coordinates[pointsBelowX[i]].y=coordinates[pointsBelowX[i]].y+direction.y-xAxis.y;
    }
    //increment along in x direction
    for(int j=0; j<=20; ++j){
        moveSuccessful=!checkIntersections();
        if(moveSuccessful) return true;
        for(int i=0; i<pointsBelowX.size();++i){
            coordinates[pointsBelowX[i]].x=coordinates[pointsBelowX[i]].x+directionInc.x;
            coordinates[pointsBelowX[i]].y=coordinates[pointsBelowX[i]].y+directionInc.y;
        }
    }

    //if not successful revert
    for(int i=0; i<pointsBelowX.size(); ++i) coordinates[pointsBelowX[i]]=initialCrdsBelowX[i]; //revert

    return false;
}

void PeriodicGeometryOptimiser::wrapCoordinates() {
    //apply periodic boundary conditions to coordinates
    for(int i=0; i<nPnts; ++i){
        coordinates[i].x=coordinates[i].x-pbcX*floor(coordinates[i].x*pbcRX);
        coordinates[i].y=coordinates[i].y-pbcY*floor(coordinates[i].y*pbcRY);
    }
    return;
}

/*###### HARMONIC APERIODIC OPTIMISER ######
 * takes single k and list of r0, no minimum image convention */

HarmonicAperiodicGO::HarmonicAperiodicGO(){
    return;
};

void HarmonicAperiodicGO::setPotentialParameters(double harmK, vector<double> harmR0) {
    //take single force constant and list of r0
    forceK=harmK;
    energyK=0.5*forceK;
    r0=harmR0;
    return;
}

void HarmonicAperiodicGO::calculateBondForces() {
    //calculate force from all bonds
    Crd2d force;
    int a, b;
    for(int i=0; i<nBonds; ++i){
        a=bonds[i].a;
        b=bonds[i].b;
        force=bondForce(coordinates[a], coordinates[b], r0[i]);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    return;
}

Crd2d HarmonicAperiodicGO::bondForce(Crd2d &cI, Crd2d &cJ, double &r0) {
    //force of single bond, F=-k(r-r0)
    Crd2d f; //force
    Vec2d fDir(cI,cJ); //force direction
    double fMag, r; //magnitude of force, distance between points
    r=sqrt(fDir.x*fDir.x+fDir.y*fDir.y);
    fMag=forceK*(r-r0)/r; //divide by r to cancel length built into fDir
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void HarmonicAperiodicGO::calculateAngleForces() {
    //no angle potential
    return;
}

double HarmonicAperiodicGO::calculateBondEnergies() {
    //calculate energies of all bonds, U=0.5*k*(r-r0)^2
    double energy=0.0;
    double r;
    Vec2d v;
    for(int i=0; i<nBonds; ++i){
        v=Vec2d(coordinates[bonds[i].a],coordinates[bonds[i].b]);
        r=sqrt(v.x*v.x+v.y*v.y);
        energy=energy+energyK*pow((r-r0[i]),2);
    }
    return energy;
}

double HarmonicAperiodicGO::calculateAngleEnergies() {
    //no angle contribution
    return 0.0;
}

/*###### HARMONIC PERIODIC OPTIMISER ######
 * takes single k and list of r0, periodic boundary conditions */

HarmonicPeriodicGO::HarmonicPeriodicGO(){
    return;
};

void HarmonicPeriodicGO::setPotentialParameters(double harmK, vector<double> harmR0) {
    //take single force constant and list of r0
    forceK=harmK;
    energyK=0.5*forceK;
    r0=harmR0;
    return;
}

void HarmonicPeriodicGO::calculateBondForces() {
    //calculate force from all bonds
    Crd2d force;
    int a, b;
    for(int i=0; i<nBonds; ++i){
        a=bonds[i].a;
        b=bonds[i].b;
        force=bondForce(coordinates[a], coordinates[b], r0[i]);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    return;
}

Crd2d HarmonicPeriodicGO::bondForce(Crd2d &cI, Crd2d &cJ, double &r0) {
    //force of single bond, F=-k(r-r0)
    Crd2d f; //force
    Vec2d fDir(cI,cJ); //force direction
    fDir.x=fDir.x-pbcX*round(fDir.x*pbcRX); //apply mic
    fDir.y=fDir.y-pbcY*round(fDir.y*pbcRY);
    double fMag, r; //magnitude of force, distance between points
    r=sqrt(fDir.x*fDir.x+fDir.y*fDir.y);
    fMag=forceK*(r-r0)/r; //divide by r to cancel length built into fDir
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void HarmonicPeriodicGO::calculateAngleForces() {
    //no angle potential
    return;
}

double HarmonicPeriodicGO::calculateBondEnergies() {
    //calculate energies of all bonds, U=0.5*k*(r-r0)^2
    double energy=0.0;
    double r;
    Vec2d v;
    for(int i=0; i<nBonds; ++i){
        v=Vec2d(coordinates[bonds[i].a],coordinates[bonds[i].b]);
        v.x=v.x-pbcX*round(v.x*pbcRX);
        v.y=v.y-pbcY*round(v.y*pbcRY);
        r=sqrt(v.x*v.x+v.y*v.y);
        energy=energy+energyK*pow((r-r0[i]),2);
    }
    return energy;
}

double HarmonicPeriodicGO::calculateAngleEnergies() {
    //no angle contribution
    return 0.0;
}

/*###### KEATING APERIODIC OPTIMISER ######
 * takes a, alpha, beta, no minimum image convention */

KeatingAperiodicGO::KeatingAperiodicGO(){
    return;
};

void KeatingAperiodicGO::setPotentialParameters(double a, double alpha, double beta) {
    //take single force constant and list of r0
    aSq=a*a;
    aSq_2=aSq/2.0;
    kBondForce=(3.0/4.0)*(alpha/aSq);
    kAngleForce=(3.0/4.0)*(beta/aSq);
    kBondEnergy=(3.0/16.0)*(alpha/aSq);
    kAngleEnergy=(3.0/8.0)*(beta/aSq);
    return;
}

void KeatingAperiodicGO::calculateBondForces() {
    //calculate force from all bonds
    Crd2d force;
    int a, b;
    for(int i=0; i<nBonds; ++i){
        a=bonds[i].a;
        b=bonds[i].b;
        force=bondForce(coordinates[a], coordinates[b]);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    return;
}

Crd2d KeatingAperiodicGO::bondForce(Crd2d &cI, Crd2d &cJ) {
    //force of single bond, F=-k(r^2-a^2)r
    Crd2d f; //force
    Vec2d fDir(cI,cJ); //force direction
    double fMag=kBondForce*(fDir.x*fDir.x+fDir.y*fDir.y-aSq); //magnitude of force
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void KeatingAperiodicGO::calculateAngleForces() {
    //calculate force due to angles
    Crd2d forceI, forceJ, forceK; //force on points i--j--k
    int a,b,c;
    for(int i=0; i<nAngles; ++i){
        a=angles[i].a;
        b=angles[i].b;
        c=angles[i].c;
        angleForce(coordinates[a],coordinates[b],coordinates[c],forceI,forceJ,forceK);
        forces[a].x=forces[a].x+forceI.x;
        forces[a].y=forces[a].y+forceI.y;
        forces[b].x=forces[b].x+forceJ.x;
        forces[b].y=forces[b].y+forceJ.y;
        forces[c].x=forces[c].x+forceK.x;
        forces[c].y=forces[c].y+forceK.y;
    }
    return;
}

void KeatingAperiodicGO::angleForce(Crd2d &cI, Crd2d &cJ, Crd2d &cK, Crd2d &fI, Crd2d &fJ, Crd2d &fK) {
    //force of single angle, f=-k(r1*r2+a^2/2)r1r2
    Vec2d vecRij(cJ,cI), vecRkj(cJ,cK); //vectors to edge points from central point, convention Rij=ri-rj
    double rij, rkj; //lengths of vectors
    vecRij.normalise(rij); //normalise and get length
    vecRkj.normalise(rkj); //normalise and get length
    double cosTheta=vecRij.x*vecRkj.x+vecRij.y*vecRkj.y; //get angle between vectors
    double fMag=-kAngleForce*(rij*rkj*cosTheta+aSq_2); //neglect r1r2 multiplication as need to divide by one in direction
    fI.x=(fMag*rkj)*(vecRkj.x-cosTheta*vecRij.x);
    fI.y=(fMag*rkj)*(vecRkj.y-cosTheta*vecRij.y);
    fK.x=(fMag*rij)*(vecRij.x-cosTheta*vecRkj.x);
    fK.y=(fMag*rij)*(vecRij.y-cosTheta*vecRkj.y);
    fJ.x=-fI.x-fK.x;
    fJ.y=-fI.y-fK.y;
    return;
}

double KeatingAperiodicGO::calculateBondEnergies() {
    //calculate energies of all bonds, U=k*(r^2-a^2)^2
    double energy=0.0;
    Vec2d v;
    for(int i=0; i<nBonds; ++i){
        v=Vec2d(coordinates[bonds[i].a],coordinates[bonds[i].b]);
        energy=energy+kBondEnergy*pow((v.x*v.x+v.y*v.y-aSq),2);
    }
    return energy;
}

double KeatingAperiodicGO::calculateAngleEnergies() {
    //calculate energies of all angles, U=k*(r1.r2+a^2/2)^2
    double energy=0.0;
    Vec2d vij, vkj; //convention Rij=ri-rj
    for(int i=0; i<nAngles; ++i){
        vij=Vec2d(coordinates[angles[i].b],coordinates[angles[i].a]);
        vkj=Vec2d(coordinates[angles[i].b],coordinates[angles[i].c]);
        energy=energy+kAngleEnergy*pow((vij.x*vkj.x+vij.y*vkj.y+aSq_2),2);
    }
    return energy;
}

/*###### KEATING PERIODIC OPTIMISER ######
 * takes a, alpha, beta, periodic boundary conditions */

KeatingPeriodicGO::KeatingPeriodicGO(){
    return;
};

void KeatingPeriodicGO::setPotentialParameters(double a, double alpha, double beta) {
    //take single force constant and list of r0
    aSq=a*a;
    aSq_2=aSq/2.0;
    kBondForce=(3.0/4.0)*(alpha/aSq);
    kAngleForce=(3.0/4.0)*(beta/aSq);
    kBondEnergy=(3.0/16.0)*(alpha/aSq);
    kAngleEnergy=(3.0/8.0)*(beta/aSq);
    return;
}

void KeatingPeriodicGO::calculateBondForces() {
    //calculate force from all bonds
    Crd2d force;
    int a, b;
    for(int i=0; i<nBonds; ++i){
        a=bonds[i].a;
        b=bonds[i].b;
        force=bondForce(coordinates[a], coordinates[b]);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    return;
}

Crd2d KeatingPeriodicGO::bondForce(Crd2d &cI, Crd2d &cJ) {
    //force of single bond, F=-k(r^2-a^2)r
    Crd2d f; //force
    Vec2d fDir(cI,cJ); //force direction
    fDir.x=fDir.x-pbcX*round(fDir.x*pbcRX); //apply mic
    fDir.y=fDir.y-pbcY*round(fDir.y*pbcRY);
    double fMag=kBondForce*(fDir.x*fDir.x+fDir.y*fDir.y-aSq); //magnitude of force
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void KeatingPeriodicGO::calculateAngleForces() {
    //calculate force due to angles
    Crd2d forceI, forceJ, forceK; //force on points i--j--k
    int a,b,c;
    for(int i=0; i<nAngles; ++i){
        a=angles[i].a;
        b=angles[i].b;
        c=angles[i].c;
        angleForce(coordinates[a],coordinates[b],coordinates[c],forceI,forceJ,forceK);
        forces[a].x=forces[a].x+forceI.x;
        forces[a].y=forces[a].y+forceI.y;
        forces[b].x=forces[b].x+forceJ.x;
        forces[b].y=forces[b].y+forceJ.y;
        forces[c].x=forces[c].x+forceK.x;
        forces[c].y=forces[c].y+forceK.y;
    }
    return;
}

void KeatingPeriodicGO::angleForce(Crd2d &cI, Crd2d &cJ, Crd2d &cK, Crd2d &fI, Crd2d &fJ, Crd2d &fK) {
    //force of single angle, f=-k(r1*r2+a^2/2)r1r2
    Vec2d vecRij(cJ,cI), vecRkj(cJ,cK); //vectors to edge points from central point, convention Rij=ri-rj
    vecRij.x=vecRij.x-pbcX*round(vecRij.x*pbcRX); //apply mic
    vecRij.y=vecRij.y-pbcY*round(vecRij.y*pbcRY);
    vecRkj.x=vecRkj.x-pbcX*round(vecRkj.x*pbcRX); //apply mic
    vecRkj.y=vecRkj.y-pbcY*round(vecRkj.y*pbcRY);
    double rij, rkj; //lengths of vectors
    vecRij.normalise(rij); //normalise and get length
    vecRkj.normalise(rkj); //normalise and get length
    double cosTheta=vecRij.x*vecRkj.x+vecRij.y*vecRkj.y; //get angle between vectors
    double fMag=-kAngleForce*(rij*rkj*cosTheta+aSq_2); //neglect r1r2 multiplication as need to divide by one in direction
    fI.x=(fMag*rkj)*(vecRkj.x-cosTheta*vecRij.x);
    fI.y=(fMag*rkj)*(vecRkj.y-cosTheta*vecRij.y);
    fK.x=(fMag*rij)*(vecRij.x-cosTheta*vecRkj.x);
    fK.y=(fMag*rij)*(vecRij.y-cosTheta*vecRkj.y);
    fJ.x=-fI.x-fK.x;
    fJ.y=-fI.y-fK.y;
    return;
}

double KeatingPeriodicGO::calculateBondEnergies() {
    //calculate energies of all bonds, U=k*(r^2-a^2)^2
    double energy=0.0;
    Vec2d v;
    for(int i=0; i<nBonds; ++i){
        v=Vec2d(coordinates[bonds[i].a],coordinates[bonds[i].b]);
        v.x=v.x-pbcX*round(v.x*pbcRX);
        v.y=v.y-pbcY*round(v.y*pbcRY);
        energy=energy+kBondEnergy*pow((v.x*v.x+v.y*v.y-aSq),2);
    }
    return energy;
}

double KeatingPeriodicGO::calculateAngleEnergies() {
    //calculate energies of all angles, U=k*(r1.r2+a^2/2)^2
    double energy=0.0;
    Vec2d vij, vkj; //convention Rij=ri-rj
    for(int i=0; i<nAngles; ++i){
        vij=Vec2d(coordinates[angles[i].b],coordinates[angles[i].a]);
        vkj=Vec2d(coordinates[angles[i].b],coordinates[angles[i].c]);
        vij.x=vij.x-pbcX*round(vij.x*pbcRX);
        vij.y=vij.y-pbcY*round(vij.y*pbcRY);
        vkj.x=vkj.x-pbcX*round(vkj.x*pbcRX);
        vkj.y=vkj.y-pbcY*round(vkj.y*pbcRY);
        energy=energy+kAngleEnergy*pow((vij.x*vkj.x+vij.y*vkj.y+aSq_2),2);
    }
    return energy;
}

/*###### KEATING APERIODIC BOND ONLY OPTIMISER ######
 * takes range of a, single alpha, only bond term not angle, no minimum image convention */

KeatingAperiodicBondOnlyGO::KeatingAperiodicBondOnlyGO(){
    return;
};

void KeatingAperiodicBondOnlyGO::setPotentialParameters(double alpha, vector<double> a) {
    //take list of minimum distances and alpha
    int n=a.size();
    aSq.resize(n);
    aSq_2.resize(n);
    kBondForce.resize(n);
    kBondEnergy.resize(n);
    for(int i=0; i<n; ++i){
        aSq[i]=a[i]*a[i];
        aSq_2[i]=0.5*aSq[i];
        kBondForce[i]=(3.0/4.0)*(alpha/aSq[i]);
        kBondEnergy[i]=(3.0/16.0)*(alpha/aSq[i]);
    }
    return;
}

void KeatingAperiodicBondOnlyGO::calculateBondForces() {
    //calculate force from all bonds
    Crd2d force;
    int a, b;
    for(int i=0; i<nBonds; ++i){
        a=bonds[i].a;
        b=bonds[i].b;
        force=bondForce(coordinates[a], coordinates[b], i);
        forces[a].x=forces[a].x-force.x;
        forces[a].y=forces[a].y-force.y;
        forces[b].x=forces[b].x+force.x;
        forces[b].y=forces[b].y+force.y;
    }
    return;
}

Crd2d KeatingAperiodicBondOnlyGO::bondForce(Crd2d &cI, Crd2d &cJ, int bond) {
    //force of single bond, F=-k(r^2-a^2)r
    Crd2d f; //force
    Vec2d fDir(cI,cJ); //force direction
    double fMag=kBondForce[bond]*(fDir.x*fDir.x+fDir.y*fDir.y-aSq[bond]); //magnitude of force
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void KeatingAperiodicBondOnlyGO::calculateAngleForces() {
    //no angle contribution
    return;
}

double KeatingAperiodicBondOnlyGO::calculateBondEnergies() {
    //calculate energies of all bonds, U=k*(r^2-a^2)^2
    double energy=0.0;
    Vec2d v;
    for(int i=0; i<nBonds; ++i){
        v=Vec2d(coordinates[bonds[i].a],coordinates[bonds[i].b]);
        energy=energy+kBondEnergy[i]*pow((v.x*v.x+v.y*v.y-aSq[i]),2);
    }
    return energy;
}

double KeatingAperiodicBondOnlyGO::calculateAngleEnergies() {
    //no angle contribution
    return 0.0;
}