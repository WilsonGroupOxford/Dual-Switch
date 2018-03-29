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

int HarmonicMinimiser::steepestDescent() {
    //perform steepest descent until energy converged or hits iteration limit

    if(!resolveInitialIntersections()) return 0;

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

bool HarmonicMinimiser::resolveInitialIntersections() {

    bool resolved=!checkIntersections();
    if(resolved) return true; //return if no intersections

    //get intersecting lines
    vector<DoublePair> intersectingLines=getIntersections();
    int nIntersectingLines=intersectingLines.size();

    //find main line common to all
    int majorIntersectionCount;
    Pair majorIntersection=findMajorIntersection(intersectingLines, majorIntersectionCount);

    //abort if not one major line
    if(majorIntersectionCount!=nIntersectingLines) return false;

    //get list of unique points
    vector<int> uniquePoints=findUniquePoints(intersectingLines,majorIntersection);

    //resolve intersections
    resolved=moveIntersectingPoints(uniquePoints, majorIntersection, nIntersectingLines);

    return resolved;
}

vector<DoublePair> HarmonicMinimiser::getIntersections() {
    vector<DoublePair> intersectingLines; //stored intersecting pairs
    intersectingLines.clear();
    for(int i=0; i<nIntersectionPairs; ++i){
        if(properIntersectionLines(coordinates[intersectionPairs[i].a],coordinates[intersectionPairs[i].b],
                                   coordinates[intersectionPairs[i].c], coordinates[intersectionPairs[i].d])){
            intersectingLines.push_back(intersectionPairs[i]);
        };
    }
    return intersectingLines;
}

Pair HarmonicMinimiser::findMajorIntersection(vector<DoublePair> &intersectingLines, int &majorCount) {
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

vector<int> HarmonicMinimiser::findUniquePoints(vector<DoublePair> &intersectingLines, Pair &majorLine) {
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

bool HarmonicMinimiser::moveIntersectingPoints(vector<int> &uniquePoints, Pair &majorIntersection,
                                               int &nIntersections) {
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

//###### Steepest Descent With PBCs ######

int HarmonicMinimiser::steepestDescent(double x, double y, double rx, double ry) {
    //perform steepest descent until energy converged or hits iteration limit, overloaded for periodic boundary conditions

    //setup periodic boundary variables
    pbcX=x;
    pbcY=y;
    pbcRX=rx;
    pbcRY=ry;
    if(!resolveInitialIntersectionsPBC()) return 0;
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

    return 1;
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
    fDir.x=fDir.x-pbcX*round(fDir.x*pbcRX); //apply mic
    fDir.y=fDir.y-pbcY*round(fDir.y*pbcRY); //apply mic
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
        coordinates[i].x=coordinates[i].x-pbcX*floor(coordinates[i].x*pbcRX);
        coordinates[i].y=coordinates[i].y-pbcY*floor(coordinates[i].y*pbcRY);
    }
    return;
}

bool HarmonicMinimiser::resolveInitialIntersectionsPBC() {

    bool resolved=!checkIntersectionsPBC();
    if(resolved) return true; //return if no intersections

    //get intersecting lines
    vector<DoublePair> intersectingLines=getIntersectionsPBC();
    int nIntersectingLines=intersectingLines.size();

    //find main line common to all
    int majorIntersectionCount;
    Pair majorIntersection=findMajorIntersection(intersectingLines, majorIntersectionCount);

    //abort if not one major line
    if(majorIntersectionCount!=nIntersectingLines) return false;

    //get list of unique points
    vector<int> uniquePoints=findUniquePoints(intersectingLines,majorIntersection);

    //resolve intersections
    resolved=moveIntersectingPointsPBC(uniquePoints, majorIntersection, nIntersectingLines);

    return resolved;
}

vector<DoublePair> HarmonicMinimiser::getIntersectionsPBC() {
    vector<DoublePair> intersectingLines; //stored intersecting pairs
    intersectingLines.clear();
    Crd2d line1a, line1b, line2a, line2b; //coordinates as minimum images to crd 1a
    for(int i=0; i<nIntersectionPairs; ++i){
        line1a=coordinates[intersectionPairs[i].a];
        line1b=minimumImageCrd(line1a,coordinates[intersectionPairs[i].b],pbcX,pbcY,pbcRX,pbcRY);
        line2a=minimumImageCrd(line1a,coordinates[intersectionPairs[i].c],pbcX,pbcY,pbcRX,pbcRY);
        line2b=minimumImageCrd(line1a,coordinates[intersectionPairs[i].d],pbcX,pbcY,pbcRX,pbcRY);
        if(properIntersectionLines(line1a,line1b,line2a,line2b)) intersectingLines.push_back(intersectionPairs[i]);
    }
    return intersectingLines;
}

bool HarmonicMinimiser::moveIntersectingPointsPBC(vector<int> &uniquePoints, Pair &majorIntersection,
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
    bool moveSuccessful=!checkIntersectionsPBC();
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
    moveSuccessful=!checkIntersectionsPBC();
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
        moveSuccessful=!checkIntersectionsPBC();
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
        moveSuccessful=!checkIntersectionsPBC();
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

//###### GETTERS ######

vector<Crd2d> HarmonicMinimiser::getMinimisedCoordinates() {
    //return minimised coordinates
    return coordinates;
}

double HarmonicMinimiser::getEnergy() {
    //return final minimised energy
    return currentEnergy;
}

int HarmonicMinimiser::getIterations() {
    //return number of iterations of minimisation
    return iterations;
}

void HarmonicMinimiser::runTests() {
    //test minimiser

}

//###### Keating Minimiser ######

KeatingMinimiser::KeatingMinimiser(vector<Crd2d> crds, double cc, double inc, int maxIt) {
    coordinates=crds;
    convergenceCriteria=cc;
    lineSearchIncrement=inc;
    maxIterations=maxIt;
    nPnts=coordinates.size();
    zeroForce=Crd2d();
    zeroForces.resize(nPnts,Crd2d());
    return;
}

void KeatingMinimiser::setParameters(double kpA, double kpAlpha, double kpBeta) {
    //initialise keating potential parameters

    a=kpA;
    alpha=kpAlpha;
    beta=kpBeta;

    aSq=a*a;
    aSq_2=aSq/2.0;
    kBondForce=(3.0/4.0)*(alpha/aSq);
    kAngleForce=(3.0/4.0)*(beta/aSq);
    kBondEnergy=(3.0/16.0)*(alpha/aSq);
    kAngleEnergy=(3.0/8.0)*(beta/aSq);

    return;
}

void KeatingMinimiser::setInteractions(vector<Pair> kpBonds, vector<Trio> kpAngles, vector<DoublePair> ints) {
    //lists of atoms forming bonds, angles, and line intersections

    bonds=kpBonds;
    angles=kpAngles;
    intersectionPairs=ints;

    nBonds=bonds.size();
    nAngles=angles.size();
    nIntersectionPairs=intersectionPairs.size();
    return;
}

int KeatingMinimiser::steepestDescent() {
    //aperiodic steepest descent optimisation

    //check if intersections, if so return without minimising
    if(checkIntersections()) return 0;

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

void KeatingMinimiser::calculateForces() {
    //calculate force on each point due to harmonic pairs - set to zero if a fixed point
    forces=zeroForces;
    //calculate force due to bonds
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
    //calculate force due to angles
    Crd2d forceI, forceJ, forceK; //force on points i--j--k
    int c;
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

Crd2d KeatingMinimiser::bondForce(Crd2d &c1, Crd2d &c2) {
    //f=-k(r^2-a^2)r
    Crd2d f, fDir; //force and direction
    double fMag; //magnitude of force
    fDir.x=c2.x-c1.x;
    fDir.y=c2.y-c1.y;
    fMag=kBondForce*(fDir.x*fDir.x+fDir.y*fDir.y-aSq);
    f.x=-fMag*fDir.x;
    f.y=-fMag*fDir.y;
    return f;
}

void KeatingMinimiser::angleForce(Crd2d &cI, Crd2d &cJ, Crd2d &cK, Crd2d &fI, Crd2d &fJ, Crd2d &fK) {
    //f=-k(r1*r2+a^2/2)r1r2
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

void KeatingMinimiser::lineSearch() {
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

double KeatingMinimiser::calculateEnergy() {
    //calculate energy of harmonic interactions
    double energy=0.0;
    if(checkIntersections()) return numeric_limits<double>::infinity();
    for(int i=0; i<nBonds; ++i) energy=energy+bondEnergy(coordinates[bonds[i].a], coordinates[bonds[i].b]);
    for(int i=0; i<nAngles; ++i) energy=energy+angleEnergy(coordinates[angles[i].a], coordinates[angles[i].b], coordinates[angles[i].c]);
    return energy;
}

double KeatingMinimiser::bondEnergy(Crd2d &c1, Crd2d &c2) {
    //u=k*(r^2-a^2)^2
    Vec2d v(c1,c2);
    return kBondEnergy*pow((v.x*v.x+v.y*v.y-aSq),2);
}

double KeatingMinimiser::angleEnergy(Crd2d &cI, Crd2d &cJ, Crd2d &cK) {
    //u=k*(r1.r2+a^2/2)^2
    Vec2d vij(cJ,cI), vkj(cJ,cK); //convention Rij=ri-rj
    return kAngleEnergy*pow((vij.x*vkj.x+vij.y*vkj.y+aSq_2),2);
}

void KeatingMinimiser::checkConvergence() {
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

bool KeatingMinimiser::checkIntersections() {
    bool intersection=false;
    for(int i=0; i<nIntersectionPairs; ++i){
        intersection=properIntersectionLines(coordinates[intersectionPairs[i].a],coordinates[intersectionPairs[i].b],
                                             coordinates[intersectionPairs[i].c], coordinates[intersectionPairs[i].d]);
        if(intersection) return true;
    }
    return intersection;
}

//###### GETTERS ######
vector<Crd2d> KeatingMinimiser::getMinimisedCoordinates() {
    //return minimised coordinates
    return coordinates;
}