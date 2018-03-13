#include "crd.h"

double crdDistance(Crd2d &c1, Crd2d &c2){
    //distance between two coordiates in 2d
    return sqrt((c2.x-c1.x)*(c2.x-c1.x)+(c2.y-c1.y)*(c2.y-c1.y));
}
double vectorLength(Crd2d &v){
    //length of vector
    return sqrt(v.x*v.x+v.y*v.y);
}

Crd2d vectorFromCrds(Crd2d &c1, Crd2d &c2){
    //vector c2-c1
    return {(c2.x-c1.x),(c2.y-c1.y)};
}
Crd2d vectorFromCrds(Crd2d &c1, Crd2d &c2, double &x, double &y, double &rX, double &rY){
    //vector c2-c1 with periodic boundary conditions
    Crd2d vec((c2.x-c1.x),(c2.y-c1.y));
    vec.x=vec.x-x*round(vec.x*rX);
    vec.y=vec.y-y*round(vec.y*rY);
    return vec;
}

//##### IMPLEMENTED USING COMPUTATIONAL GEOMETRY IN C #######
double signedAreaSqTriangle(Crd2d &a, Crd2d &b, Crd2d &c){
    //returns square of the signed area of a triangle defined by three points
    return a.x*b.y-a.y*b.x+a.y*c.x-a.x*c.y+b.x*c.y-c.x*b.y;
}
bool leftTriangle(Crd2d &a, Crd2d &b, Crd2d &c){
    //returns true if the signed area squared of a triangle is positive
    return signedAreaSqTriangle(a,b,c)>0.0;
}
bool collinearPoints(Crd2d &a, Crd2d &b, Crd2d &c, double threshold){
    //returns true if points are colliear within given threshold
    return fabs(signedAreaSqTriangle(a,b,c))<threshold;
}
bool properIntersectionLines(Crd2d &l1a, Crd2d &l1b, Crd2d &l2a, Crd2d &l2b){
    //returns true if lines properly intersect i.e. no points collinear
    if(collinearPoints(l1a,l1b,l2a) || collinearPoints(l1a,l1b,l2b) || collinearPoints(l2a,l2b,l1a) || collinearPoints(l2a,l2b,l1b)) return false;
    return (leftTriangle(l1a,l1b,l2a)^leftTriangle(l1a,l1b,l2b)) && (leftTriangle(l2a,l2b,l1a)^leftTriangle(l2a,l2b,l1b));
}
// ########## END ##########

vector<double> leastSquaresLinearRegression(vector<Crd2d> &data){
    //perform least squares linear regression and return gradient, intercept and r-squared
    vector<double> coefficients(3); //gradient, intercept and r-squared
    int n=data.size();
    double sumX=0.0, sumY=0.0, sumXY=0.0, sumXX=0.0, sumYY=0.0;
    for(int i=0; i<n; ++i){
        sumX=sumX+data[i].x;
        sumY=sumY+data[i].y;
        sumXY=sumXY+data[i].x*data[i].y;
        sumXX=sumXX+data[i].x*data[i].x;
        sumYY=sumYY+data[i].y*data[i].y;
    }
    double sqSumX=sumX*sumX, sqSumY=sumY*sumY;
    double sdX=sqrt(n*sumXX-sqSumX), sdY=sqrt(n*sumYY-sqSumY);
    double r=(n*sumXY-sumX*sumY)/(sdX*sdY);
    coefficients[0]=(r*sdY)/sdX;
    coefficients[1]=(sumY-coefficients[0]*sumX)/n;
    coefficients[2]=r*r;
    return coefficients;
}


