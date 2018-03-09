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


