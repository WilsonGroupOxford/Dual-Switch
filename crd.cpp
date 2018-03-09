#include "crd.h"

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


