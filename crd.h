//coordinate stuc and manipulations

#ifndef DUAL_SWITCH_CRD_H
#define DUAL_SWITCH_CRD_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

//2 dimensional coordinate
struct Crd2d{
    double x, y;
    Crd2d() : x(0.0), y(0.0) {}
    Crd2d(double xInit, double yInit) : x(xInit), y(yInit) {}
};

//constants
double const pi=M_PI;

//coordinate manipulations
double crdDistance(Crd2d &c1, Crd2d &c2);
double vectorLength(Crd2d &v);
Crd2d vectorFromCrds(Crd2d &c1, Crd2d &c2);
Crd2d vectorFromCrds(Crd2d &c1, Crd2d &c2, double &x, double &y, double &rX, double &rY);

//computational geometry in c
double signedAreaSqTriangle(Crd2d &a, Crd2d &b, Crd2d &c);
bool leftTriangle(Crd2d &a, Crd2d &b, Crd2d &c);
bool collinearPoints(Crd2d &a, Crd2d &b, Crd2d &c, double threshold=0.00001);
bool properIntersectionLines(Crd2d &l1a, Crd2d &l1b, Crd2d &l2a, Crd2d &l2b);

//other useful maths functions
vector<double> leastSquaresLinearRegression(vector<Crd2d> &data);

#endif //DUAL_SWITCH_CRD_H
