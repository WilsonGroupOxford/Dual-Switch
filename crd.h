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

struct Crd2d{
    double x, y;
    Crd2d() : x(0.0), y(0.0) {}
    Crd2d(double xInit, double yInit) : x(xInit), y(yInit) {}
};

double const pi=M_PI;

vector<double> leastSquaresLinearRegression(vector<Crd2d> &data);

#endif //DUAL_SWITCH_CRD_H
