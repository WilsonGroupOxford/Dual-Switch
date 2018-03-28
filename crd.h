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
    Crd2d();
    Crd2d(double xInit, double yInit);
};

//2 dimensional integer pair
struct Pair{
    int a, b;
    Pair();
    Pair(int aInit, int bInit);
};

//4 dimensional integer pair
struct DoublePair{
    int a, b, c, d;
    DoublePair();
    DoublePair(int aInit, int bInit, int cInit, int dInit);
    void sort();
    bool checkPairs();
    string getID();
};

//2 dimensional integer pair
struct Triangle{
    int a, b, c;
    Triangle();
    Triangle(int aInit, int bInit, int cInit);
    void sort();
    string getID();
};

//2 dimensional vector
struct Vec2d{
    double x, y;
    Vec2d();
    Vec2d(double xInit, double yInit);
    Vec2d(Crd2d c1, Crd2d c2);
    Vec2d(Crd2d c1, Crd2d c2, double pbcX, double pbcY, double rX, double rY);
    double length();
    void normalise();
    void rotate90();
    void scale(double s);
    void invert();
    void addVector(Vec2d v);
};





//constants
double const pi=M_PI;

double crdDistance(Crd2d &c1, Crd2d &c2);
Crd2d recentreCrdByCrd(Crd2d &c1, Crd2d &c2);
Crd2d recentreCrdByCrd(Crd2d &c1, Crd2d &c2, double &x, double &y, double &rX, double &rY);
Crd2d crdFromVectorAndCrd(Vec2d &v, Crd2d &c);
Crd2d minimumImageCrd(Crd2d &c1, Crd2d &c2, double &pbcX, double &pbcY, double &rPbcX, double &rPbcY);
void applyPeriodicBoundary(Crd2d &c1, double &pbcX, double &pbcY, double &rPbcX, double &rPbcY);
double vectorDotProduct(Vec2d &v1, Vec2d &v2);
Crd2d crdCentreOfMass(vector<Crd2d> crds);

//computational geometry in c
double signedAreaSqTriangle(Crd2d &a, Crd2d &b, Crd2d &c);
bool leftTriangle(Crd2d &a, Crd2d &b, Crd2d &c);
bool collinearPoints(Crd2d &a, Crd2d &b, Crd2d &c, double threshold=0.00001);
bool properIntersectionLines(Crd2d &l1a, Crd2d &l1b, Crd2d &l2a, Crd2d &l2b);
bool improperIntersectionLines(Crd2d &l1a, Crd2d &l1b, Crd2d &l2a, Crd2d &l2b);
bool betweenPoints(Crd2d &a, Crd2d &b, Crd2d &c);
bool lineIntersection(Crd2d &l1a, Crd2d &l1b, Crd2d &l2a, Crd2d &l2b);

//other useful maths functions
vector<double> leastSquaresLinearRegression(vector<Crd2d> &data);

#endif //DUAL_SWITCH_CRD_H
