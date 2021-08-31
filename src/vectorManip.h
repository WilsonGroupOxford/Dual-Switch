//some useful vector manipulations

#ifndef DUAL_SWITCH_VECTORMANIP_H
#define DUAL_SWITCH_VECTORMANIP_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

bool checkVectorForValue(vector<int> &v, int &value);
vector<int> getCommonValuesBetweenVectors(vector<int> v1, vector<int> v2);
void removeValueFromVectorByRef(vector<int> &vec, int &value);
void removeValuesFromVectorByRef(vector<int> &vec, vector<int> &values);
void addValuesToVectorByRef(vector<int> &vec, vector<int> &values);
void sortVectorByRef(vector<int> &v);
vector<int> sortVectorRemoveDuplicates(vector<int> v1);
#endif //DUAL_SWITCH_VECTORMANIP_H
