#include "vectorManip.h"

bool checkVectorForValue(vector<int> &v, int &value){
    //check to see if a vector contains a value
    if(find(v.begin(), v.end(), value) != v.end()) {
        return true;
    }
    else return false;
}

vector<int> getCommonValuesBetweenVectors(vector<int> v1, vector<int> v2) {
    //get common values between two vectors
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    vector<int> common;
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(common));
    return common;
}

void removeValueFromVectorByRef(vector<int> &vec, int &value){
    //remove vector of values from vector
    vec.erase(remove(vec.begin(), vec.end(), value), vec.end());
}

void removeValuesFromVectorByRef(vector<int> &vec, vector<int> &values){
    //remove vector of values from vector
    for (int i=0; i<values.size();++i){
        vec.erase(remove(vec.begin(), vec.end(), values[i]), vec.end());
    }
}

void addValuesToVectorByRef(vector<int> &vec, vector<int> &values){
    //add vector of values to vector
    for (int i=0;i<values.size();++i) vec.push_back(values[i]);
}

void sortVectorByRef(vector<int> &v){
    //sort
    sort(v.begin(), v.end());
    return;
}

vector<int> sortVectorRemoveDuplicates(vector<int> v1){
    //sort and remove repeated values
    sortVectorByRef(v1);
    vector<int> v2;
    v2.clear();
    v2.push_back(v1[0]);
    for (int i=1;i<v1.size();++i)if (v1[i]!=v2.rbegin()[0]) v2.push_back(v1[i]);
    return v2;
}