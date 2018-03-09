#include "vectorManip.h"

bool checkVectorForValue(vector<int> &v, int &value){
    //check to see if a vector contains a value
    if(find(v.begin(), v.end(), value) != v.end()) {
        return true;
    }
    else return false;
}