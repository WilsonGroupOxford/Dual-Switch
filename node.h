//node in network, with coordinates, connections

#ifndef DUAL_SWITCH_NODE_H
#define DUAL_SWITCH_NODE_H

#include <vector>
#include "crd.h"

using namespace std;

class Node {
private:


public:

    //Node properties
    int size, sizeIndex, minSize, maxSize; //number of connections, index of size, min and max sizes
    bool edge; //whether on edge of network
    Crd2d coordinate; //2d coordinate
    vector<int> connections; //to other nodes


    //Constructors
    Node();
    Node(int s, int minS, int maxS, bool e=false); //set initial size (and edge)

    //Connection manipulation
    void addConnection(int cnx); //add without changing size
    void makeConnection(int cnx); //make increments size
    void breakConnection(int cnx); //break decrements size


};


#endif //DUAL_SWITCH_NODE_H
