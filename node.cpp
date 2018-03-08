#include "node.h"

Node::Node() {
    //default constructor of 6 connected node at origin
    size=6;
    minSize=6;
    maxSize=10;
    sizeIndex=size-minSize;
    edge=false;
    coordinate=Crd2d();
    connections.clear();
    return;
}

Node::Node(int s, int minS, int maxS, bool e=false) {
    //constructor to set initial size
    size=s;
    minSize=minS;
    maxSize=maxS;
    sizeIndex=size-minSize;
    edge=e;
    coordinate=Crd2d();
    connections.clear();
    return;
}

void Node::addConnection(int cnx) {
    //add connection without changing size
    connections.push_back(cnx);
    return;
}

bool Node::makeConnection(int cnx) {
    //make new connection and increase node size if doesn't violate limits
    bool flag=true;
    if(size+1<=maxSize){
        connections.push_back(cnx);
        size=++size;
        sizeIndex=++sizeIndex;
    }
    else flag=false;
    return flag;
}

bool Node::breakConnection(int cnx) {
    //break old connection and decrease node size if doesn't violate limits
    bool flag=true;
    if(size-1>=minSize){
        connections.erase(remove(connections.begin(), connections.end(), cnx), connections.end());
        size=--size;
        sizeIndex=--sizeIndex;
    }
    else flag=false;
    return flag;
}
