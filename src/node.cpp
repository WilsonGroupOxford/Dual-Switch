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

Node::Node(int s, int minS, int maxS, bool e) {
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

void Node::makeConnection(int cnx) {
    //make new connection and increase node size/sizeIndex
    connections.push_back(cnx);
    size=++size;
    sizeIndex=++sizeIndex;
    return;
}

void Node::breakConnection(int cnx) {
    //break old connection and decrease node size/sizeIndex
    connections.erase(remove(connections.begin(), connections.end(), cnx), connections.end());
    size=--size;
    sizeIndex=--sizeIndex;
    return;
}

Vertex::Vertex() {
    //default constructor
    coordination=0;
    coordinate=Crd2d();
    for(int i=0; i<3; ++i){
        connections[i]=-1;
        dualNodes[i]=-1;
    }
    return;
}

Vertex::Vertex(int node0, int node1, int node2) {
    //constructor with dual nodes set
    coordination=0;
    coordinate=Crd2d();
    for(int i=0; i<3; ++i) connections[i]=-1;
    dualNodes[0]=node0;
    dualNodes[1]=node1;
    dualNodes[2]=node2;
    return;
}

void Vertex::addConnection(int cnx) {
    //add connection and increase coordination, if not already present
    if(connections[0]!=cnx && connections[1]!=cnx && connections[2]!=cnx) {
        connections[coordination] = cnx;
        coordination = ++coordination;
    }
    return;
}

Ring::Ring() {
    //default constructor
    size=6;
    id=-1;
    chain=new int[size+1](); //+1 as one node repeated at start and end
    return;
}

Ring::Ring(int s, vector<int> c, int d) {
    //construct of given size with provided nodes
    size=s;
    id=d;
    chain=new int[size+1]();
    for(int i=0; i<size; ++i) chain[i]=c[i];
    chain[size]=c[0];
    return;
}

vector<int> Ring::path() {
    //return chain, minus last element, as vector
    vector<int> path;
    path.clear();
    for(int i=0; i<size; ++i) path.push_back(chain[i]);
    return path;
}

Rdf::Rdf() {
    //default constructor
    id="None";
    binwidth=0.1;
    extent=10.0;
    rBinwidth=1.0/binwidth;
    nBins=floor(extent*rBinwidth);
    histogram=new int[nBins]();
    return;
}

Rdf::Rdf(double bw, double ext, string name) {
    //initialise with bin widths and extent
    id=name;
    binwidth=bw;
    extent=ext;
    rBinwidth=1.0/binwidth;
    nBins=int(floor(extent*rBinwidth));
    histogram=new int[nBins]();
    return;
}

void Rdf::addValue(double value) {
    //add value to histogram
    int bin=int(floor(value*rBinwidth));
    if(bin<nBins) histogram[bin]=++histogram[bin];
    return;
}

vector<int> Rdf::vectorHistogram() {
    //return vector form of histogram
    vector<int> vecHistogram(nBins);
    for(int i=0; i<nBins; ++i) vecHistogram[i]=histogram[i];
    return vecHistogram;
}

void Rdf::print() {
    //write to screen
    cout<<id<<endl;
    for(int i=0; i<nBins; ++i) cout<<histogram[i]<<endl;
    return;
}