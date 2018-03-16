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
    for(int i=0; i<nBins; ++i) cout<<histogram[i]<<endl;
    return;
}