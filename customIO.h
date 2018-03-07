//simplify reading and writing to files/console

#ifndef DUAL_SWITCH_CUSTOMIO_H
#define DUAL_SWITCH_CUSTOMIO_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <sstream>

using namespace std;

//functions which write to the console
void consoleValue(int value);
void consoleVector(vector<int> &vec);
void consoleVector(vector<double> &vec);
void consoleVector(vector<bool> &vec);
void consoleMatrix(vector< vector<int> > &mat);
void consoleMatrix(vector< vector<double> > &mat);

//functions which read from file
void readFileSkipLines(ifstream &file, int nLines=1);
void readFileValue(ifstream &file, string &value);
void readFileValue(ifstream &file, int &value);
void readFileValue(ifstream &file, bool &value);
void readFileValue(ifstream &file, double &value);
void readFileRowVector(ifstream &file, vector<int> &vec, int nCols);
void readFileRowVector(ifstream &file, vector<double> &vec, int nCols);
void readFileArbitraryColsEOF(ifstream &file, vector< vector<double> > &readVector, int skipCol);
void readFileArbitraryColsEOF(ifstream &file, vector< vector<int> > &readVector, int skipCol);
void readFileColVector(ifstream &file, vector<int> &vec, int nRows);
void readFileMatrix(ifstream &file, vector< vector<double> > &mat, int nRows, int nCols);
void readFileAdaptiveMatrix(ifstream &file, vector< vector<int> > &mat, int nRows);

//functions which write to file
void writeFileLine(ofstream &openFile, string text, bool lineBreak=true);
void writeFileLine(ofstream &openFile, int value, bool lineBreak=true);
void writeFileLine(ofstream &openFile, int value, int w, string text, bool lineBreak=true);
void writeFileValue(ofstream &openFile, double value);
void writeFileRowVector(ofstream &openFile, vector<int> vec);
void writeFileRowVector(ofstream &openFile, vector<double> vec);
void writeFileColVector(ofstream &openFile, vector<int> vec);
void writeFileColVector(ofstream &openFile, vector<double> vec);
void writeFileMatrix(ofstream &openFile, vector< vector<int> > mat);
void writeFileMatrix(ofstream &openFile, vector< vector<double> > mat);
void writeFileVectorAndMatrix(ofstream &openFile, vector<int> vec, vector< vector<int> > mat);
void writeFilePassFail(ofstream &openFile, bool flag);
void writeFilePassFail(ofstream &openFile, string condition, bool flag);
void writeFileTrueFalse(ofstream &openFile, string condition, bool tf);
void writeFileDashes(ofstream &openFile, int n=70);

#endif //DUAL_SWITCH_CUSTOMIO_H
