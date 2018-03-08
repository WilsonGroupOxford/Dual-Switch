
#include "customIO.h"

void consoleValue(int value){
    //write value to console
    cout<<value<<endl;
    return;
}

void consoleVector(vector<int> &vec){
    //write vector to console
    for (int i=0;i<vec.size();++i){
        cout<<setw(15)<<left<<vec[i];
    }
    cout<<endl;
    return;
}
void consoleVector(vector<double> &vec){
    //write vector to console
    for (int i=0;i<vec.size();++i){
        cout<<setw(15)<<left<<vec[i];
    }
    cout<<endl;
    return;
}
void consoleVector(vector<bool> &vec){
    //write vector to console
    for (int i=0;i<vec.size();++i){
        cout<<setw(15)<<left<<vec[i];
    }
    cout<<endl;
    return;
}

void consoleMatrix(vector< vector<int> > &mat){
    //write matrix to console
    for (int i=0;i<mat.size();++i){
        for (int j=0; j<mat[i].size();++j){
            cout<<setw(15)<<left<<mat[i][j];
        }
        cout<<endl;
    }
    return;
}
void consoleMatrix(vector< vector<double> > &mat){
    //write matrix to console
    for (int i=0;i<mat.size();++i){
        for (int j=0; j<mat[i].size();++j){
            cout<<setw(15)<<left<<mat[i][j];
        }
        cout<<endl;
    }
    return;
}

void readFileSkipLines(ifstream &file, int nLines){
    //skip over line
    string line;
    for(int i=0; i<nLines;++i) getline(file,line);
    return;
}

void readFileValue(ifstream &file, string &value){
    //read value from first column
    string line;
    getline(file,line);
    istringstream ss(line);
    ss >> value;
}
void readFileValue(ifstream &file, int &value){
    //read value from first column
    string line;
    getline(file,line);
    istringstream ss(line);
    ss >> value;
}
void readFileValue(ifstream &file, bool &value){
    //read value from first column
    string line;
    getline(file,line);
    istringstream ss(line);
    ss >> value;
}
void readFileValue(ifstream &file, double &value){
    //read value from first column
    string line;
    getline(file,line);
    istringstream ss(line);
    ss >> value;
}

void readFileRowVector(ifstream &file, vector<int> &vec, int nCols){
    //read in vector of n columns
    vec.clear();
    string line;
    getline(file,line);
    istringstream ss(line);
    int val;
    for (int i=0; i<nCols; ++i){
        ss >> val;
        vec.push_back(val);
    }
}
void readFileRowVector(ifstream &file, vector<double> &vec, int nCols){
    //read in vector of n columns
    vec.clear();
    string line;
    getline(file,line);
    istringstream ss(line);
    double val;
    for (int i=0; i<nCols; ++i){
        ss >> val;
        vec.push_back(val);
    }
}
void readFileColVector(ifstream &file, vector<int> &vec, int nRows){
    //read in vector of n rows
    vec.clear();
    string line;
    int value;
    for(int i=0; i<nRows;++i){
        getline(file,line);
        istringstream ss(line);
        ss >> value;
        vec.push_back(value);
    }
    return;
}

void readFileMatrix(ifstream &file, vector< vector<double> > &mat, int nRows, int nCols){
    //read in vector of specified rows and columns
    mat.clear();
    string line;
    vector<double> row;
    for (int i=0; i<nRows; ++i){
        row.clear();
        getline(file,line);
        istringstream ss(line);
        double val;
        for (int j=0; j<nCols; ++j){
            ss >> val;
            row.push_back(val);
        }
        mat.push_back(row);
    }
}

void readFileAdaptiveMatrix(ifstream &file, vector< vector<int> > &mat, int nRows){
    //read in vector of specified rows and adaptive number of columns
    mat.clear();
    string line;
    vector<int> row;
    for (int i=0; i<nRows; ++i){
        row.clear();
        getline(file,line);
        istringstream ss(line);
        int value;
        while (ss >> value){
            row.push_back(value);
        }
        mat.push_back(row);
    }
}

vector< vector<int> > readFileRowsCols(ifstream &file, int nRows, int nCols){
    //read number of rows and columns
    vector< vector<int> > data;
    vector<int> row;
    string line;
    int val;
    data.clear();
    for (int i=0; i<nRows; ++i){
        getline(file,line);
        istringstream ss(line);
        row.clear();
        for (int j=0; j<nCols; ++j){
            ss >> val;
            row.push_back(val);
        }
        data.push_back(row);
    }
    return data;
}

void readFileArbitraryColsEOF(ifstream &file, vector< vector<double> > &readVector, int skipCol){
    readVector.clear();
    string line;
    vector<double> row;
    double val;
    int col;
    while(getline(file, line)){
        istringstream ss(line);
        row.clear();
        col=0;
        while (ss >> val) {
            if (col != skipCol) row.push_back(val);
            col = ++col;
        }
        readVector.push_back(row);
    }
}
void readFileArbitraryColsEOF(ifstream &file, vector< vector<int> > &readVector, int skipCol){
    readVector.clear();
    string line;
    vector<int> row;
    int val;
    int col;
    while(getline(file, line)){
        istringstream ss(line);
        row.clear();
        col=0;
        while (ss >> val) {
            if (col != skipCol) row.push_back(val);
            col = ++col;
        }
        readVector.push_back(row);
    }
}

void writeFileLine(ofstream &openFile, string text, bool lineBreak){
    //write line to file with or without break
    if(lineBreak) openFile<<text<<endl;
    else openFile<<text;
}
void writeFileLine(ofstream &openFile, int value, bool lineBreak){
    //write line to file with or without break
    if(lineBreak) openFile<<value<<endl;
    else openFile<<value;
}
void writeFileLine(ofstream &openFile, int value, int w, string text, bool lineBreak){
    //write line to file with or without break
    if(lineBreak) openFile<<setw(w)<<left<<value<<setw(w)<<left<<text<<endl;
    else openFile<<setw(w)<<left<<value<<setw(w)<<left<<text<<setw(w)<<left<<" ";
}

void writeFileValue(ofstream &openFile, double value){
    openFile<<value<<endl;
}

void writeFileRowVector(ofstream &openFile, vector<int> vec){
    //write vector to file
    for(int i=0; i<vec.size();++i) openFile<<setw(15)<<left<<vec[i];
    openFile<<endl;
}
void writeFileRowVector(ofstream &openFile, vector<double> vec){
    //write vector to file
    for(int i=0; i<vec.size();++i) openFile<<setw(15)<<left<<vec[i];
    openFile<<endl;
}
void writeFileColVector(ofstream &openFile, vector<int> vec){
    //write vector to file
    for(int i=0; i<vec.size();++i) openFile<<setw(15)<<left<<vec[i]<<endl;
}
void writeFileColVector(ofstream &openFile, vector<double> vec){
    //write vector to file
    for(int i=0; i<vec.size();++i) openFile<<setw(15)<<left<<vec[i]<<endl;
}

void writeFileMatrix(ofstream &openFile, vector< vector<int> > mat){
    //write matrix to file
    for(int i=0; i<mat.size();++i) {
        for(int j=0; j<mat[i].size();++j){
            openFile<<setw(10)<<left<<mat[i][j];
        }
        openFile<<endl;
    }
}
void writeFileMatrix(ofstream &openFile, vector< vector<double> > mat){
    //write matrix to file
    for(int i=0; i<mat.size();++i) {
        for(int j=0; j<mat[i].size();++j){
            openFile<<setw(15)<<left<<mat[i][j];
        }
        openFile<<endl;
    }
}
void writeFileVectorAndMatrix(ofstream &openFile, vector<int> vec, vector< vector<int> > mat){
    //write matrix to file
    for(int i=0; i<mat.size();++i) {
        openFile<<setw(15)<<left<<vec[i];
        for(int j=0; j<mat[i].size();++j){
            openFile<<setw(15)<<left<<mat[i][j];
        }
        openFile<<endl;
    }
}

void writeFilePassFail(ofstream &openFile, bool flag){
    if(flag) openFile<<"passed"<<endl;
    else openFile<<"failed"<<endl;
    return;
}
void writeFilePassFail(ofstream &openFile, string condition, bool flag){
    openFile<<condition<<": ";
    if(flag) openFile<<"passed"<<endl;
    else openFile<<"failed"<<endl;
    return;
}

void writeFileTrueFalse(ofstream &openFile, string condition, bool tf){
    openFile<<condition<<": ";
    if(tf) openFile<<"true"<<endl;
    else openFile<<"false"<<endl;
    return;
}

void writeFileDashes(ofstream &openFile, int n){
    //write dashed line
    for (int i=0; i<n;++i){
        openFile<<"-";
    }
    openFile<<endl;
}

void writeFileCrd(ofstream &openFile, Crd2d crd){
    //write coordinate to file
    openFile<<setw(15)<<left<<crd.x<<setw(15)<<left<<crd.y<<endl;
}