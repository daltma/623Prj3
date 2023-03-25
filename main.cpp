#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <fstream>
#include <sstream>

//#include "massMatrix.h"
#include "functions.h"
#include "read_gri_p3.h"
#include "processMesh.h"

using namespace std;

int main() {
    cout << "Input P: " << endl;
    int p;
    cin >> p;
    string fileName = "/home/daltma/Desktop/623/prj3/multiQ3.gri";

    int nnode = getnnode(fileName);
    int nbedge = getnbedge(fileName);
    int nelem = getnelem(fileName);     
    int niedge = 0.5 * (3 * nelem - nbedge); // number of interior edges

    vector<vector<double>> nodes = getnodes(fileName);
    vector<vector<double>> elem = getelement(fileName);
    vector<vector<double>> bounds = getboundaries(fileName);

// Process Structures for Curved Mesh
    int ncelem = 0;
    for (int i = 0; i < nelem; i++){
        if (elem[i].size() > 3){
            ncelem += 1;
        }
    }

    vector<vector<double>> curveNodes(nelem,vector<double>(10));

    ofstream outputFile;

    for (int i = 0; i < nelem; i++){
        if (elem[i].size() > 3){
            for (int j = 0; j < 10;j++){
                curveNodes[i][j] = elem[i][j];
            }
            elem[i] = {elem[i][0],elem[i][3],elem[i][9]};
        }
        else {
            for (int j = 0; j < 3;j++){
                curveNodes[i][j] = elem[i][j];
                curveNodes[i][j+7] = 0;
            }
        }
    } 
//
    outputFile.close();

    vector<vector<double>> interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
    vector<vector<double>> I2E = genI2E(niedge, nelem, interiorFaces, elem);
    vector<vector<double>> B2E = genB2E(nbedge, nelem, bounds, elem);
    
    vector<vector<double>> In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    vector<vector<double>> Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    vector<double> Area = genArea(nelem, elem, nodes);

    vector<double> speed(nelem);

    int prank = (p+1)*(p+2)/2;
    vector<vector<double>> U(4, vector<double>(nelem*prank)); // initialize u
    for (int i = 0; i < nelem*prank; i++){
        U[0][i]=1;
        U[1][i]=.25*cos(8*3.14159/180);
        U[2][i]=.25*sin(8*3.14159/180);
        U[3][i]=1/(.4*1.4)+.25*.25/2;
    }

    double CFL = 1;
/*
    vector<vector<double>> resi = calcRes(p, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, U,elem, speed,curveNodes);

    double L1 = l1norm(resi, prank, nelem);
    cout << "Residual L1 Norm: "<<L1 << endl;
*/
    rk4out results;
    results = dg(p, 100000, CFL, Area, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, U, elem, speed, curveNodes);
    double sum;

    outputFile.open("/home/daltma/Desktop/623/prj3/p0_multiQ3_residual.txt");
    for (int i = 0; i < results.residuals.size(); i++){
        outputFile << results.residuals[i];
        outputFile<<endl;
    }
    outputFile.close();

    outputFile.open("/home/daltma/Desktop/623/prj3/p0_multiQ3_state.txt");
    for (int i=0;i<prank*nelem;i++){ 
        for (int k = 0; k<4;k++){
            outputFile << results.U[k][i] <<" ";
        }
        outputFile << endl;
    }
    outputFile.close();
}