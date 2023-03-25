//
//  read_gri_p3.h
//  p3-AEROSP-623
//
//  Created by Jake Yeaman on 3/21/23.
//

#ifndef read_gri_p3_h
#define read_gri_p3_h

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

vector<vector<double>> getnodes(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
    int trash;
    inputFile >> trash;
//        int b4;
//    inputFile >> b4; // farfield boundary
//        int b1;
//    inputFile >> b1; // slat boundary
//        int b2;
//    inputFile >> b2; // main boundary
//        int b3;
//    inputFile >> b3; // flap boundary

    // Create a vector of vectors to store nodes
    vector<vector<double>> nodes(numnodes, vector<double>(2));

    // Read in nodes from the file
    for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < 2; j++) {
            inputFile >> nodes[i][j];
        }
    }

    // Close the input file
    inputFile.close();

return nodes;

}



vector<vector<double>> getelement(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
    
    int trash;
    inputFile >> trash;
    
    int skip = numnodes*2;

    
    // skip the values that I dont want
     vector<double> skipvec(skip);

    for (int i = 0; i < skip; i++) {
     inputFile >> skipvec[i];
    
    }
    
    // skip boundaries
    
    int discard;
    string discardString;
    inputFile >> discard; // number of boundaries
    
    int b1, b2, b3, b4;
    
    inputFile >> b4 >> discard >> discardString;
    
    // skip boundary 4
    for(int i = 0; i < b4; i++){
        inputFile >> discard >> discard;
    }
    
    inputFile >> b1 >> discard >> discardString;
    
    // skip boundary 1
    for(int i = 0; i < b1; i++){
        inputFile >> discard >> discard;
    }
    
    inputFile >> b2 >> discard >> discardString;
    
    // skip boundary 2
    for(int i = 0; i < b2; i++){
        inputFile >> discard >> discard;
    }
    
    inputFile >> b3 >> discard >> discardString;
    
    // skip boundary 3
    for(int i = 0; i < b3; i++){
        inputFile >> discard >> discard;
    }
    
    // elements consisting of three nodes
    int n3Elem;
    inputFile >> n3Elem >> discard >> discardString;
    
    vector<double> currElem(3);
    vector<vector<double>> elem;
    
    for(int i = 0; i<n3Elem; i++){
        inputFile >> currElem[0] >> currElem[1] >> currElem[2];
        elem.push_back(currElem);
    }
    
    // elements consisting of three nodes
    int n10Elem;
    inputFile >> n10Elem >> discard >> discardString;
    
    vector<double> currElem2(10);
    
    for(int i = 0; i<n10Elem; i++){
        inputFile >> currElem2[0] >> currElem2[1] >> currElem2[2] >> currElem2[3] >> currElem2[4] >> currElem2[5] >> currElem2[6] >> currElem2[7] >> currElem2[8] >> currElem2[9];
        elem.push_back(currElem2);
    }
    


//
//    // Read in elem from the file
//    for (int i = 0; i < numelem; i++) {
//        for (int j = 0; j < 3; j++) {
//            inputFile >> elem[i][j];
//        }
//    }

    // Close the input file
    inputFile.close();

return elem;

}

vector<vector<double>> getboundaries(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
        int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
    int trash;
    inputFile >> trash;
    
//        int b4;
//    inputFile >> b4; // farfield boundary
//        int b1;
//    inputFile >> b1; // slat boundary
//        int b2;
//    inputFile >> b2; // main boundary
//        int b3;
//    inputFile >> b3; // flap boundary

    int skip=numnodes*2;

   
    // skip the values that I dont want
     vector<double> skipvec(skip);

    for (int i = 0; i < skip; i++) {
     inputFile >> skipvec[i];

    }
//
//    // Create a vector of vectors to store bounds
//    vector<vector<double>> bounds((b1+b2+b3+b4), vector<double>(3));
    
    int discard;
    string discardString;
    inputFile >> discard; // number of boundaries
    vector<vector<double>> bounds;
    vector<double> currBounds(3);
    
    int b1, b2, b3, b4;
    
    inputFile >> b4 >> discard >> discardString;
    
    // skip boundary 4
    for(int i = 0; i < b4; i++){
        inputFile >> currBounds[0] >> currBounds[1];
        currBounds[2] = 4;
        bounds.push_back(currBounds);
    }
    
    inputFile >> b1 >> discard >> discardString;
    
    // skip boundary 1
    for(int i = 0; i < b1; i++){
        inputFile >> currBounds[0] >> currBounds[1];
        currBounds[2] = 1;
        bounds.push_back(currBounds);
    }
    
    inputFile >> b2 >> discard >> discardString;
    
    // skip boundary 2
    for(int i = 0; i < b2; i++){
        inputFile >> currBounds[0] >> currBounds[1];
        currBounds[2] = 2;
        bounds.push_back(currBounds);
    }
    
    inputFile >> b3 >> discard >> discardString;
    
    // skip boundary 3
    for(int i = 0; i < b3; i++){
        inputFile >> currBounds[0] >> currBounds[1];
        currBounds[2] = 3;
        bounds.push_back(currBounds);
    }
    
    
    // Close the input file
    inputFile.close();

return bounds;

}

int getnnode(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int nnode;
    inputFile >> nnode; // nodes
    
    inputFile.close();

return nnode;

}

int getnelem(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int nnode;
    inputFile >> nnode; // nodes
        int nelem;
    inputFile >> nelem; // elemnts
        
    inputFile.close();

return nelem;

}

int getnbedge(string fileName){

    // Open the input file
    ifstream inputFile(fileName);

    // Read the first line
    int numnodes;
    inputFile >> numnodes; // nodes
        int numelem;
    inputFile >> numelem; // elemnts
    
    int trash;
    inputFile >> trash;
    
    int skip = numnodes*2;

    
    // skip the values that I dont want
     vector<double> skipvec(skip);

    for (int i = 0; i < skip; i++) {
     inputFile >> skipvec[i];
    
    }
    
    // skip boundaries
    
    int discard;
    string discardString;
    inputFile >> discard; // number of boundaries
    
    int b1, b2, b3, b4;
    
    inputFile >> b4 >> discard >> discardString;
    
    // skip boundary 4
    for(int i = 0; i < b4; i++){
        inputFile >> discard >> discard;
    }
    
    inputFile >> b1 >> discard >> discardString;
    
    // skip boundary 1
    for(int i = 0; i < b1; i++){
        inputFile >> discard >> discard;
    }
    
    inputFile >> b2 >> discard >> discardString;
    
    // skip boundary 2
    for(int i = 0; i < b2; i++){
        inputFile >> discard >> discard;
    }
    
    inputFile >> b3 >> discard >> discardString;
    
    // skip boundary 3
    for(int i = 0; i < b3; i++){
        inputFile >> discard >> discard;
    }

    int nbedge = b1+b2+b3+b4;
    
    inputFile.close();

return nbedge;

}


#endif /* read_gri_p3_h */
