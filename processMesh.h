//
//  processMesh.h
//  AEROSP-523-Project-1
//
//  Created by Jake Yeaman on 1/21/23.
//

#ifndef processMesh_h
#define processMesh_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <cmath>

using namespace std;

// OBTAINING SIZES OF NODE, BOUNDARY, AND ELEMENT MATRICES FOR PREALLOCATION PURPOSES
void getMatSizes(string fileName, int &nnode, int &nbedge, int &nelem){
	// Open the file for reading
	ifstream file(fileName);
	
	nnode = 0;
	nbedge = 0;
	nelem = 0;

	// If there is a file in the directory matching the input parameter "fileName", open it
	if (file.is_open()) {
		// String to store each line of the file
		string currLine;
		
		// Parse each line of the file
		while (getline(file, currLine)) {
			
			// COUNT THE NUMBER OF NODES IN MESH BY COUNTING THE NUMBER OF LINES STARTING WITH "GRID"
			if (currLine.find("GRID") == 0) {
				nnode++;
			}
			
			// COUNT THE NUMBER OF BOUNDARY EDGES IN MESH BY COUNTING THE NUMBER OF LINES STARTING WITH "CBAR"
			else if (currLine.find("CBAR") == 0) {
				nbedge++;
			}
			// COUNT THE NUMBER OF ELEMENTS IN MESH BY COUNTING THE NUMBER OF LINES STARTING WITH "CTRIA3"
			else if (currLine.find("CTRIA3") == 0) {
				nelem++;
			}
		}
		// Close the file
		file.close();
	}
	// If there is not a file in the directory matching the input parameter "fileName", cout an error statement
	else {
		cout << "Error: Unable to open file " << fileName << "\n";
	}
	
}

// READ GMSH .BDF OUTPUT FILE INTO NODE, BOUNDARY, AND ELEMENT MATRICES
// User input a Gmsh ".bdf" file type
/// PASS IN THE PARAMETERS WE NEED BY REFERENCE /////
void readGmshFile(string fileName,int nnode,int nbedge,int nelem, vector<vector<double>> &nodes,vector<vector<double>> &elem,vector<vector<double>> &bounds) {
    // Open the file for reading
    ifstream file(fileName);
	
	
	int discard; // Dummy variable for parameters we do not care about

    // If there is a file in the directory matching the input parameter "fileName", open it
    if (file.is_open()) {
        // String to store each line of the file
        string currLine;
		
		// INITIALIZE DATA STRUCTURE SIZES

		string type;
		
		while (getline(file, currLine)) {
			
			// READ IN THE NODE DATA (Coordinates and row indices represent node numbers)
			if (currLine.find("GRID") == 0) {
				
				size_t ID;
				stringstream ss(currLine);
				ss >> type >> ID >> discard;
				string xStr = currLine.substr(24, 8);
				string yStr = currLine.substr(32, 8);
				nodes[ID-1][0] = stod(xStr);
				nodes[ID-1][1] = stod(yStr);
			}
			
			// READ IN THE BOUNDARY
			else if (currLine.find("CBAR") == 0) {
				size_t ID;

				stringstream ss(currLine);
				ss >> type >> ID;
				ID = ID - 1;
				ss >> bounds[ID][2] >> bounds[ID][0] >> bounds[ID][1];// CBAR >> bLineID >> boundaryGroup >> node1 >> node2
			}
			// COUNT THE NUMBER OF ELEMENTS IN MESH BY COUNTING THE NUMBER OF LINES STARTING WITH "CTRIA3"
			else if (currLine.find("CTRIA3") == 0) {
				size_t ID;
				stringstream ss(currLine);
				ss >> type >> ID >> discard;
				ID = ID - nbedge - 1;
				ss >> elem[ID][0] >> elem[ID][1] >> elem[ID][2];
			}
			else if (currLine.find("ENDDATA") == 0) {
				break;
			}
		}
		
        // Close the file
        file.close();
    }
    // If there is not a file in the directory matching the input parameter "fileName", cout an error statement
    else {
        cout << "Error: Unable to open file " << fileName << "\n";
    }
}

// CREATE A VECTOR WHICH CONTAINS THE AMOUNT OF EDGES IN EACH BOUNDARY GROUP
vector<int> boundSizes(vector<vector<double>> &bounds,int nbedge){
	
	unordered_set<int> unique_values;
	
	for(int i=0;i < nbedge; i++){
		unique_values.insert(bounds[i][2]);
	}
	
	size_t nBGroup = unique_values.size();
	vector<int> nBFace(nBGroup,0); // Vector which reports the # of bound faces in each bound group. Vector index represents bGroup
	
	for(int i=0;i < nbedge; i++){
		nBFace[bounds[i][2] - 1]++; // incrementing group index by one
	}
	
	return nBFace;
}

// FROM THE BASIC MESH DATA STRUCTURES (NODES, ELEMENTS, AND BOUNDARIES) GENERATE A MESH OUTPUT FILE IN .GRI FORMAT
void generateGri(string fileName, int nnode,int nbedge,int nelem, vector<vector<double>> nodes, vector<vector<double>> elem, vector<vector<double>> bounds){
	ofstream file(fileName);
	if (file.is_open()) {
			
		file << nnode << " " << nelem << " 2\n";
		
		// Writing Node Data Section
		for(int i = 0; i < nnode; i++){
			file << nodes[i][0] << " " << nodes[i][1] << "\n";
		}
		
		// Number of Boundary Groups
		vector<int> nBFace = boundSizes(bounds,nbedge);
		int nBGroup = int(nBFace.size());
		
		// Writing Boundary Section
		file << nBGroup << "\n";
		for(int i = 0; i < nBGroup; i++){
			file << nBFace[i] << " 2 " << (i + 1) << "\n"; // file << # bound. faces in i+1th group << 2 = linear nodes << BGroup Title == i+1th index
			
			int startingIndex = 0;
			for(int j = 0; j < nbedge; j++){
				// Find the first occurance of the current boundary group in the bounds matrix
				// Based on the format of the .bdf format, boundaries are given in order of boundary group so this makes sense
				if(bounds[j][2] == (i + 1)){
					startingIndex = j;
					break;
				}
			}
			
			// Writing the boundary node values section for the given boundary group in .bdf file
			for (int j = 0; j < nBFace[i]; j++){
				file << bounds[startingIndex + j][0] << " " << bounds[startingIndex + j][1] << "\n";
			}
		}
		
		// Writing Element Section
		// ASSUMING THERE IS ONLY ONE ELEMENT GROUP, COULD BE WRONG BUT I THINK THIS IS THE CASE
		
		file << nelem << " 1 TriLagrange\n"; // nElem Order Basis
		for (int i = 0; i < nelem; i++){
			file << elem[i][0] << " " <<  elem[i][1] << " " << elem[i][2] << "\n";
		}
		
		// Close the file
		file.close();
		}
	else {
		cout << "Error: Unable to open file.\n";
	}
}

// TASK 2

vector<vector<double>> genInteriorFaceVec(int niedge, int nelem, int nbedge, vector<vector<double>> &bounds, vector<vector<double>> &elem){
	
	vector<vector<double>> interiorFaces;
	int facePairs[3][2] = { {0, 1}, {1, 2}, {2, 0} };
	
	// Creating a set of all nodes that lie on a boundary
	unordered_set<double> bNodeSet;
	for(int iB = 0; iB < nbedge; iB++){
		bNodeSet.insert(bounds[iB][0]);
		bNodeSet.insert(bounds[iB][1]);
	}
	
	// Parsing through all faces and storing those on the interior into the interiorFaces vector
	for(int iElem = 0; iElem < nelem; iElem++){
		for (int iFace = 0; iFace < 3; iFace++){
			// Parsing through each face of the current element
			double faceNode1 = elem[iElem][facePairs[iFace][0]];
			double faceNode2 = elem[iElem][facePairs[iFace][1]];
			// Check to see if both nodes do not lie on a boundary
			if(bNodeSet.find(faceNode1)==bNodeSet.end() || bNodeSet.find(faceNode2)==bNodeSet.end()){
				// Ensure face is not already stored
				bool inVec = false;
				int indexMatchingFace = -1;
				for(int i = 0; i < interiorFaces.size(); i++){
					// Is already accounted for
					if( (faceNode1 == interiorFaces[i][0] && faceNode2 == interiorFaces[i][1]) || (faceNode2 == interiorFaces[i][0] && faceNode1 == interiorFaces[i][1]) ){
						indexMatchingFace = i;
						inVec = true;
					}
				}
				// If boundary face not yet stored, insert it into interiorFaces vector
				if(inVec == false){
					interiorFaces.push_back({faceNode1,faceNode2, double(iElem) + 1});
				}
				else{ // If boundary face is already stored, add the second element that is adjacent to face
					interiorFaces[indexMatchingFace].push_back(iElem + 1);
				}
				// Reset the flag
				inVec = false;
			}
		}
	}
		
	return interiorFaces;
}

vector<vector<double>> genI2E(int niedge, int nelem, vector<vector<double>> &interiorFaces, vector<vector<double>> &elem){
	vector<vector<double>> I2E(niedge, vector<double>(4));
	
	int intEdgeNum = int(interiorFaces.size());
	int counter = 0;
	for(int i = 0; i < intEdgeNum; i++){
		double node1 = interiorFaces[i][0];
		double node2 = interiorFaces[i][1];
		double elemL = interiorFaces[i][2];
		double elemR = interiorFaces[i][3];
		// Calculate local face number for elemL
		double faceL = -1;
		vector<double> elemNodesL({elem[elemL-1][0],elem[elemL-1][1],elem[elemL-1][2]});
		for(int j = 0; j < 3; j++){
			if(elemNodesL[j] != node1 && elemNodesL[j] != node2){
				// faceL = elemNodesL[j];
				faceL = j+1;
				break;
			}
		}
		// Calculate local face number for elemR
		double faceR = -1;
		vector<double> elemNodesR({elem[elemR-1][0],elem[elemR-1][1],elem[elemR-1][2]});
		for(int j = 0; j < 3; j++){
			if(elemNodesR[j] != node1 && elemNodesR[j] != node2){
				// faceR = elemNodesR[j];
				faceR = j+1;
				break;
			}
		}
		
		I2E[i][0] = elemL;
		I2E[i][1] = faceL;
		I2E[i][2] = elemR;
		I2E[i][3] = faceR;
		counter = i;
		//cout<<I2E[i][0]<<" "<<I2E[i][1]<<" "<<I2E[i][2]<<" "<<I2E[i][3]<<" iteration: "<<counter<<"\n";
	}
	
	return I2E;
}

// NEEDS TO VERIFY
vector<vector<double>> genB2E(int nbedge, int nelem, vector<vector<double>> &bounds,vector<vector<double>> &elem){
	vector<vector<double>> B2E(nbedge, vector<double>(3));
	
	for(int i = 0; i < nbedge; i++){
		double n1 = bounds[i][0];
		double n2 = bounds[i][1];
		// Saving bgroup
		B2E[i][2] = bounds[i][2];
		// find element containing nodes
		for(int j = 0; j < nelem; j++){
			
			double val1 = elem[j][0];
			double val2 = elem[j][1];
			double val3 = elem[j][2];
			
			if ((val1 == n1 && val2 == n2) || (val1 == n2 && val2 == n1) || (val2 == n1 && val3 == n2) || (val2 == n2 && val3 == n1) || (val1 == n1 && val3 == n2) || (val1 == n2 && val3 == n1)) {
				
				B2E[i][0] = double(j + 1);
				
				if(val1 != n1 && val1 != n2){
					// B2E[i][1] = elem[j][0];
					B2E[i][1] = 1;
				}
				else if(val2 != n1 && val2 != n2){
					// B2E[i][1] = elem[j][1];
					B2E[i][1] = 2;
				}
				else{
					// B2E[i][1] = elem[j][2];
					B2E[i][1] = 3;
				}
				
				break;
				
			}
		}
	}
	
	return B2E;
}

vector<vector<double>> genIn(int niedge, vector<vector<double>> &interiorFaces, vector<vector<double>> &elem, vector<vector<double>> &nodes, vector<vector<double>> &I2E){
	
	vector<vector<double>> In(niedge, vector<double>(2));
	
	for(int i = 0; i < niedge; i++){
		// Nodes on  boundary
		
		double n1 = interiorFaces[i][0], n2 = interiorFaces[i][1];
		// local face value
		size_t iFaceLocal = size_t(I2E[i][3]) - 1;
		size_t iElem = interiorFaces[i][3];
		double fLocal = elem[iElem - 1][iFaceLocal]; // faceR (the node)
		size_t ifLocal = fLocal - 1; // index of local node
		// Obtain coordinate of local node and boundary node
		double x1 = nodes[n1 - 1][0], x2 = nodes[n2 - 1][0];
		double y1 = nodes[n1 - 1][1], y2 = nodes[n2 - 1][1];
		double xLocal = nodes[ifLocal][0], yLocal = nodes[ifLocal][1];
		
		double magnitude = sqrt(pow(x2-x1, 2) + pow(y2-y1,2));
		
		double xComponent = (1/magnitude)*(x2-x1);
		double yComponent = (1/magnitude)*(y2-y1);

		// TODO: Now compute the normal by rotating the vector pointing from n1 to n2 away from localNode and normalize
		
		// Seeing where xLocal and yLocal lie relative to the extention of the face to determine which way to rotate the face vector
		double slope = (y2-y1)/(x2-x1);
		double xNorm = 0, yNorm = 0;
		
		if(isinf(slope)){ // if vertical line compare xLocal and x1
			double xLine = x1;
			if(xLine < xLocal){
				if(y2 > y1){ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
				else{ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
			}
			else{ // xLine > xLocal
				if(y2 > y1){ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
				else{ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
			}
		}
		else{
			double yLine = slope*(xLocal - x1) + y1;
			// Rotation counter clockwise implemented as follows, need to figure out how to check to see if we must instead rotate clockwise.
			if(yLine < yLocal){
				if(x2 > x1){ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
				else{ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
			}
			else{ // yLine > yLocal
				if(x2 > x1){ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
				else{ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
			}
		}
		In[i] = vector<double> {xNorm,yNorm};
	}
	
	return In;
}

vector<vector<double>> genBn(int nbedge, vector<vector<double>> &bounds, vector<vector<double>> &nodes, vector<vector<double>> &elem, vector<vector<double>> &B2E){
	
	vector<vector<double>> Bn(nbedge, vector<double>(2));
	
	for(int i = 0; i < nbedge; i++){
		// Nodes on  boundary
		double n1 = bounds[i][0], n2 = bounds[i][1];
		// local face value
		size_t iFaceLocal = size_t(B2E[i][1]) - 1;
		size_t iElem = B2E[i][0];

		double fLocal = elem[iElem - 1][iFaceLocal]; // faceR (the node)
		size_t ifLocal = fLocal - 1; // index of local node
		// Obtain coordinate of local node and boundary node
		double x1 = nodes[n1 - 1][0], x2 = nodes[n2 - 1][0];
		double y1 = nodes[n1 - 1][1], y2 = nodes[n2 - 1][1];
		double xLocal = nodes[ifLocal][0], yLocal = nodes[ifLocal][1];

		double magnitude = sqrt(pow(x2-x1, 2) + pow(y2-y1,2));
		
		double xComponent = (1/magnitude)*(x2-x1);
		double yComponent = (1/magnitude)*(y2-y1);
		
		// TODO: Now compute the normal by rotating the vector pointing from n1 to n2 away from localNode and normalize
		
		// Seeing where xLocal and yLocal lie relative to the extention of the face to determine which way to rotate the face vector
		double slope = (y2-y1)/(x2-x1);
		double xNorm = 0, yNorm = 0;
		
		if(isinf(slope)){ // if vertical line compare xLocal and x1
			double xLine = x1;
			if(xLine < xLocal){
				if(y2 > y1){ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
				else{ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
			}
			else{ // xLine > xLocal
				if(y2 > y1){ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
				else{ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
			}
		}
		else{
			double yLine = slope*(xLocal - x1) + y1;
			// Rotation counter clockwise implemented as follows, need to figure out how to check to see if we must instead rotate clockwise.
			if(yLine < yLocal){
				if(x2 > x1){ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
				else{ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
			}
			else{ // yLine > yLocal
				if(x2 > x1){ // CW
					xNorm = yComponent;
					yNorm = -xComponent;
				}
				else{ // CCW
					xNorm = -yComponent;
					yNorm = xComponent;
				}
			}
		}
		
		Bn[i] = vector<double> {xNorm,yNorm};

	}
	
	return Bn;
}

vector<double> genArea(int nelem, vector<vector<double>> &elem, vector<vector<double>> &nodes){
	vector<double> AREA(nelem);
	
	for(int i = 0; i < nelem; i++){
		
		double n1, n2, n3;
		double x1, y1, x2, y2, x3, y3;
		double a, b, c, s;
		
		n1 = elem[i][0];
		n2 = elem[i][1];
		n3 = elem[i][2];
		
		x1 = nodes[n1 - 1][0];
		x2 = nodes[n2 - 1][0];
		x3 = nodes[n3 - 1][0];
		
		y1 = nodes[n1 - 1][1];
		y2 = nodes[n2 - 1][1];
		y3 = nodes[n3 - 1][1];
		
		// Calculate the lengths of the three sides
		a = sqrt(pow(x2 - x3, 2) + pow(y2 - y3, 2));
		b = sqrt(pow(x1 - x3, 2) + pow(y1 - y3, 2));
		c = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
		
		// Calculate the semi-perimeter of the triangle
		s = (a + b + c) / 2;
		
		AREA[i] = sqrt(s * (s - a) * (s - b) * (s - c));;
		
	}
	
	return AREA;
}

// VERIFICATION
// Calculates the error on each element for mesh verification
vector<double> verification(int nelem, vector<vector<double>> &I2E, vector<vector<double>> &B2E, vector<vector<double>> &In, vector<vector<double>> &Bn, vector<vector<double>> &nodes, vector<vector<double>> &elem, vector<vector<double>> &bounds){
	
	vector<vector<double>> error(nelem, vector<double>(2));
	vector<double> err_mag(nelem);
	int I2E_size = int(I2E.size());
	int B2E_size = int(B2E.size());

	// loop through interior edges
	for (int i = 0; i < I2E_size; i++){
		int elemL = I2E[i][0] - 1;
		int elemR = I2E[i][2] - 1;
		int face = I2E[i][1];
		double length = -1;
		vector<double> coord1(2);
		vector<double> coord2(2);
		

		// calculate length of face
		if (face == 1){
			int node1 = elem[elemL][1] - 1;
			int node2 = elem[elemL][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 2){
			int node1 = elem[elemL][0] - 1;
			int node2 = elem[elemL][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 3){
			int node1 = elem[elemL][0] - 1;
			int node2 = elem[elemL][1] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		
		// compute error of interior edges with length and normal vectors
		error[elemL][0] -= In[i][0]*length;
		error[elemL][1] -= In[i][1]*length;
		error[elemR][0] += In[i][0]*length;
		error[elemR][1] += In[i][1]*length;
		
	}
	// loop through boundary edges
	for (int i = 0; i < B2E_size; i++){
		int Belem = B2E[i][0] - 1;
		int face = B2E[i][1];
		// int bgroup = B2E[i][2];
		double length = -1;
		vector<double> coord1(2);
		vector<double> coord2(2);

		// calculate length of face
		if (face == 1){
			int node1 = elem[Belem][1] - 1;
			int node2 = elem[Belem][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 2){
			int node1 = elem[Belem][0] - 1;
			int node2 = elem[Belem][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 3){
			int node1 = elem[Belem][0] - 1;
			int node2 = elem[Belem][1] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		
		// compute error of boundary edges with length and normal vectors
		error[Belem][0] += Bn[i][0]*length;
		error[Belem][1] += Bn[i][1]*length;

	}
	// compute magnitude of error for each element
	for (int i = 0; i < nelem; i++){
		err_mag[i] = sqrt(pow(error[i][0],2) + pow(error[i][1],2));
	}
	return err_mag;
}



#endif /* processMesh_h */
