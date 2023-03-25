#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <memory>
#include <array>

#include "quadatures.h"
#include "fluxes.h"
#include "shape.h"
#include "invMatrix.h"

using namespace std;

struct outNorm{
    vector<double> nv;
    double jEdge;
};
struct rk4stage{
    vector<vector<double>> F;
    vector<vector<double>> U;
};
struct rk4out{
    vector<vector<double>> U;
    vector<double> residuals;
};

// converts 1D reference space to 2D coordinates, unrolled array
// [ x1 x2 x3 y1 y2 y3 ...]
unique_ptr<array<double,2>> ref2glob(int edge, vector<double> &s, const int direct){
    int size = s.size();

    auto xy = make_unique<array<double, 2>>(array<double, 2>{});

    if (direct == 1){
        for (int i = 0; i < size; i++){
            s[i] = 1 - s[i];
        }
    }

    switch (edge){
        case 1:
        for (int i = 0; i < size; i++){
            (*xy)[i] = 1-s[i];
            (*xy)[i+size] = s[i];
        }
        break;
        case 2:
        for (int i = 0; i < size; i++){
            (*xy)[i] = 0;
            (*xy)[i+size] = 1-s[i];
        }
        break;
        case 3:
        for (int i = 0; i < size; i++){
            (*xy)[i] = s[i];
            (*xy)[i+size] = 0;
        }
        break;
    }  
    return xy;
}
vector<double> timestep(int nelem, double CFL, vector<double> Area, vector<double> speed){
    vector<double> time(nelem);

    for (int i = 0; i < nelem; i++){
        time[i] = 2*CFL*Area[i]/speed[i];
    }
    return time;
}
vector<double> glagrange(vector<double> xn, int j, vector<double> x) {
    int n = xn.size();
    vector<int> nj;
    for(int i = 0; i < n; i++) {
        if(i != j-1) {
            nj.push_back(i);
        }
    }
    double den = 1.0;
    for(int i = 0; i < nj.size(); i++) {
        den *= xn[j-1] - xn[nj[i]];
    }
    vector<double> gphi(x.size(), 0.0);
    if(nj.size() == 0) {
        return gphi;
    } else if(nj.size() == 1) {
        for(int i = 0; i < x.size(); i++) {
            gphi[i] = 1.0;
        }
    } else {
        for(int k = 0; k < x.size(); k++) {
            double num = 0.0;
            for(int i = 0; i < nj.size(); i++) {
                vector<int> nij;
                for(int j = 0; j < nj.size(); j++) {
                    if(j != i) {
                        nij.push_back(nj[j]);
                    }
                }
                vector<double> xnij;
                for(int j = 0; j < nij.size(); j++) {
                    xnij.push_back(xn[nij[j]]);
                }
                double prod = 1.0;
                for(int j = 0; j < xnij.size(); j++) {
                    prod *= x[k] - xnij[j];
                }
                num += prod;
            }
            gphi[k] = num / den;
        }
    }
    return gphi;
}
vector<vector<double>> gbasis(vector<double> xn, vector<double> x){
    int n = x.size();
    int order = xn.size() - 1;
    vector<vector<double>> gPhi(order+1,vector<double>(n));
    for (int i = 0; i < order+1; i++){
        vector<double> B = glagrange(xn,i+1,x);
        gPhi[i] = B;
    }
    return gPhi;
}
outNorm normal(vector<vector<double>> xedge, vector<double> quad,int q, int q_pt){
    outNorm output;
    
    vector<double> n(2);

    int Nq = q + 1;

    double dx = 0;
    double dy = 0;

    vector<double> xref{0,0.333333333333333,0.666666666666667,1};

    vector<vector<double>> gphi = gbasis(xref,quad);
    
    for (int j = 0; j < Nq; j++)  {  
        dx += xedge[j][0]*gphi[j][q_pt];
        dy += xedge[j][1]*gphi[j][q_pt];
    }

    double magdx = sqrt(dx*dx + dy*dy);
    
    double tx = dx/magdx;
    double ty = dy/magdx;

    n = {ty, -tx};
    
    output.nv = n;
    output.jEdge = magdx;

    return output;
}
void jacobian(double **quad, vector<vector<double>> coord, vector<vector<double>> &invJ, int q, int Nq, double* detJ){

    double* xref = new double[2];
    int gphi_size = (q+1)*(q+2);
   
    double* xy = *quad;

    double* gphi = new double[gphi_size];

    for (int j = 0; j < Nq; j++){
        xref[0] = xy[j];
        xref[1] = xy[j+Nq];
        
        gradientL(xref, q, &gphi);
        vector<vector<double>> matJ(2, vector<double>(2));
        double det=0;
        
        matJ[0][0] = 0;
        matJ[0][1] = 0;
        matJ[1][0] = 0;
        matJ[1][1] = 0;

        for (int i = 0; i < gphi_size/2; i++){
            matJ[0][0] = matJ[0][0] + coord[i][0]*gphi[i];
            matJ[0][1] = matJ[0][1] + coord[i][0]*gphi[i+(gphi_size/2)];
            matJ[1][0] = matJ[1][0] + coord[i][1]*gphi[i];
            matJ[1][1] = matJ[1][1] + coord[i][1]*gphi[i+(gphi_size/2)]; 
        }

        det = matJ[0][0]*matJ[1][1] - matJ[0][1]*matJ[1][0];

        detJ[j] = det;
        invJ[0][j] = matJ[1][1]/det;
        invJ[1][j] = -matJ[0][1]/det;
        invJ[2][j] = -matJ[1][0]/det;
        invJ[3][j] = matJ[0][0]/det;
    }
    delete[] gphi;
    delete[] xref;
}
void edgeRes(int nelem, int p, vector<vector<double>> interiorFaces, vector<vector<double>> bounds, vector<vector<double>> I2E, vector<vector<double>> B2E, vector<vector<double>> Bn, vector<vector<double>> In, vector<vector<double>> nodes, vector<vector<double>> &residual, vector<vector<double>> U, vector<double> &speed, vector<vector<double>> curveNodes, vector<vector<double>> elem){
    int niedge = interiorFaces.size(); // initialize number of interior edges
    int nbedge = bounds.size(); // initialize number of 
    quadOut out1D;
    quadOut out2D;
    quadOut out3D;
    quadOut out4D;

    int q1 = 1; // geometry order of interior edge
    int o = 2*p + 1;
    
    out1D = quadature(o, 1);

    int N1 = out1D.Nq;
    vector<double> quad1(N1);
    vector<double> w1(N1);

    for (int i = 0; i < N1; i++){
        double xdata = out1D.x[i];
        double wdata =out1D.wq[i];
        w1[i] = wdata;
        quad1[i] = xdata;
    }

    int prank = (p+1)*(p+2)/2;
    double *pphi =  new double[prank];
    double *pphi1 =  new double[prank];
    double *xref = new double[2];
    double *xref2 = new double[2];
    double ephi[3][N1*prank];
    double ephiCC[3][N1*prank];

    structFlux Flux;

    for (int e = 1; e < 4; e++){

        auto xy= ref2glob(e, quad1, 0);
        auto xyCC= ref2glob(e, quad1, 1);

        for (int j = 0; j < N1; j++){
            xref[0] = (*xy)[j];
            xref[1] = (*xy)[j+N1];

            shapeL(xref, p, &pphi);
            for (int k = 0; k < prank; k++){
                ephi[e-1][k+(prank*j)] = pphi[k];
            }

            xref2[0] = (*xyCC)[j];
            xref2[1] = (*xyCC)[j+prank];

            shapeL(xref2, p, &pphi1);
            for (int k = 0; k < prank; k++){
                ephiCC[e-1][k+(prank*j)] = pphi1[k];
            }
        }
    }

    // Loop over Interior Edges
    for (int ie = 0; ie < niedge; ie++){
        vector<double> n(2);
        int elemL = I2E[ie][0] - 1;
        int elemR = I2E[ie][2] - 1;
        int edgeL = I2E[ie][1] - 1;
        int edgeR = I2E[ie][3] - 1;
        int node1 = interiorFaces[ie][0] - 1;
        int node2 = interiorFaces[ie][1] - 1;

        // Compute edge length for residual integration
        double length = sqrt(pow(nodes[node1][0]-nodes[node2][0],2) + pow(nodes[node1][1]-nodes[node2][1],2));

        n[0] = In[ie][0];
        n[1] = In[ie][1];
       
        // Loop over 1D Quadature points
        for (int i = 0; i < N1; i++){
            vector<double> uL = {0,0,0,0};
            vector<double> uR = {0,0,0,0};
            // Loop over entries of state vector
            for (int k = 0; k < 4; k++){
                // Loop over basis functions
                for (int j = 0; j < prank; j++){
                    // Calculate left and right states

                    uL[k] = uL[k] + ephi[edgeL][j+(prank*i)]*U[k][j+(elemL*prank)];
                    uR[k] = uR[k] + ephiCC[edgeR][j+(prank*i)]*U[k][j+(elemR*prank)];
                }
            }

            Flux = roe(uL,uR,1.4,n);
            
            speed[elemL] = Flux.smag * length;
            speed[elemR] = Flux.smag * length;

            // Loop over entries of residual/state vector
            for (int k = 0; k < 4; k++){ 
                // Loop over basis functions 
                for (int j = 0; j < prank; j++){
                    // Add to left, subtract from right
                    residual[k][j+(elemL*prank)] = residual[k][j+(prank*elemL)] + ephi[edgeL][j+(i*prank)]*w1[i]*length*Flux.F[k];
                    residual[k][j+(elemR*prank)] = residual[k][j+(prank*elemR)] - ephi[edgeR][j+(i*prank)]*w1[i]*length*Flux.F[k];
                }
            }

        }
    } // end interior edge loop

    int q2 = 3; // geometry order of curved boundary edge (change for curved meshes)
    int o1; 
    int q3 =1; // geometry order of farfield edge
    int geomO;

    if(q2 == 1){
        o1 = 2*p + 1;
        o = 2*p + 1;
        geomO = 1;
    }
    else if(q2 > 1){
        o1 = 2*p + 1;
        o = 2*p + 1 + (q2 - 1);
        geomO = 2;
    }

    out2D = quadature(o, 1);
    out3D = quadature(o1, 1);

    int N2 = out2D.Nq; //curved edge 1d quad points
    vector<double> quad2(N2); // curved edge 1d quad points
    vector<double> w2(N2); // curved edge 1d weights

    int N3 = out3D.Nq;
    vector<double> quad3(N3); // linear edge 1d quad points
    vector<double> w3(N3);

    for (int i = 0; i < N2; i++){
        double xdata = out2D.x[i];
        w2[i]=out2D.wq[i];
        quad2[i] = xdata;
    }
    for (int i = 0; i < N3; i++){
        double xdata1 = out3D.x[i];
        w3[i]=out3D.wq[i];
        quad3[i] = xdata1;
    }

    double ephi2[3][N2*prank];
    vector<vector<double>> ephi3(3,vector<double>(N3*prank));
    double egphi2[3][N2*prank*2];
    double *pgphi = new double[prank*2];
    
    outNorm normv;
    for (int e = 1; e < 4; e++){
        double *xy1= ref2glob(e, quad2, 0);
        double *xy2= ref2glob(e, quad3, 0);

        // curve edge orientation
        for (int j = 0; j < N2; j++){
            xref[0] = xy1[j];
            xref[1] = xy1[j+N2];
            shapeL(xref, p, &pphi);
            gradientL(xref,p,&pgphi);
            for (int k = 0; k < prank; k++){
                ephi2[e-1][k+(prank*j)] = pphi[k]; 
            }
            for (int k = 0; k < 2*prank; k++){
                egphi2[e-1][k+(2*prank*j)] = pgphi[k];
            } 
        }
        // farfield edge orientation
        for (int j = 0; j < N3; j++){
            xref2[0] = xy2[j];
            xref2[1] = xy2[j+N3];
            shapeL(xref2, p, &pphi);
            for (int k = 0; k < prank; k++){
                ephi3[e-1][k+(prank*j)] = pphi[k];
            }
        }
    }
    // Loop over boundary edges
    for (int be = 0; be < nbedge; be++){
        vector<double> n(2);
        double jEdge;

        int elems = B2E[be][0] - 1;
        int edge= B2E[be][1] - 1;
        int group = B2E[be][2];

        int nnode;
        int elemNode;
        
        if (group == 4){
            nnode = 2;
            elemNode = 3;
        }
        else {
            nnode = 2;
            elemNode = (q2+1)*(q2+2)/2;
        }
        vector<vector<double>> elemXY(elemNode, vector<double>(2)); 

        double length;
        // Store node indeces for boundary edge

        vector<vector<double>> invJ(4,vector<double>(N2));
        vector<vector<double>> coord(nnode, vector<double>(2));
        vector<vector<double>> xedge(4,vector<double>(2));

        if (group == 4){
            for (int i = 0; i < nnode; i++){
                coord[i][0] = nodes[bounds[be][i]-1][0];
                coord[i][1] = nodes[bounds[be][i]-1][1];
            }
            for (int i = 0; i < elemNode; i++){
                elemXY[i][0] = nodes[elem[elems][i]-1][0];
                elemXY[i][1] = nodes[elem[elems][i]-1][1];
            }
        }
        else {
            for (int i = 0; i < nnode; i++){
                //coord[i][0] = nodes[bounds[be][i]-1][0];
                //coord[i][1] = nodes[bounds[be][i]-1][1];
            }
            for (int i = 0; i < elemNode; i++){
                elemXY[i][0] = nodes[curveNodes[elems][i]-1][0];
                elemXY[i][1] = nodes[curveNodes[elems][i]-1][1];
            }
            if (edge == 0){
                for (int j = 0; j < 2; j++){
                    xedge[0][j] = nodes[curveNodes[elems][3]-1][j];
                    xedge[1][j] = nodes[curveNodes[elems][6]-1][j];
                    xedge[2][j] = nodes[curveNodes[elems][8]-1][j];
                    xedge[3][j] = nodes[curveNodes[elems][9]-1][j];
                }
            }
            else if (edge == 1){
                for (int j = 0; j < 2; j++){
                    xedge[0][j] = nodes[curveNodes[elems][9]-1][j];
                    xedge[1][j] = nodes[curveNodes[elems][7]-1][j];
                    xedge[2][j] = nodes[curveNodes[elems][4]-1][j];
                    xedge[3][j] = nodes[curveNodes[elems][0]-1][j];
                }
            }
            else if (edge == 2){
                for (int j = 0; j < 2; j++){
                    xedge[0][j] = nodes[curveNodes[elems][0]-1][j];
                    xedge[1][j] = nodes[curveNodes[elems][1]-1][j];
                    xedge[2][j] = nodes[curveNodes[elems][2]-1][j];
                    xedge[3][j] = nodes[curveNodes[elems][3]-1][j];
                }
            }
        }
        vector<double> uR = {0,0,0,0};
        
        // If farfield (not curved), the normal and jEdge are the same as linear case
        if (bounds[be][2] == 4){
            n[0] = Bn[be][0];
            n[1] = Bn[be][1];
            
            int node1 = bounds[be][0] - 1;
            int node2 = bounds[be][1] - 1;
            jEdge = sqrt(pow(nodes[node1][0]-nodes[node2][0],2) + pow(nodes[node1][1]-nodes[node2][1],2));
            length =  jEdge;

            // Loop over Quadature points 
            for (int i = 0; i < N3; i++){
                vector<double> uL = {0,0,0,0};
                // Loop over entries in state vector
                for (int k = 0; k < 4; k++){
                    // Loop over basis functions
                    for (int j = 0; j < prank; j++){
                        // Compute left state
                        uL[k] = uL[k] + ephi3[edge][j+(prank*i)]*U[k][j+(elems*prank)];
                    }
                    // set right state to free stream
                    uR[0]=1;
                    uR[1]=.25*cos(8*3.14159/180);
                    uR[2]=.25*sin(8*3.14159/180);
                    uR[3]=1/(.4*1.4)+.25*.25/2;
                }

                Flux = wallFlux(uL,n,1.4);
                //Flux = roe(uL,uR,1.4,n); // free stream preservation

                speed[elems] += Flux.smag * length;

                // Loop over entries in residual/state vector
                for (int k = 0; k < 4; k++){
                    // Loop over basis functions 
                    for (int j = 0; j < prank; j++){
                        residual[k][j+(prank*elems)] = residual[k][j+(prank*elems)] - ephi3[edge][j+(i*prank)]*w3[i]*jEdge*Flux.F[k];
                    }
                }
            }
        }    
        // If airfoil boundary (curved), compute normal and jEdge
        else {
            double *xy2d = ref2glob(edge, quad2, 0);
            double *det = jacobian(&xy2d, elemXY, invJ, q2, N2);

            // Loop over 2D Quadature points 

            for (int i = 0; i < N2; i++){
                vector<double> uL = {0,0,0,0};

                vector<double> jMat = {det[i]*invJ[3][i],-det[i]*invJ[1][i],-det[i]*invJ[2][i],det[i]*invJ[0][i]};

                normv = normal(xedge, quad2, q2, i);
                
                int node1 = bounds[be][0] - 1;
                int node2 = bounds[be][1] - 1;
                length = sqrt(pow(nodes[node1][0]-nodes[node2][0],2) + pow(nodes[node1][1]-nodes[node2][1],2));
                double jDet;

                n[0] = normv.nv[0];
                n[1] = normv.nv[1];
                jDet = normv.jEdge;                    

                // Loop over entries in state vector
                for (int k = 0; k < 4; k++){
                    // Loop over basis functions
                    for (int j = 0; j < prank; j++){
                        // Compute left state
                        uL[k] = uL[k] + ephi2[edge][j+(prank*i)]*U[k][j+(elems*prank)]; 
                        
                    }
                    // set right state to free stream
                    uR[0]=1;
                    uR[1]=.25*cos(8*3.14159/180);
                    uR[2]=.25*sin(8*3.14159/180);
                    uR[3]=1/(.4*1.4)+.25*.25/2;
                }

                Flux = wallFlux(uL,n,1.4);
                //Flux = roe(uL,uR,1.4,n); // free stream preservation
                
                speed[elems] += Flux.smag * length;

                // Loop over entries in residual/state vector
                for (int k = 0; k < 4; k++){
                    // Loop over basis functions 
                    for (int j = 0; j < prank; j++){
                        residual[k][j+(prank*elems)] = residual[k][j+(prank*elems)] + ephi2[edge][j+(i*prank)]*w2[i]*jDet*Flux.F[k];  
                    }
                }
            }   
            delete xy2d;
            delete det;
        }
    }
}

double* getMassMat(vector<vector<double>> &quad, int &p, vector<double> &Wq, int &Nq, double *Jmat, int &phi_size){
    double* xref = new double[2];
    double* M = new double[phi_size*phi_size]();
    double *phi = new double[phi_size];

    for (int i = 0; i < Nq; i++){
        xref[0] = quad[i][0];
        xref[1] = quad[i][1]; 

        shapeL(xref, p, &phi);

        for (int j = 0; j < phi_size; j++){
            for (int k = 0; k < phi_size; k++){
                //M[j] = phi[k]*phi[j]*Jmat[j]*Wq[j];
                M[k+(j*phi_size)] += phi[j]*phi[k]*Jmat[i]*Wq[i];
            }
        }
    }
    delete[] xref;
    delete[] phi;
    return M;
}

void interiorres(int p, int nelem, int nbasisfunctions, vector<vector<double>> &residual,vector<vector<double>> U, vector<vector<double>> &elem,  vector<vector<double>> &nodes,  vector<vector<double>> &curveNodes) {
//from integration order get quadrature points
//and also get basis functions
//evaluate basis functions at all quadrature points
//interpolate state to quad points using the matrix
//evaluate flux at each quad point
    int q;
    for (int i=0;i<nelem;++i){
        //interpolate state to quad points using the matrix
        
        int pts;
        if (curveNodes[i][3] == 0){
            q = 1;
        }
        else{
            q = 3;
        }

        quadOut out2D;

        out2D = quadature((2*p + 1)+2*(q - 1), 2);
        int N2 = out2D.Nq;

        double *quad2 = new double(N2*2);
        vector<double> w2(N2);
        double test[N2*2];
        for (int v = 0; v < 2*N2; v++){
            test[v] = out2D.x[v];
        }
        quad2 = test;
        for (int v=0;v<N2;++v){
            w2[v]=out2D.wq[v];
        } 
        //double *w2 = out2D.wq;
        
        double *pphi =  new double[nbasisfunctions];
        //double *xref = new double[2];
        double *pgphi =  new double[2*nbasisfunctions];

        outNorm normv;
        
        int size;
        if (curveNodes[i][3]!=0){
            size = 10;
        }
        else{
            size = 3;
        }
 
        vector<vector<double>> phi2(N2,vector<double> (nbasisfunctions));
        double gphi2[N2*nbasisfunctions];
        vector<vector<double>> coord(size,vector<double>(2));
        
        for (int j = 0; j < size; j++){   
            for (int k = 0; k < 2; k++){
                coord[j][k] = nodes[curveNodes[i][j]-1][k];
            }
        }
        
        for (int j = 0; j < N2; j++){
            double *xref = new double[2];
            xref[0] = quad2[j];
            xref[1] = quad2[j+N2];
            
            shapeL(xref, p, &pphi);
            
            for (int k = 0; k < nbasisfunctions; k++){
                phi2[j][k] = pphi[k];
            }
        }

        vector<vector<double>> invJ(4,vector<double>(N2));
        double* DetJ = jacobian(&quad2, coord, invJ,q,N2);

        vector<vector<double>> interp(N2,vector<double> (4));
        //vector<double> elemstate{1,.25*cos(8*3.14159/180),.25*sin(8*3.14159/180),1/(.4*1.4)+.25*.25/2};
        vector<vector<double>> elemstate(nbasisfunctions,vector<double> (4));
        for (int a=0;a<nbasisfunctions;++a){
            for(int b=0;b<4;++b){
                elemstate[a][b]=U[b][a+i*nbasisfunctions];
            }
        }
        for (int a=0;a<N2;++a){
            for (int b=0;b<4;++b){
                double sum=0;
                for (int c=0;c<nbasisfunctions;++c){
                    sum=sum+phi2[a][c]*elemstate[c][b];
                }
                interp[a][b]=sum;
            }
        }
        for (int j=0;j<nbasisfunctions;++j) {
            vector<double> rescont{0,0,0,0};
            double *xref = new double[2];
            
            for (int r=0;r<N2;r++) {
                xref[0] = quad2[r];
                xref[1] = quad2[r+N2];

                //grad is based on the basis function the residual is on, and the quad point
                double *gphi =  new double[2*nbasisfunctions];
                gradientL(xref, p, &gphi);
                double *grad = new double[2];
                grad[0]=gphi[j];
                grad[1]=gphi[j+nbasisfunctions];
                
                vector<double> hold{grad[0]*invJ[0][r]+grad[1]*invJ[2][r],grad[0]*invJ[1][r]+grad[1]*invJ[3][r]};
                
                //vector<double> hold{(grad[0]*J[1][1]-grad[1]*J[1][0])/(J[0][0]*J[1][1]-J[0][1]*J[1][0]),(-grad[0]*J[0][1]+grad[1]*J[0][0])/(J[0][0]*J[1][1]-J[0][1]*J[1][0])};
                double detJ=DetJ[r];
                double w=w2[r];
                hold[0]=hold[0]*detJ*w;
                hold[1]=hold[1]*detJ*w;

                vector<double> state=interp[r];
                double press=.4*(state[3]-.5*(state[1]*state[1]+state[2]*state[2])/state[0]);
                double H=(state[3]+press)/state[0];
                vector<double> Fx{state[1],state[1]*state[1]/state[0]+press,state[1]*state[2]/state[0],state[1]*H};
                vector<double> Fy{state[2],state[1]*state[2]/state[0],state[2]*state[2]/state[0]+press,state[2]*H};
                    
                rescont[0]=rescont[0]+hold[0]*Fx[0]+hold[1]*Fy[0];
                rescont[1]=rescont[1]+hold[0]*Fx[1]+hold[1]*Fy[1];
                rescont[2]=rescont[2]+hold[0]*Fx[2]+hold[1]*Fy[2];
                rescont[3]=rescont[3]+hold[0]*Fx[3]+hold[1]*Fy[3];
            }

            for (int k = 0; k < 4; k++){
                residual[k][j+nbasisfunctions*i] = rescont[k];
            }   

        }
    }  
}
vector<vector<double>> calcRes(int p, vector<vector<double>> interiorFaces, vector<vector<double>> bounds, vector<vector<double>> I2E, vector<vector<double>> B2E, vector<vector<double>> Bn, vector<vector<double>> In, vector<vector<double>> nodes, vector<vector<double>> U,vector<vector<double>> elem, vector<double> &speed, vector<vector<double>> curveNodes){
    int prank = (p+1)*(p+2)/2;
    int nelem = U[0].size()/prank;

    
    vector<vector<double>> residualE(4, vector<double>(prank*nelem));
    vector<vector<double>> residualInt(4, vector<double>(prank*nelem));
    vector<vector<double>> residual(4, vector<double>(prank*nelem));

    edgeRes(nelem, p, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, residualE, U, speed,curveNodes,elem);
    interiorres(p, nelem, prank, residualInt, U, elem, nodes,curveNodes);

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < prank*nelem; j++){
            residual[i][j] = residualE[i][j] - residualInt[i][j];
        }
    }

    return residual;
}
double l1norm(vector<vector<double>> residual, int prank, int nelem){
    double sumRes = 0;
    int size = residual[0].size();

    for (int i = 0; i < size; i++){
        for (int j = 0; j < 4; j++){
            sumRes += abs(residual[j][i]);
        }
    }
    return sumRes;
}
rk4stage rk4calc(int p, int nelem, int phi_size, int stage, vector<vector<double>> curveNodes,vector<vector<double>> nodes, vector<vector<double>> R, vector<vector<double>> U, vector<double> dt){
    quadOut outQuad;
    rk4stage output;

    vector<vector<double>> F(4,vector<double>(nelem*phi_size));
    vector<vector<double>> U_temp(4,vector<double>(nelem*phi_size));

    for (int i = 0; i < nelem; i++){
        int pts;
        int q;
        int prank = (p+1)*(p+2)/2;

        if (curveNodes[i][3] == 0){
            q = 1;
        }
        else{
            q = 3;
        }
        
        int o = 2*p + 1 + 2*(q - 1);
        outQuad = quadature(o, 2);
        int Nq = outQuad.Nq;
        double *quad2 = new double(Nq*2);
        vector<double> wq(Nq);
        double test[Nq*2];

        vector<vector<double>> quadV(Nq,vector<double>(2));
        vector<vector<double>> invJ(4,vector<double>(Nq));;

        for (int v = 0; v < 2*Nq; v++){
            test[v] = outQuad.x[v];
        }
        quad2 = test;
        
        for (int v=0;v<Nq;++v){
            wq[v]=outQuad.wq[v];
        }
        for (int j = 0; j < Nq; j++){
            quadV[j][0] = quad2[j];
            quadV[j][1] = quad2[j+Nq];

        }
        if (curveNodes[i][3]==0){
            pts = 3;
        }
        else{
            pts = 10;
        }
        vector<vector<double>> coord(pts,vector<double>(2));
        for (int j = 0; j < pts; j++){   
            for (int k = 0; k < 2; k++){
                coord[j][k] = nodes[curveNodes[i][j]-1][k];
            }
        }         

        double* Jmat = jacobian(&quad2, coord,invJ,q,Nq);
        double* massMat = getMassMat(quadV, p, wq, Nq, Jmat, phi_size);
        double* invMassMat = new double[prank*prank];
        invmat(massMat, invMassMat, phi_size);

        for (int j = 0; j < phi_size; j++){
            for (int k = 0; k < phi_size; k++){
                F[0][k+(i*phi_size)] = -invMassMat[k+(j*phi_size)]*R[0][k+(i*phi_size)];
                F[1][k+(i*phi_size)] = -invMassMat[k+(j*phi_size)]*R[1][k+(i*phi_size)];
                F[2][k+(i*phi_size)] = -invMassMat[k+(j*phi_size)]*R[2][k+(i*phi_size)];
                F[3][k+(i*phi_size)] = -invMassMat[k+(j*phi_size)]*R[3][k+(i*phi_size)];
            }
        }

        if (stage != 3){
            for (int j = 0; j < phi_size; j++){
                U_temp[0][j+(i*phi_size)] = U[0][j+(i*phi_size)] + 0.5*dt[i]*F[0][j+(i*phi_size)];
                U_temp[1][j+(i*phi_size)] = U[1][j+(i*phi_size)] + 0.5*dt[i]*F[1][j+(i*phi_size)];
                U_temp[2][j+(i*phi_size)] = U[2][j+(i*phi_size)] + 0.5*dt[i]*F[2][j+(i*phi_size)];
                U_temp[3][j+(i*phi_size)] = U[3][j+(i*phi_size)] + 0.5*dt[i]*F[3][j+(i*phi_size)];
            }
        }
        else if (stage == 3){
            for (int j = 0; j < phi_size; j++){
                U_temp[0][j+(i*phi_size)] = U[0][j+(i*phi_size)] + dt[i]*F[0][j+(i*phi_size)];
                U_temp[1][j+(i*phi_size)] = U[1][j+(i*phi_size)] + dt[i]*F[1][j+(i*phi_size)];
                U_temp[2][j+(i*phi_size)] = U[2][j+(i*phi_size)] + dt[i]*F[2][j+(i*phi_size)];
                U_temp[3][j+(i*phi_size)] = U[3][j+(i*phi_size)] + dt[i]*F[3][j+(i*phi_size)];
            }
        }
    }
    output.F = F;
    output.U = U_temp;

    return output;
}
rk4out dg(int p, int Nt, double CFL, vector<double> Area, vector<vector<double>> interiorFaces, vector<vector<double>> bounds, vector<vector<double>> I2E, vector<vector<double>> B2E, vector<vector<double>> Bn, vector<vector<double>> In, vector<vector<double>> nodes, vector<vector<double>> U, vector<vector<double>> elem, vector<double> &speed,vector<vector<double>> curveNodes) {
    rk4stage output0;
    rk4stage output1;
    rk4stage output2;
    rk4stage output3;

    rk4out results;

    int nelem = elem.size();

    int phi_size = (p+1)*(p+2)/2;

    vector<vector<double>> coord(phi_size,vector<double>(2));

    vector<vector<double>> F1(4,vector<double>(nelem*phi_size));
    vector<vector<double>> F2(4,vector<double>(nelem*phi_size));
    vector<vector<double>> F3(4,vector<double>(nelem*phi_size));
    vector<vector<double>> U_temp(4,vector<double>(nelem*phi_size));
    vector<vector<double>> Un(4,vector<double>(nelem*phi_size));
    vector<vector<double>> U1(4,vector<double>(nelem*phi_size));
    vector<vector<double>> R(4,vector<double>(nelem*phi_size));

    vector<double> dt(nelem);
    vector<double> resi;

    for (int t = 0; t < Nt; t++){
        // initial F calculation
        R = calcRes(p, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, U,elem, speed,curveNodes);
        dt = timestep(nelem, CFL, Area, speed);
        output0 = rk4calc(p, nelem, phi_size, 0, curveNodes, nodes,  R, U, dt);

        // 1st stage
        R = calcRes(p, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, output0.U, elem, speed,curveNodes);
        //dt = timestep(nelem, CFL, Area, speed);
        output1 = rk4calc(p, nelem, phi_size, 1, curveNodes, nodes,  R, U, dt);

        // 2nd stage
        R = calcRes(p, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, output1.U, elem, speed,curveNodes);
        //dt = timestep(nelem, CFL, Area, speed);
        output2 = rk4calc(p, nelem, phi_size, 2, curveNodes, nodes,  R, U, dt);

        // 3rd stage
        R = calcRes(p, interiorFaces, bounds, I2E, B2E, Bn, In, nodes, output2.U, elem, speed,curveNodes);
        //dt = timestep(nelem, CFL, Area, speed);
        output3 = rk4calc(p, nelem, phi_size, 3, curveNodes, nodes,  R, U, dt);

        // update state
        double Rl1 = l1norm(R, phi_size, nelem);
        resi.push_back(Rl1);
        
        for (int i = 0; i < nelem; i++){
            for (int j = 0; j < phi_size; j++){
                U[0][j+(i*phi_size)] = U[0][j+(i*phi_size)] + (1/6)*dt[j+(i*phi_size)]*(output0.F[0][j+(i*phi_size)] + 2*output1.F[0][j+(i*phi_size)] + 2*output2.F[0][j+(i*phi_size)] + output3.F[0][j+(i*phi_size)]);
                U[1][j+(i*phi_size)] = U[1][j+(i*phi_size)] + (1/6)*dt[j+(i*phi_size)]*(output0.F[1][j+(i*phi_size)] + 2*output1.F[1][j+(i*phi_size)] + 2*output2.F[1][j+(i*phi_size)] + output3.F[1][j+(i*phi_size)]);
                U[2][j+(i*phi_size)] = U[2][j+(i*phi_size)] + (1/6)*dt[j+(i*phi_size)]*(output0.F[2][j+(i*phi_size)] + 2*output1.F[2][j+(i*phi_size)] + 2*output2.F[2][j+(i*phi_size)] + output3.F[2][j+(i*phi_size)]);
                U[3][j+(i*phi_size)] = U[3][j+(i*phi_size)] + (1/6)*dt[j+(i*phi_size)]*(output0.F[3][j+(i*phi_size)] + 2*output1.F[3][j+(i*phi_size)] + 2*output2.F[3][j+(i*phi_size)] + output3.F[3][j+(i*phi_size)]);
            }
        }
        
        if ((t+1) % 250 == 0){
            cout << "Iteration: " << t+1 << endl;
        }
        if (Rl1 < 1e-5){
            //cout << "Iteration: " << t << endl;
            //break;
        }
       // cout<< t << endl;
    }
    for (int k = 0; k < 4; k++){
        for (int i = 0; i < nelem; i++){
            for (int j = 0; j < phi_size; j++){
                Un[k][j+(i*phi_size)]= U[k][j+(i*phi_size)];
            }
        }
    }
    results.residuals = resi;
    results.U = Un;

    return results;
}