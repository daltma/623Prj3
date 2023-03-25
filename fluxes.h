/* This function calculates the 2D Roe Flux
Inputs: UL, UR = left, right states, as 4x1 vectors
		gamma = ratio of specific heats (e.g. 1.4)
		n = left-to-right unit normal 2x1 vector
Outputs: F = numerical normal flux (4x1 vector)
		 smag = max wave speed estimate
		 */

#include <iostream>
#include <cmath>
#include <vector>
#include "tools.h"
#pragma once
using namespace std;

struct structFlux {
	vector<double> F;
	double smag;
};

structFlux roe(vector<double>& UL, vector<double>& UR, double gamma, vector<double> n) {
	const double gmi = gamma - 1.0;

	// process left state
	double rL, uL, vL, unL, qL, pL, rHL, HL, cL;

	rL = UL[0];
	uL = UL[1] / rL;
	vL = UL[2] / rL;
	unL = uL * n[0] + vL * n[1];
	qL = sqrt(pow(UL[1], 2) + pow(UL[2], 2)) / rL;
	pL = (gamma - 1) * (UL[3] - 0.5 * rL * pow(qL, 2));
	if (pL < 0 || rL < 0) {
		cout << "Non-physical state! A" << "\n";
		cout << UL[0] << "\n";
		cout << UL[1] << "\n";
	}
	if (pL < 0) {
		pL = -pL;
	}
	if (rL < 0) {
		rL = -rL;
	}
	rHL = UL[3] + pL;
	HL = rHL / rL;
	cL = sqrt(gamma * pL / rL);

	// left flux
	vector<double> FL(4, 0.0);
	FL[0] = rL * unL;
	FL[1] = UL[1] * unL + pL * n[0];
	FL[2] = UL[2] * unL + pL * n[1];
	FL[3] = rHL * unL;

	// process right state
	double rR, uR, vR, unR, qR, pR, rHR, HR, cR;

	rR = UR[0];
	uR = UR[1] / rR;
	vR = UR[2] / rR;
	unR = uR * n[0] + vR * n[1];
	qR = sqrt(pow(UR[1], 2) + pow(UR[2], 2)) / rR;
	pR = (gamma - 1) * (UR[3] - 0.5 * rR * pow(qR, 2));
	if (pR < 0 || rR < 0) {
		cout << "Non-physical state! B" << "\n";
	}
	if (pR < 0) {
		pR = -pR;
	}
	if (rR < 0) {
		rR = -rR;
	}
	rHR = UR[3] + pR;
	HR = rHR / rR;
	cR = sqrt(gamma * pR / rR);

	// right flux
	vector<double> FR(4, 0.0);
	FR[0] = rR * unR;
	FR[1] = UR[1] * unR + pR * n[0];
	FR[2] = UR[2] * unR + pR * n[1];
	FR[3] = rHR * unR;

	// difference in states
	vector<double> du;
	du = subtractVectors(UR, UL);

	// Roe average
	double di, d1, ui, vi, Hi, af, ucp, c2, ci, ci1;
	di = sqrt(rR / rL);
	d1 = 1.0 / (1.0 + di);

	ui = (di * uR + uL) * d1;
	vi = (di * vR + vL) * d1;
	Hi = (di * HR + HL) * d1;

	af = 0.5 * (ui * ui + vi * vi);
	ucp = ui * n[0] + vi * n[1];
	c2 = gmi * (Hi - af);
	if (c2 < 0) {
		cout << "Non-physical state! C" << "\n";
		c2 = -c2;
	}
	ci = sqrt(c2);
	ci1 = 1.0 / ci;

	// eigenvalues
	vector<double> l(3, 0.0);
	l[0] = ucp + ci; l[1] = ucp - ci; l[2] = ucp;

	// entropy fix
	double epsilon, l3;
	epsilon = ci * .1;
	for (int i = 0; i < 3; i++) {
		if ((l[i] < epsilon) && (l[i] > -epsilon)) {
			l[i] = 0.5 * (epsilon + l[i] * l[i] / epsilon);
		}
	}

	l = absolute(l); l3 = l[2];

	double s1, s2, G1, G2, C1, C2;
	// average and half-difference of 1st and 2nd eigs
	s1 = 0.5 * (l[0] + l[1]);
	s2 = 0.5 * (l[0] - l[1]);

	// left eigenvector product generators
	G1 = gmi * (af * du[0] - ui * du[1] - vi * du[2] + du[3]);
	G2 = -ucp * du[0] + du[1] * n[0] + du[2] * n[1];

	// required functions of G1 and G2(again, see Theory guide)
	C1 = G1 * (s1 - l3) * ci1 * ci1 + G2 * s2 * ci1;
	C2 = G1 * s2 * ci1 + G2 * (s1 - l3);

	// flux assembly
	vector<double> F(4, 0.0);
	F[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (l3 * du[0] + C1);
	F[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (l3 * du[1] + C1 * ui + C2 * n[0]);
	F[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (l3 * du[2] + C1 * vi + C2 * n[1]);
	F[3] = 0.5 * (FL[3] + FR[3]) - 0.5 * (l3 * du[3] + C1 * Hi + C2 * ucp);

	// max wave speed
	double smag;
	smag = max(l);

	structFlux output;
	output.F = F;
	output.smag = smag;
	return output;
}
structFlux wallFlux(const vector<double> &u,const vector<double> &n,const double &gam){
	vector<double> v(2);
    vector<double> vb(2);
   	vector<double> F(4);
    double pb;
    double smag;

	v = { u[1] / u[0], u[2] / u[0] };

	vb[0] = v[0] - (v[0] * n[0] + v[1] * n[1]) * n[0];
	vb[1] = v[1] - (v[0] * n[0] + v[1] * n[1]) * n[1];
	pb = (gam - 1) * (u[3] - .5 * u[0] * (pow(vb[0], 2) + pow(vb[1], 2)));
	if ((pb < 0) || (u[0] < 0)) { cout << "Non - physical state! D" << "/n"; }
	if (pb < 0) { pb = -pb; };
	F = { 0,pb * n[0],pb * n[1],0 };
	smag = abs(v[0] * n[0] + v[1] * n[1] + sqrt(gam * pb / u[0]));
    
    structFlux output;
	output.F = F;
	output.smag = smag;
	return output;
}