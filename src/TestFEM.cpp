/*
    LIFE: Lattice boltzmann-Immersed boundary-Finite Element
    Copyright (C) 2019 Joseph O'Connor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Includes
#include "../inc/FEMBody.h"


class TestFEMClass : public FEMBodyListenerClass
{
	public:

	TestFEMClass() {
		nNodes = static_cast<int>(floor(length / Dx())) + 1;
		cout << "nNodes = " << nNodes << endl;

		for (int i = 0; i < nNodes; i++) {
			pos.push_back({i * Dx(), 0});
			vel.push_back({0,0});
		}
		time = 0;

		fBody = new FEMBodyClass(this, {0,0}, {length, width}, 0, "10", "CLAMPED", rho, ym);
		fBody->dynamicFEM();
		ofstream of("pos.dat");
		printNodePositions(of);
	}

	void printNodePositions(ostream &os) {
		os << "x" << "\t" << "y" << endl;
		for (const auto &posi : pos) {
			os << posi[0] << "\t" << posi[1] << endl;
		}
	}

	double Dt() {return 0.00001;};
	double t() {return time;};
	double Dx() {return 1./10;};
	double Dm() {return 1.0;};
	int bodyID() {return 0;};
	void setNodePosition(int nodeID, array<double, dims> newpos) {
		cout << "Setting node " << nodeID << " position to " << newpos[0] << "," << newpos[1] << endl; 
		pos[nodeID] = newpos;
	};
	void setNodeVelocity(int nodeID, array<double, dims> newvel) {vel[nodeID] = newvel;};
	int numNodes() {return nNodes;};
	double epsilon(int /*nodeID*/) {return 2;};
	array<double, dims> force(int nodeID) {
		if (nodeID == nNodes-1) {
			const double Mcr = 2*M_PI*ym*secMomArea/length;
			const double tipMoment = 0.05 * Mcr;
			const double tipForce = tipMoment / length;
			const double tipForceDensity = tipForce * Dx();

			cout << "tipForceDensity = " << tipForceDensity << endl;

			return array<double, dims>({0,tipForce});
		}
		else {
			return array<double, dims>({0,0});
		}
	}

	void advanceTime() {time += Dt();};
		
	const double ym = 1.4e6;
	const double rho = 10000.0;
	const double length = 1.0;
	const double width = 0.02;
	const double secMomArea = pow(width,4)/12;

	array<double, dims> getNodePosition(int nodeID) {return pos[nodeID];};

	double getTime() {return time;};

	private:
	FEMBodyClass  *fBody;

	vector<array<double, dims>> pos;
	vector<array<double, dims>> vel;
	int nNodes;
	double time;
};

// ***** Main function ***** //
int main() {

	cout << "*** INITIALISING TestFEMMain ***" << endl;

	TestFEMClass testFEM;

	// FEMBodyClass  fBody(&testFEM, {0,0}, {0.35,0.02}, 0, "20", "CLAMPED", 10000.0, 1.4e6);
	
	// for (int i = 0; i < 20; i++) {
	// 	auto pos = testFEM.getNodePosition(87);
	// 	cout << "Tip position at time " << testFEM.getTime() << " : " << pos[0] << " " << pos[1] << endl;
	// 	fBody.dynamicFEM();
	// 	testFEM.advanceTime();
	// }

// Try the built-in gravity weight thing. Also change force to force density via 		F = Tsub * (((-epsilon * 1.0 * forceScale) * force) + weight); in FEMElement.cpp

	cout << "TestFEMMain done." << endl;
}
