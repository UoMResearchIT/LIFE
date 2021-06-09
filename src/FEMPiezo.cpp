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
#include "../inc/FEMPiezo.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"


// Newton raphson iterator
void FEMPiezoClass::newtonRaphsonDynamic() {

	// Build global matrices
	buildGlobalMatrices();

	// Apply Newmark scheme (using Newmark coefficients)
	setNewmark();

	// Solve linear system using LAPACK library
	delX = Utils::solveLAPACK(Kp, Fp, bcDOFs);

	// Add deltaU to X
	X = X + delX;

	// Update U
	fPtr->U = X;
	fPtr->U.pop_back();

	// Update FEM positions
	fPtr->updateFEMValues();
}


// Check convergence of Newton Raphson iterator
inline double FEMPiezoClass::checkNRConvergence () {

	// Get the norm of delU
	return sqrt(delX * delX) / (ref_L * sqrt(static_cast<double>(delX.size())));
}


// Build global matrices
void FEMPiezoClass::buildGlobalMatrices() {

	// Set matrices to zero
	fill(Mp.begin(), Mp.end(), 0.0);
	fill(Dp.begin(), Dp.end(), 0.0);
	fill(Kp.begin(), Kp.end(), 0.0);
	fill(Fp.begin(), Fp.end(), 0.0);

	// Build global matrices Mp
	vector<double> Mm = fPtr->M;
		// Copy of Mm
	int dim = fPtr->bodyDOFs;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			Mp[i * piezoDOFs + j] = Mm[i * dim + j];
		}
	}
		// Add piezo characteristic
	Mp[piezoDOFs * piezoDOFs -1] = L;

	// Build global matrices Dp
	Dp[piezoDOFs * piezoDOFs -1] = Rohm;

	// Build global matrices Kp
	vector<double> Km = fPtr->K;
		// Build C
	double C = dielec_cst * fPtr->L0 / hp; // Non dimensionalised by the Width
		// Build K1 and K2 (with the hyp of one patch over all the length on each side and for a single flag)
	vector<double> K1;
	vector<double> K2;
	K1.resize(dim, 0.0);
	K2.resize(dim, 0.0);

	K1[0] = -1*piezo_cst; // Non dimensionalised by the Width
	K1[2] = piezo_cst * (h+hp)/2; // Non dimensionalised by the Width
	K1[dim-3] = piezo_cst; // Non dimensionalised by the Width
	K1[dim-1] = -1*piezo_cst * (h+hp)/2; // Non dimensionalised by the Width

	K2[0] = piezo_cst; // Non dimensionalised by the Width
	K2[2] = piezo_cst * (h+hp)/2; // Non dimensionalised by the Width
	K2[dim-3] = -1*piezo_cst; // Non dimensionalised by the Width
	K2[dim-1] = -1*piezo_cst * (h+hp)/2; // Non dimensionalised by the Width
		// Build Kp
	vector<double> Kp;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			Kp[i * piezoDOFs + j] = Km[i * dim + j] + (K1[i]*K1[j] + K2[i]*K2[j])/C;
		}
	}
	for (size_t ij = 0; ij < dim; ij++) {
		Kp[(piezoDOFs-1) * piezoDOFs + ij] = (K1[ij] + K2[ij])/C;
		Kp[ij * piezoDOFs + (piezoDOFs-1)] = (K1[ij] + K2[ij])/C;
	}
	Kp[piezoDOFs*piezoDOFs-1] = 2/C;

	// Build global matrices Fp
	vector<double> F = fPtr->F;
		// Copy of F
	for (size_t i = 0; i < F.size(); i++) {
		Fp[i] = F[i];
	}

	// Build global matrices Rp
	vector<double> R = fPtr->R;
		// Copy of R
	for (size_t i = 0; i < R.size(); i++) {
		Rp[i] = R[i];
	}
}


// Set Newmark
void FEMPiezoClass::setNewmark() {

	// Newmark-beta method for time integration
	double Dt = fPtr->iPtr->oPtr->gPtr->Dt;
	double a0, a2, a3, a11, a9, a12, a13;
	a0 = 1.0 / (alpha * SQ(Dt));
	a2 = 1.0 / (alpha * Dt);
	a3 = 1.0 / (2.0 * alpha) - 1.0;
	a9 = Dt * (1.0 - delta);
	a11 = delta / (alpha * Dt);
	a12 = delta / alpha;
	a13 = Dt * delta * a3;

	// Calculate effective load vector
	Fp = Rp - Fp + Utils::MatMultiply(Mp, a0 * X + a2 * Xdot + a3 * Xdotdot)
				 + Utils::MatMultiply(Dp, (1 - a12) * Xdot + (a9 - a13) * Xdotdot - a11 * X);

	// Calculate effective stiffness matrix
	Kp = Kp + a0 * Mp + a11 * Dp;
}


// Finish Newmark
void FEMPiezoClass::finishNewmark() {

	// Get timestep
	double Dt = fPtr->iPtr->oPtr->gPtr->Dt;

	// Newmark coefficients
	double a6 = 1.0 / (alpha * SQ(Dt));
	double a7 = -1.0 / (alpha * Dt);
	double a8 = -(1.0 / (2.0 * alpha) - 1.0);
	double a9 = Dt * (1.0 - delta);
	double a10 = delta * Dt;

	// Update velocities and accelerations
	Xdotdot = a6 * (X - X_n) + a7 * Xdot_n + a8 * Xdotdot_n;
	Xdot = Xdot_n + a9 * Xdotdot_n + a10 * Xdotdot;
}


// Reset the start of time step values
void FEMPiezoClass::resetValues() {

	// Reset start of time step values
	X_nm2.swap(X_nm1);
	X_nm1.swap(X_n);
	X_n.swap(X);
	Xdot_n.swap(Xdot);
	Xdotdot_n.swap(Xdotdot);
}


// Custom constructor for building FEM piezo
FEMPiezoClass::FEMPiezoClass(FEMBodyClass *fBodyPtr, double h, double hp, double piezo_cst, double dielec_cst, double Rohm, double L) {

	// Set pointer
	fPtr = fBodyPtr;

	// Set sub residaul, numerator and denominator to initial value
	subRes = 0.0;
	subNum = 0.0;
	subDen = 0.0;

	// Get number of DOFs required for BC
	bcDOFs = fPtr->bcDOFs;

	// Unpack geometry
	

	

	// Get number of DOFs in piezobody
	piezoDOFs = fPtr->bodyDOFs + 1; // Add in function of the number of flags

	// Size the matrices
	Mp.resize(piezoDOFs * piezoDOFs, 0.0);
	Dp.resize(piezoDOFs * piezoDOFs, 0.0);
	Kp.resize(piezoDOFs * piezoDOFs, 0.0);
	Rp.resize(piezoDOFs, 0.0);
	Fp.resize(piezoDOFs, 0.0);
	X.resize(piezoDOFs, 0.0);
	delX.resize(piezoDOFs, 0.0);
	Xdot.resize(piezoDOFs, 0.0);
	Xdotdot.resize(piezoDOFs, 0.0);
	X_n.resize(piezoDOFs, 0.0);
	X_nm1.resize(piezoDOFs, 0.0);
	X_nm2.resize(piezoDOFs, 0.0);
	Xdot_n.resize(piezoDOFs, 0.0);
	Xdotdot_n.resize(piezoDOFs, 0.0);
	X_km1.resize(piezoDOFs, 0.0);
	Rp_k.resize(piezoDOFs, 0.0);
	Rp_km1.resize(piezoDOFs, 0.0);
}