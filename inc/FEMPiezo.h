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

#ifndef FEMPIEZO_H // FEMPIEZO_H
#define FEMPIEZO_H

// Includes
#include "defs.h"
#include "FEMBody.h"


// Forward declarations
class FEMBodyClass;

// FEM piezo class
class FEMPiezoClass {

    // Friend classes


// Default constructor and destructor
public:
	FEMPiezoClass() {fPtr = NULL; piezoDOFs = 0; bcDOFs = 0; itNR = 0; resNR = 0.0; subRes = 0.0; subNum = 0.0; subDen = 0.0;};
	~FEMPiezoClass() {};

	// Custom constructor for building corotational FEM body
	FEMPiezoClass(FEMBodyClass *fBodyPtr);


    // Private members
private:

    // Pointer to FEM body
	FEMBodyClass *fPtr;

    // Geometry

    // System values
	int piezoDOFs;				// DOFs for whole body
	int bcDOFs;					// Number of DOFs removed when applying BCs
	int itNR;					// Number of iterations for Newton-Raphson solver
    double resNR;				// Residual Newton-Raphson solver reached


    // System matrices
    vector<double> Mp;						// Mass matrix
    vector<double> Dp;						// Damping matrix
	vector<double> Kp;						// Stiffness matrix
    vector<double> Rp;						// Load vector
    vector<double> Fp;						// Vector of internal forces
    vector<double> X;						// Vector of unknown
    vector<double> delX;					// Vector of incremental displacements
	vector<double> Xdot;					// Vector velocities
	vector<double> Xdotdot;					// Vector of accelerations
	vector<double> X_n;						// Vector of displacements at start of timestep
	vector<double> X_nm1;					// Vector of displacements at last time step
	vector<double> X_nm2;					// Vector of displacements at two time steps ago
	vector<double> Xdot_n;					// Vector velocities at start of timestep
	vector<double> Xdotdot_n;				// Vector of accelerations at start of timestep
	vector<double> X_km1;					// Vector of displacements at last iteration

    // Parameters for Aitken relaxation
	double subRes;				// Residual
	double subNum;				// Numerator
	double subDen;				// Denominator
    vector<double> Rp_k;			// Vector of nodal residuals
	vector<double> Rp_km1;		// Vector of nodal residuals at last iteration

    // Private methods
private:

    // FEM methods


};

#endif // FEMPIEZO_H