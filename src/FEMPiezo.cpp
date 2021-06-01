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








// Custom constructor for building FEM piezo
FEMPiezoClass::FEMPiezoClass(FEMBodyClass *fBodyPtr) {

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
	piezoDOFs = fPtr->bodyDOFs; // + qqc lié avec Q



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