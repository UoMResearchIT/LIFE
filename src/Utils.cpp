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
#include "../inc/Utils.h"

using namespace LIFE;

// Write header at start of run
void Utils::writeHeader() {

	// Write out version number and release date
	cout << endl << "*** LIFE " << version << " (" << date << ") ***" << endl << endl;

	// Write out license header
	cout << "LIFE: Lattice boltzmann-Immersed boundary-Finite Element" << endl;
	cout << "Copyright (C) 2019 Joseph O'Connor" << endl << endl;

	cout << "This program is free software: you can redistribute it and/or modify" << endl;
	cout << "it under the terms of the GNU General Public License as published by" << endl;
	cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
	cout << "(at your option) any later version." << endl << endl;

	cout << "This program is distributed in the hope that it will be useful," << endl;
	cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
	cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
	cout << "GNU General Public License for more details." << endl << endl;

	cout << "You should have received a copy of the GNU General Public License" << endl;
	cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
}

// Create directories at start
bool Utils::createDirectories() {

	// First check if restart file exists
	if (boost::filesystem::exists("Results/Restart/Fluid.restart"))
		return true;

	// Remove directories
	if (boost::filesystem::exists("Results")) {
		if (!boost::filesystem::remove_all("Results"))
			ERROR("Problem removing results directory...exiting");
	}

	// Create results directory
	if (!boost::filesystem::create_directory("Results"))
		ERROR("Problem creating results directory...exiting");

	// Create VTK directory
#ifdef VTK
	if (!boost::filesystem::create_directory("Results/VTK"))
		ERROR("Problem creating vtk directory...exiting");
#endif

	// Create ASCII directory
#ifdef ASCII
	if (!boost::filesystem::create_directory("Results/dat"))
		ERROR("Problem creating dat directory...exiting");
#endif

	// Create YAML directory
#ifdef YAML
	if (!boost::filesystem::create_directory("Results/YAML"))
		ERROR("Problem creating YAML directory...exiting");
#endif

	// Only create restart directory if we need to
	if (tRestart > 0) {
		if (!boost::filesystem::create_directory("Results/Restart"))
			ERROR("Problem creating restart directory...exiting");
	}

	// If we make it here then this is not a restart case
	return false;
}


// Get number of omp threads (as built in doesn't work on GCC)
int Utils::omp_thread_count() {

	// Thread counters
	int n = 0;

	// Do parallel reduction on n
#pragma omp parallel reduction(+:n)
	n += 1;

	// Return total thread count
	return n;
}

// Convert seconds to hours:minutes:seconds
array<int, 3> Utils::secs2hms(double seconds) {

	// Round to nearest second
	seconds = round(seconds);

	// Get number of hours and number of minutes
	double hours = seconds / (60.0 * 60.0);
	double minutes = seconds / (60.0);

	// Declare vector
	array<int, 3> hms;

	// Get hours:minutes:seconds
	hms[0] = static_cast<int>(floor(hours));
	hms[1] = static_cast<int>(floor(minutes - (static_cast<double>(hms[0]) * 60.0)));
	hms[2] = static_cast<int>(floor(static_cast<double>(seconds) - static_cast<double>(hms[1]) * 60.0 - static_cast<double>(hms[0]) * 60.0 * 60.0));

	// Return
	return hms;
}

// Get string for boundary condition
string Utils::getBoundaryString(eLatType BCType) {

	// String
	string str;

	// Check against possible options
	if (BCType == eFluid)
		str = "Periodic BC";
	else if (BCType == eVelocity)
		str = "Velocity BC";
	else if (BCType == ePressure)
		str = "Pressure BC";
	else if (BCType == eWall)
		str = "Wall BC";
	else if (BCType == eFreeSlip)
		str = "Free Slip BC";
	else if (BCType == eConvective)
		str = "Convective BC";

	// Return
	return str;
}

// Solve linear system using LAPACK routines
vector<double> Utils::solveLAPACK(vector<double> A, vector<double> b, int BC) {

	// Set up the correct values
	char trans = 'T';
	int dim = static_cast<int>(sqrt(static_cast<int>(A.size())));
	int row = dim - BC;
	int col = dim - BC;
	int offset = BC * dim + BC;
    int nrhs = 1;
    int LDA = dim;
    int LDB = dim;
    int info;
	std::vector<int> ipiv(row, 0);

    // Factorise and solve
	dgetrf_(&row, &col, A.data() + offset, &LDA, ipiv.data(), &info);
	dgetrs_(&trans, &row, &nrhs, A.data() + offset, &LDA, ipiv.data(), b.data() + BC, &LDB, &info);

	// Set return values not included to zero
	fill(b.begin(), b.begin() + BC, 0.0);

	// Return RHS
	return b;
}
