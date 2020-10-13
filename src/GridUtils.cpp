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
#include "../inc/GridUtils.h"
#include "../inc/Utils.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"

// Read in restart file
void GridUtils::readRestart(GridClass &grid) {

	// Starting main algorithm
	cout << endl << endl << endl << "*** READ IN RESTART DATA ***";

	// Reading in restart files
	cout << endl << endl << "Reading in restart files...";

	// Read grid restart data
	grid.readRestart();

	// Read bodies restart data
	grid.oPtr->readRestart();

	// Delete VTK files that were written after last restart written
	GridUtils::deleteVTKs(grid);

	// Write out header
	cout << "finished";
}

// Write out restart file
void GridUtils::writeRestart(GridClass &grid) {

	// Wrting restart files
	cout << endl << "Writing restart data...";

	// Write grid info
	grid.writeRestart();

	// Write IBM info
	grid.oPtr->writeRestart();

	// Write header
	cout << "finished";
}

// Write info
void GridUtils::writeInfo(GridClass &grid) {

	// Write grid info
	grid.writeInfo();

	// Write IBM info
	grid.oPtr->writeInfo();

	// Write out forces on bodies
#ifdef FORCES
	grid.oPtr->writeTotalForces();
#endif

	// Write out forces on bodies
#ifdef TIPS
	grid.oPtr->writeTips();
#endif
}

// Write log
void GridUtils::writeLog(GridClass &grid) {

	// Write grid info
	grid.writeLog();

	// Write IBM info
	grid.oPtr->writeLog();
}

// Write VTK
void GridUtils::writeVTK(GridClass &grid) {

	// Write VTK
#ifdef VTK

	// Write header
	cout << endl << "Writing VTK data...";

	// Write grid VTK
	grid.writeVTK();

	// Write IBM VTK
	grid.oPtr->writeVTK();

	// Write header
	cout << "finished";
#endif
}

// Delete future VTKs
void GridUtils::deleteVTKs(GridClass &grid) {

	// Results path
	string path = "Results/VTK";

	// Check if it exists
	if (boost::filesystem::exists(path)) {

		// Loop through three different file types
		for (int i = 0; i < 3; i++) {

			// File strings
			string fileStr, extStr;

			// Get types
			if (i == 0) {
				fileStr = "Fluid";
				extStr = ".vti";
			}
			else if (i == 1) {
				fileStr = "IBM";
				extStr = ".vtp";
			}
			else if (i == 2) {
				fileStr = "FEM";
				extStr = ".vtp";
			}

			// Get directory iterator
			boost::filesystem::directory_iterator endit;

			// Loop through all fluid vtk files
			for (boost::filesystem::directory_iterator it(path); it != endit; it++) {

				// Get string
				string fStr = it->path().stem().string();

				// Check to make sure it is the files we want
				if (fStr.find(fileStr) != string::npos && it->path().extension() == extStr) {

					// Get time step
					string tStepStr = fStr.substr(fileStr.size() + 1, fStr.size());

					// Convert to int
					int tStep = stoi(tStepStr);

					// If tStep is bigger than current time then delete
					if (tStep >= grid.t) {
						if (!boost::filesystem::remove(path + "/" + fStr + extStr))
							ERROR("Problem removing future VTK files...exiting");
					}
				}
			}
		}
	}
}
