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

#ifndef GRIDUTILS_H // UTILS_H
#define GRIDUTILS_H

// Includes
#include <array>
#include "params.h"

// Forward declarations
class GridClass;

// Utility namespace for helper functions
namespace GridUtils {

// DECLARATIONS

// Read in restart file
void readRestart(GridClass &grid);

// Read in restart file
void writeRestart(GridClass &grid);

// Write out info
void writeInfo(GridClass &grid);

// Write out log
void writeLog(GridClass &grid);

// Write out VTK
void writeVTK(GridClass &grid);

// Delete future VTKs
void deleteVTKs(GridClass &grid);

} // namespace GridUtils

#endif // UTILS_H
