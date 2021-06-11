
#ifndef PARAMS_H	// PARAMS_H
#define PARAMS_H

// Include defs
#include "defs.h"

// Set number of OMP threads (if commented then it will use system max)
//#define THREADS 1

// Set resFactor (for easily changing mesh resolution)
const int resFactor = 3;

// Simulation options
//#define CENTRAL_MOMENTS			// Use central moments collision operator (BGK otherwise)
//#define INLET_RAMP 2.0			// Inlet velocity ramp-up
//#define WOMERSLEY 5.0				// Womersley number for oscillating pressure gradients
//#define UNI_EPSILON				// Calculate epsilon over all IBM bodies at once
#define ORDERED						// For deterministic reduction operations
//#define INITIAL_DEFLECT 0.01		// Set an initial deflection (fraction of L)
#define PIEZO_EFFECT

// Outputs
#define VTK								// Write out VTK
//#define VTK_FEM						// Write out the FEM VTK
#define FORCES							// Write out forces on structures
#define TIPS							// Write out tip positions
#define ENERGIES						// Write out Strain Energies

// Domain setup (lattice)
const int Nx = resFactor * 10 * 10 + 1;   	// Number of lattice sites in x-direction
const int Ny = resFactor * 10 * 5 + 1;		// Number of lattice sites in y-direction

// Domain setup (physical)
const double height_p = 0.2;			// Domain height (m)
const double rho_p = 1.2047;			// Fluid density (kg/m^3)
const double nu_p = 0.5e-2;			// Fluid kinematic viscosity (m^2/s)

// Initial conditions
const double ux0_p = 10.0;				// Initial x-velocity (m/s)
const double uy0_p = 0.0;				// Initial y-velocity (m/s)

// Gravity and pressure gradient
const double gravityX = 0.0;			// Gravity component in x-direction (m/s^2)
const double gravityY = 10.0;			// Gravity component in y-direction (m/s^2)
const double dpdx = 0.0;				// Pressure gradient in x-direction (Pa/m)
const double dpdy = 0.0;				// Pressure gradient in x-direction (Pa/m)

// Boundary conditions (set to eFluid for periodic)
#define WALL_LEFT	eVelocity			// Boundary condition at left wall
#define WALL_RIGHT	eConvective			// Boundary condition at right wall
#define WALL_BOTTOM	eFluid			// Boundary condition at bottom wall
#define WALL_TOP	eFluid			// Boundary condition at top wall
//#define PROFILE 	eParabolic			// Profile shape (uncomment for uniform)

// Inlet conditions
const double uxInlet_p = ux0_p;			// Inlet x-velocity (m/s)
const double uyInlet_p = 0.0;			// Inlet y-velocity (m/s)

// FEM Newmark integration parameters
const double alpha = 0.25;				// Alpha parameter
const double delta = 0.5;				// Delta parameter

// FSI Coupling parameters
const double relaxMax = 1.0;			// Relaxation parameter (initial guess)
const double subTol = 1e-5;				// Tolerance for subiterations

// Set omega by three different methods (comment/uncomment)
// Set omega directly
//const double omega = 1.875;
//const double tStep = pow(1.0 / sqrt(3.0), 2.0) * pow(height_p / (Ny - 1), 2.0) * (1.0 - 0.5 * omega) / (nu_p * omega);

// Set time step
const double tStep = 0.00005 / (resFactor * resFactor);
const double omega = 1.0 / (nu_p * tStep / (pow(1.0 / sqrt(3.0), 2.0) * pow(height_p / (Ny - 1), 2.0)) + 0.5);

// Set lattice velocity
//const double uLB = (2.0 / 3.0) * 0.2 * 1.0 / sqrt(3.0);
//const double omega = 1.0 / (uLB * (Ny-1) / (uxInlet_p * height_p / (3.0 * nu_p)) + 0.5);

// Number of time steps and how often to write out
const int nSteps = 100;	//static_cast<int>(round(1 / tStep));	// Number of timesteps
const int tinfo = 1;	//nSteps / 1000;						// Frequency to write out info and logs
const int tVTK = 1;		//nSteps / 1000;						// Frequency to write out VTK
const int tRestart = nSteps / 1;								// Frequency to write out restart files

// Reference values
const double ref_nu = nu_p;			// Reference kinematic viscosity (m^2/s)
const double ref_rho = rho_p;		// Reference density (kg/m^3)
const double ref_P = 0.0;			// Reference pressure (Pa)
const double ref_L = 0.05;			// Reference length (m)
const double ref_U = uxInlet_p;		// Reference velocity (m/s)

// Precision for ASCII output
const int PRECISION = 10;

#endif	// PARAMS_H
