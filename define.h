//   Test code
//  ===========================================================================================
//   Version 1.0 - Jan 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>       			// clock_t, clock(), CLOCKS_PER_SEC
#include "lib/gnuplot-iostream.h"		// Graphical output (linux only)
#include <unistd.h>					// Para usleep()



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Mesh discretization
#define ZRlevel 2.0			  // Reference bed elevation

#define simTime 600      // Simulation time (s)
#define Toutput 60       // Output time (s)
#define nCells 500		// Number of cells
#define dx 2.0		    // Spatial increment [m]
#define a 0.5			// Advection velocity [m/s]
#define b 0.02			// Diffusion coefficient
#define C1 0.01
#define C2 200



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simulation setup
#define Niter 10       // Number of iterations for screen output
#define DISPLAY 1

#define rhow 1000.0   // Water density
#define mu 0.001			// Water dynamic viscosity
#define hmin 0.001    // Minimum depth for advective flux



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Phisical-mathematical parameters
#define PI 3.14159265358979323846
#define g 9.81
#define tol6 1e-6
#define tol9 1e-9
#define tol12 1e-12



// Macros
#define MIN(x,y) (((x)<(y)) ? (x) : (y))			// Macro to find the minimum of two numbers
#define MAX(a,b) (((a)>(b)) ? (a) : (b))            // Macro to find the maximum of two numbers


