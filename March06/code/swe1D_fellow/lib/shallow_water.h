//   Shallow Water Equations 1D Simulator (SWE1D)
//  ===========================================================================================
//   Rectangular cross section
//  ===========================================================================================
//   Version 1.0 - Jan 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../define.h"


int read_geometry_file(char *tempfile, int *nCells, int *nWalls, 
	double *DeltaX, double *Xmin, double *Width, double *Length );

int read_bed_configuration(char *tempfile, int *iniBed,
	double *ZBupstream, double *ZBslope, 
	double *XBstep, double *ZBL, double *ZBR,
	double *XBcoast, double *ZBbase, double *slopeCoast);

int read_flow_configuration(char *tempfile, int *iniCond,
	double *Uconst, double *Hconst, 
	double *Xdam, double *ZSL, double *ZSR,
	double *Xstep, double *HL, double *UL, double *HR, double *UR);	

int read_simulation_setup(char *tempfile, double *simTime, 
	double *Toutput, double *CFL, double *n1);	

int read_boundary_configuration(char *tempfile, 
	double *QIN, double *HIN, 
	double *HOUT, double *ZSOUT);	

int h_initialize_variables(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr, 
	double *B, double *P, 
	double *n, double *zb, 	
	double *DU1, double *DU2,
	double *lSW);

int h_set_bed_initial_conditions(int nCells,
	double *x, double *dist, double *zb, double *zr,
	int iniBed,
	double ZBupstream, double ZBslope, 
	double XBstep, double ZBL, double ZBR,
	double XBcoast, double ZBbase, double slopeCoast);	

int h_set_flow_initial_conditions(int nCells,
	double *x, double *zb, double *h, double *u, 
	int iniCond,
	double Uconst, double Hconst, 
	double Xdam, double ZSL, double ZSR,
	double Xstep, double HL, double UL, double HR, double UR);

int h_set_inlet_initial_conditions(int nCells,
	double *h, double *u, double *M,
	double QIN, double HIN);

int h_set_outlet_initial_conditions(int nCells,
	double *h, double *u, double *M, double *zb,
	double HOUT, double ZSOUT);

int h_compute_initial_flow_variables(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr, 
	double *B, double *P, 
	double *M, double *zb,
	double *n, double n1);

int h_compute_water_mass(int nCells,
	double *A, double *dx, double *massWt);

int h_compute_flow_time_step(int nWalls, 
	double *A, double *u, double *B, double *dx, 
	double *lSW, double *dtSW);

int h_compute_wall_fluxes(int nWalls,
	double *A, double *Q, double *h, double *q, double *u, 
  	double *c, double *Fr,
	double *B, double *P,
  	double *n, double *zb,
	double *DU1, double *DU2, 
	double *dx);

int h_check_depth_positivity(int nCells, 
	double *A, double *DU1, double *dx, double *dt);

int h_update_cells(int nCells, 
	double *A, double *Q, double *h, double *q, double *u, 
  	double *c, double *Fr,
	double *B, double *P,
	double *DU1, double *DU2, 
	double *M, double *dx,
  	double dt);

int h_set_inlet_conditions(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr,
	double *B, double *P, 
	double *M,
	double QIN, double HIN);


int h_set_outlet_conditions(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr,
	double *B, double *P, 
	double *M, double *zb,
	double HOUT, double ZSOUT);