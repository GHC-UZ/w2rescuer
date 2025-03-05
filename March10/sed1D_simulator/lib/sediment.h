/* Mud flow + Exner 1D Simulator
  ===========================================
   Weakly Coupled Model
  ===========================================================================================
   Sergio Martinez-Aranda    -    Dic 2020
  =========================================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../define.h"


int read_bedload_configuration(char *tempfile, int *bedloadCap,
	double *rhos, double *poros, double *ds, 
	double *shieldsC, double *nManSed,
	double *Ggrass, double *BetaB);

int read_bedload_boundary_configuration(char *tempfile,
	double *QSIN, double *ZBIN);

int h_initialize_bedload_variables( int nCells, int nWalls,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *DQs, double *lB);

int h_compute_initial_bedload_variables( int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr);

int h_set_bedload_inlet_initial_conditions(int nCells,
	double *h, double *M, double *Qs, double *qs,
	double QSIN);

int h_compute_bedload_mass(int nCells,
	double *Ab, double *dx, double *massBt);

int h_compute_bedload_time_step(int nWalls,
	double *A, double *u, double *zb, double *qs, 
	double *dx, 
	double ds, double poros,
	double *lB, double *dtB);	

int h_compute_bedload_wall_fluxes(int nWalls, 
    double *h, double *u, double dt,
    double *Ab, double *zb, double *Qs, double *qs, 
    double *lB, double *DQs,
	double ds, double poros,
    double *dx, int *solidWall);

int h_update_bedload_cells(int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr,
    double *DQs,
	double *dx, double dt);

int h_set_bedload_inlet_conditions( int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr,
    double dt, 
	double QSIN, double ZBIN);

int h_set_bedload_outlet_conditions( int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr);    
