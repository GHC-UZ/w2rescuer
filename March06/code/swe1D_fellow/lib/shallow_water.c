//   Shallow Water Equations 1D Simulator (SWE1D)
//  ===========================================================================================
//   Rectangular cross section
//  ===========================================================================================
//   Version 1.0 - Jan 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 


#include "shallow_water.h"


int read_geometry_file(char *tempfile, int *nCells, int *nWalls, 
	double *DeltaX, double *Xmin, double *Width, double *Length){

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;
	long iaux;
	double daux;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Channel geometry definition
    while (fgets(line, sizeof(line), fp)) {
        
		item = "nCells:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*nCells) = (int)strtol(position, &endPtr, 10); //integer
        }

		item = "Xmin:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Xmin) = (double)strtod(position, &endPtr); //double 
        }

		item = "DeltaX:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*DeltaX) = (double)strtod(position, &endPtr); //double 
        }		

		item = "Width:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Width) = (double)strtod(position, &endPtr); //double 
        }

    }
	fclose(fp);	

	(*nWalls) = (*nCells)-1; //Number of walls 
	(*Length) = (*Xmin) + (*DeltaX)*(*nCells); //Channel length

	return 1;
}



int read_bed_configuration(char *tempfile, int *iniBed,
	double *ZBupstream, double *ZBslope, 
	double *XBstep, double *ZBL, double *ZBR,
	double *XBcoast, double *ZBbase, double *slopeCoast){		

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;
	long iaux;
	double daux;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {
        
		item = "iniBed:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*iniBed) = (int)strtol(position, &endPtr, 10); //integer
        }

		item = "ZBupstream:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZBupstream) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZBslope:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZBslope) = (double)strtod(position, &endPtr); //double 
        }		

		item = "XBstep:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*XBstep) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZBL:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZBL) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZBR:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZBR) = (double)strtod(position, &endPtr); //double 
        }

		item = "XBcoast:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*XBcoast) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZBbase:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZBbase) = (double)strtod(position, &endPtr); //double 
        }

		item = "slopeCoast:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*slopeCoast) = (double)strtod(position, &endPtr); //double 
        }										

    }
	fclose(fp);

	return 1;
}



int read_flow_configuration(char *tempfile, int *iniCond,
	double *Uconst, double *Hconst, 
	double *Xdam, double *ZSL, double *ZSR,
	double *Xstep, double *HL, double *UL, double *HR, double *UR){		

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;
	long iaux;
	double daux;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {
        
		item = "iniCond:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*iniCond) = (int)strtol(position, &endPtr, 10); //integer
        }

		item = "Uconst:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Uconst) = (double)strtod(position, &endPtr); //double 
        }

		item = "Hconst:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Hconst) = (double)strtod(position, &endPtr); //double 
        }		

		item = "Xdam:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Xdam) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZSL:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZSL) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZSR:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZSR) = (double)strtod(position, &endPtr); //double 
        }

		item = "Xstep:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Xstep) = (double)strtod(position, &endPtr); //double 
        }

		item = "HL:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*HL) = (double)strtod(position, &endPtr); //double 
        }

		item = "UL:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*UL) = (double)strtod(position, &endPtr); //double 
        }										

		item = "HR:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*HR) = (double)strtod(position, &endPtr); //double 
        }

		item = "UR:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*UR) = (double)strtod(position, &endPtr); //double 
        }	


    }
	fclose(fp);

	return 1;
}



int read_simulation_setup(char *tempfile, double *simTime, 
	double *Toutput, double *CFL, double *n1 ){		

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;
	long iaux;
	double daux;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {

		item = "simTime:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*simTime) = (double)strtod(position, &endPtr); //double 
        }

		item = "Toutput:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Toutput) = (double)strtod(position, &endPtr); //double 
        }		

		item = "CFL:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*CFL) = (double)strtod(position, &endPtr); //double 
        }

		item = "n1:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*n1) = (double)strtod(position, &endPtr); //double 
        }									

    }
	fclose(fp);

	return 1;
}



int read_boundary_configuration(char *tempfile, 
	double *QIN, double *HIN, 
	double *HOUT, double *ZSOUT){		

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;
	long iaux;
	double daux;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {

		item = "QIN:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*QIN) = (double)strtod(position, &endPtr); //double 
        }

		item = "HIN:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*HIN) = (double)strtod(position, &endPtr); //double 
        }		

		item = "HOUT:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*HOUT) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZSOUT:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZSOUT) = (double)strtod(position, &endPtr); //double 
        }									

    }
	fclose(fp);

	return 1;
}



int h_initialize_variables(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr, 
	double *B, double *P, 
	double *n, double *zb, 	
	double *DU1, double *DU2,
	double *lSW){

    int ic;
	for(ic=0; ic<nCells; ic++){
		A[ic] = 0.0;
		Q[ic] = 0.0;
		h[ic] = 0.0;
		q[ic] = 0.0;
		u[ic] = 0.0;

		c[ic] = 0.0;
		Fr[ic] = 0.0;	

		B[ic] = 0.0;
		P[ic] = 0.0;

		n[ic] = 0.0;
		zb[ic] = 0.0;	

		DU1[ic] = 0.0;	
		DU2[ic] = 0.0;

		lSW[ic] = 0.0;		
	}

    return 1;    

}

int h_set_bed_initial_conditions(int nCells,
	double *x, double *dist, double *zb, double *zr,
	int iniBed,
	double ZBupstream, double ZBslope, 
	double XBstep, double ZBL, double ZBR,
	double XBcoast, double ZBbase, double slopeCoast){

    int ic;
	for(ic=0; ic<nCells; ic++){

		if(iniBed==1){

			zb[ic] = ZBupstream + ZBslope*dist[ic];
			if(zb[ic] < zr[ic]) zb[ic] = zr[ic];

		}else if(iniBed==2){
			
			if(x[ic] < XBstep) {
				zb[ic] = ZBL;
			} else {
				zb[ic] = ZBR;
			}
			if(zb[ic] < zr[ic]) zb[ic] = zr[ic];
		
		}else if(iniBed==3){

			if(x[ic] < XBcoast) {
				zb[ic] = ZBbase;
			} else {
				zb[ic] = ZBbase + slopeCoast*(x[ic]-XBcoast);
			}
			if(zb[ic] < zr[ic]) zb[ic] = zr[ic];

		}		
		
	}

    return 1;

}



int h_set_flow_initial_conditions(int nCells,
	double *x, double *zb, double *h, double *u, 
	int iniCond,
	double Uconst, double Hconst, 
	double Xdam, double ZSL, double ZSR,
	double Xstep, double HL, double UL, double HR, double UR){	

    int ic;
	for(ic=0; ic<nCells; ic++){

		if(iniCond==1){

			h[ic] = Hconst;
			u[ic] = Uconst;
		
		}else if(iniCond==2){

			if(x[ic] < Xdam) {
				h[ic] = ZSL-zb[ic];
			} else {
				h[ic] = ZSR-zb[ic];
			}
			if(h[ic] < 0.0) h[ic] = 0.0;
			u[ic] = 0.0;

		}else if(iniCond==3){

			if(x[ic] < Xstep) {
				h[ic] = HL;
				u[ic] = UL;
			} else {
				h[ic] = HR;
				u[ic] = UR;
			}

		}
		
	}

    return 1;

}



int h_set_inlet_initial_conditions(int nCells,
	double *h, double *u, double *M,
	double QIN, double HIN){	

	double aux1;

	int idx=0;
	
	if(HIN > 0.0) {
		h[idx] = HIN;
	}

	if(QIN > 0.0) {
		aux1 = M[idx]*h[idx];	
		if(aux1 > 0.0) {
			u[idx] = QIN/aux1;	
		}
	}

    return 1;

}



int h_set_outlet_initial_conditions(int nCells,
	double *h, double *u, double *M, double *zb,
	double HOUT, double ZSOUT ){

	int idx=nCells-1;

    if(HOUT > 0.0) { 
        h[idx] = HOUT; 
    }

    if(ZSOUT > 0.0) { 
        h[idx] = ZSOUT-zb[idx]; 
    }	
	if(h[idx]<0.0) h[idx]=0.0;

    return 1;

}



int h_compute_initial_flow_variables(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr, 
	double *B, double *P, 
	double *M, double *zb,
	double *n, double n1){

    int ic;
	for(ic=0; ic<nCells; ic++){		

		if(h[ic] >= tol9){

			//Geometry variables
			B[ic] = M[ic]; 
			P[ic] = M[ic] + 2.*h[ic]; 	

			// Conservative variables			
			A[ic] = B[ic]*h[ic];
			Q[ic] = A[ic]*u[ic];

			c[ic] = sqrt(g*A[ic]/B[ic]);

			// Minimum depth control
			if(h[ic] >= hmin){
				q[ic] = h[ic]*u[ic];			
				Fr[ic] = fabs(u[ic])/c[ic];
			}else{
				Q[ic] = 0.0;
				u[ic] = 0.0;
				q[ic] = 0.0;			
				Fr[ic] = 0.0;		
			}		

		} else {
			B[ic] = M[ic];
			P[ic] = M[ic];

			A[ic] = 0.0;
			Q[ic] = 0.0;

			h[ic] = 0.0;
			c[ic] = 0.0;

			u[ic] = 0.0;
			q[ic] = 0.0;			
			Fr[ic] = 0.0;
		}

		// Cell Manning Coefficient
        n[ic] = n1;

	}			

    return 1;

}   



int h_compute_water_mass(int nCells,
	double *A, double *dx, double *massWt){

    int ic;

	(*massWt) = 0.0;
	for(ic=0; ic<nCells; ic++){

		(*massWt) += A[ic]*dx[ic];		
	
	}

    return 1;

}



int h_compute_flow_time_step(int nWalls, 
	double *A, double *u, double *B, double *dx, 
	double *lSW, double *dtSW){
    
	double uROE, cROE;
	double dt=1e6;

	int ic;
	for(ic=0; ic<nWalls; ic++){		
		
		cROE = sqrt( g*(A[ic]+A[ic+1]) / (B[ic]+B[ic+1]) );
		uROE = (sqrt(A[ic])*u[ic]+sqrt(A[ic+1])*u[ic+1]) / (sqrt(A[ic])+sqrt(A[ic+1]));
		
		lSW[ic]=0.0;
		if(fabs(uROE)+cROE > tol9){
			lSW[ic] = fabs(uROE)+cROE;

			dt = MIN(dt, dx[ic]/lSW[ic] ); 
		}

	}

	(*dtSW) = dt;

    return 1;

}  

int h_compute_wall_fluxes(int nWalls,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr,
	double *B, double *P,
	double *n, double *zb,
  	double *DU1, double *DU2, 
	double *dx){

	double uROE, cROE;                		// ROE averaged variables
	double lambda1, lambda2;				// Eigenvalues
	double e1[3], e2[3];					// Eigenvectors

	double deltaA0, deltaQ0;				// Espatial increments of the conserved variables
	double alpha1, alpha2;					// Wave-strenghts

	double So, Fp, Fw, b1P, b2P;			// Hydrostatic pressure force at bed
	double nav;								// 1D averaged Manning Coefficient
	double Sf, Ff, b1F, b2F;				// Fupwind force at bed
	double beta1, beta2;					// Source term coefficients

	double Astar,Ax,Axx,Qx,Qstart;			// Non-negative area values control - Beta limitation

	double lb1l, lb1r, lb2l, lb2r, lambda1plus, lambda1minus, lambda2plus, lambda2minus;		// Entropy correction

	double Aav, Bav, Pav, Bbav, hav;				// Geometry parameters averaged at the walls

	double aux1,aux2,aux3;	

	int ic;
	for(ic=0; ic<nWalls; ic++){	
		
		if(h[ic] >= tol9 || h[ic+1] >= tol9){			// Wet wall

			// Averaged quantities at the walls
			Aav = (A[ic]+A[ic+1])/2.;
			Bav = (B[ic]+B[ic+1])/2.;
			Pav = (P[ic]+P[ic+1])/2;
			hav = 0.5*(h[ic]+h[ic+1]);
			
			// ROE averaged variables at the walls
			cROE = ;
			uROE = ;
			
			// Eigenvalues 
			lb1l = u[ic]-c[ic];
			lb1r = u[ic+1]-c[ic+1];
			lambda1 = uROE-cROE;
			lb2l = u[ic]+c[ic];
			lb2r = u[ic+1]+c[ic+1];
			lambda2 = uROE+cROE;

			// Eigenvectors
			e1[0] = ;      e1[1] = ;
			e2[0] = ;      e2[1] = ;

			// Wave-strenght coefficients
			deltaA0 = A[ic+1]-A[ic];
			deltaQ0 = Q[ic+1]-Q[ic];
			alpha1 = ;
			alpha2 = ;


			// Source term coefficients
			// Bed slope coefficiente
			So = -(zb[ic+1]-zb[ic])/dx[ic];
			Fp = g*Aav*So*dx[ic]; 	
			b1P = (-1./(2.*cROE)) * (Fp);
			b2P = -b1P;

			// Bed friction coefficient
			nav = 0.5*(n[ic]+n[ic+1]);
			if(friction == 1) {
				// Manning
				if(Fmodel == 1) { Sf = nav*nav*fabs(uROE)*uROE/pow(hav,4./3.); } // Unit formulation
				//Viscous
				if(Fmodel == 2) { Sf = 3.*mu*uROE / (rhow*g*hav); } // Unit formulation
			} else {
				Sf = 0.0;
			}

			Ff = -g*Aav*Sf*dx[ic];

			b1F = (-1./(2.*cROE)) * Ff;
			b2F = -b1F;
			
			// Friction fix
			Qstart = Q[ic] + alpha1*lambda1 - b1P;
			Qx = Qstart - b1F;
			if(Qx*Qstart < 0.0){
			b1F = Qstart;
			b2F = -b1F;
			} 			

			// Global source term coefficients
			beta1 = b1P + b1F;
			beta2 = b2P + b2F;


			// Positivity fix
			Astar = A[ic]+alpha1;		
			if(Astar<tol9) Astar = 0.0;			

			if( (A[ic]>0.0 && A[ic+1]>0.0) &&
				(Astar>0.0) &&
				(lambda1*lambda2 < 0.0) ){

				Ax = A[ic]+alpha1-beta1/lambda1;
				Axx = A[ic+1]-alpha2+beta2/lambda2;
				if(Ax<tol9) Ax = 0.0;
				if(Axx<tol9) Axx = 0.0;

				if(Ax < 0.0){
					beta1 = Astar*lambda1;
					beta2 = -beta1;
				}

				if(Axx < 0.0){
					beta1 = Astar*lambda2;
					beta2 = -beta1;
				}				
			}
			
			// Upddate contributions		
			Ax = A[ic]+alpha1-beta1/lambda1;
			Axx = A[ic+1]-alpha2+beta2/lambda2;
			if(fabs(Ax)<tol12) Ax = 0.0;
			if(fabs(Axx)<tol12) Axx = 0.0;		

			// First wave
			if(A[ic]<tol9 && Ax < 0.0){	// dry-wet wall
				DU1[ic+1] += (lambda1*alpha1-beta1)*e1[0];
				DU2[ic+1] += 0.0;

			} else if(A[ic+1]<tol9 && Axx < 0.0){ // wet-dry wall
				DU1[ic] += (lambda1*alpha1-beta1)*e1[0];
				DU2[ic] += 0.0;

			} else if(lb1l < 0.0 && lb1r > 0.0){
				lambda1minus = lb1l*(lb1r-lambda1)/(lb1r-lb1l);
				DU1[ic] += (lambda1minus*alpha1-beta1)*e1[0];
				DU2[ic] += (lambda1minus*alpha1-beta1)*e1[1];
				lambda1plus = lb1r*(lambda1-lb1l)/(lb1r-lb1l);
	    		DU1[ic+1] += (lambda1plus*alpha1)*e1[0];
				DU2[ic+1] += (lambda1plus*alpha1)*e1[1];

			} else if(lambda1 < 0.0){
				DU1[ic] += (lambda1*alpha1-beta1)*e1[0];
				DU2[ic] += (lambda1*alpha1-beta1)*e1[1];

			} else if(lambda1 >= 0.0){
				DU1[ic+1] += (lambda1*alpha1-beta1)*e1[0];	
				DU2[ic+1] += (lambda1*alpha1-beta1)*e1[1];
			}
			
			// Second wave
			if(A[ic]<tol9 && Ax < 0.0){ 		// dry-wet wall
				DU1[ic+1] += (lambda2*alpha2-beta2)*e2[0];
				DU2[ic+1] += 0.0;

			} else if(A[ic+1]<tol9 && Axx < 0.0){	// wet-dry wall
				DU1[ic] += (lambda2*alpha2-beta2)*e2[0];
				DU2[ic] += 0.0;

			} else if(lb2l < 0.0 && lb2r > 0.0){
				lambda2minus=lb2l*(lb2r-lambda2)/(lb2r-lb2l);
	    		DU1[ic] += (lambda2minus*alpha2)*e2[0];
	    		DU2[ic] += (lambda2minus*alpha2)*e2[1];
				lambda2plus=lb2r*(lambda2-lb2l)/(lb2r-lb2l);
	    		DU1[ic+1] += (lambda2plus*alpha2-beta2)*e2[0];
	    		DU2[ic+1] += (lambda2plus*alpha2-beta2)*e2[1];

			} else if(lambda2 < 0.0){
				DU1[ic] += (lambda2*alpha2-beta2)*e2[0];
				DU2[ic] += (lambda2*alpha2-beta2)*e2[1];

			} else if(lambda2 >= 0.0){
				DU1[ic+1] += (lambda2*alpha2-beta2)*e2[0];
				DU2[ic+1] += (lambda2*alpha2-beta2)*e2[1];
			}

		} // End wet walls

	} // End of fluxes calculation - Wall loop				

    return 1;

}


int h_check_depth_positivity(int nCells, 
	double *A, double *DU1, double *dx, double *dt){
    
	double dt0;
	double aux1;

	dt0 = (double) (*dt);
	
	int ic;
	for(ic=0; ic<nCells; ic++){		
		
		aux1 = A[ic] - (*dt)*DU1[ic]/dx[ic];		
		if(fabs(aux1) < tol9) aux1 = 0.0;

		if( aux1 < 0.0) {
			while(aux1 < 0.0 && (*dt) > 0.01*dt0) {
				(*dt) = 0.90 * (*dt);

				aux1 = A[ic] - (*dt)*DU1[ic]/dx[ic];		
				if(fabs(aux1) < tol9) aux1 = 0.0;
			}
		}

	}			

    return 1;

}  

int h_update_cells(int nCells, 
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr,
	double *B, double *P,
	double *DU1, double *DU2, 
	double *M, double *dx,
  	double dt){
	
	int ic;
	for(ic=0; ic<nCells; ic++){		
		
		A[ic] = A[ic] - DU1[ic]*dt/dx[ic];
		Q[ic] = Q[ic] - DU2[ic]*dt/dx[ic];

		// Geometric variable actualization
		B[ic] = M[ic];
		P[ic] = M[ic] + 2.*h[ic];
		
		// Flow variables
		h[ic] = A[ic] / B[ic];		
		if(h[ic] >= tol9){
			c[ic] = sqrt(g*A[ic]/B[ic]);

			// Minimum depth control
			if(h[ic] >= hmin){
				u[ic] = Q[ic]/A[ic];
				q[ic] = h[ic]*u[ic];			
				Fr[ic] = fabs(u[ic])/c[ic];
			}else{
				Q[ic] = 0.0;
				u[ic] = 0.0;
				q[ic] = 0.0;			
				Fr[ic] = 0.0;		
			}				
		} else {
			A[ic] = 0.0;
			Q[ic] = 0.0;
			B[ic] = M[ic];
			P[ic] = M[ic];
			h[ic] = 0.0;
			c[ic] = 0.0;
			u[ic] = 0.0;
			q[ic] = 0.0;			
			Fr[ic] = 0.0;
		}

		// Reset U fluxes differences
		DU1[ic]=0.0;
		DU2[ic]=0.0;
		
	}			

    return 1;

}  

int h_set_inlet_conditions(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr,
	double *B, double *P, 
	double *M,
	double QIN, double HIN){

	int idx=0;

	if(QIN > 0.0){
		Q[idx] = QIN;
	}

	if(HIN > 0.0) {
		h[idx] = HIN;

		B[idx] = M[idx];
		P[idx] = B[idx] + 2.*h[idx];

		A[idx] = B[idx]*h[idx];
		c[idx] = sqrt(g*A[idx]/B[idx]);
	}	

	u[idx] = Q[idx]/A[idx];	
		
	Fr[idx] = fabs(u[idx])/c[idx];
	q[idx] = h[idx]*u[idx];

    return 1;

}

int h_set_outlet_conditions(int nCells,
	double *A, double *Q, double *h, double *q, double *u, 
	double *c, double *Fr,
	double *B, double *P, 
	double *M, double *zb,
	double HOUT, double ZSOUT){

	int idx=nCells-1;

	if(HOUT > 0.0) {
		h[idx] = HOUT;

		B[idx] = M[idx];
		P[idx] = B[idx] + 2.*h[idx];

		A[idx] = B[idx]*h[idx];
		c[idx] = sqrt(g*A[idx]/B[idx]);
	}

	if(ZSOUT > 0.0) {
		h[idx] = ZSOUT-zb[idx];

		B[idx] = M[idx];
		P[idx] = B[idx] + 2.*h[idx];

		A[idx] = B[idx]*h[idx];
		c[idx] = sqrt(g*A[idx]/B[idx]);
	}

	u[idx] = Q[idx]/A[idx];			
	Fr[idx] = fabs(u[idx])/c[idx];
	q[idx] = h[idx]*u[idx];

    return 1;

}














