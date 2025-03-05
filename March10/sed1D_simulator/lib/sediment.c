/* Mud flow + Exner 1D Simulator
  ===========================================
   Weakly Coupled Model
  ===========================================================================================
   Sergio Martinez-Aranda    -    Dic 2020
  =========================================================================================== */

#include "sediment.h"


int read_bedload_configuration(char *tempfile, int *bedloadCap,
	double *rhos, double *poros, double *ds, 
	double *shieldsC, double *nManSed,
	double *Ggrass, double *BetaB){		

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

		item = "rhos:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*rhos) = (double)strtod(position, &endPtr); //double 
        }

		item = "poros:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*poros) = (double)strtod(position, &endPtr); //double 
        }		

		item = "ds:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ds) = (double)strtod(position, &endPtr); //double 
        }

		item = "shieldsC:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*shieldsC) = (double)strtod(position, &endPtr); //double 
        }

		item = "nManSed:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*nManSed) = (double)strtod(position, &endPtr); //double 
        }	

		item = "bedloadCap:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*bedloadCap) = (int)strtol(position, &endPtr, 10); //integer
        }

		item = "Ggrass:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*Ggrass) = (double)strtod(position, &endPtr); //double 
        }			

		item = "BetaB:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*BetaB) = (double)strtod(position, &endPtr); //double 
        }	

    }
	fclose(fp);

	return 1;
}



int read_bedload_boundary_configuration(char *tempfile,
	double *QSIN, double *ZBIN){		

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

		item = "QSIN:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*QSIN) = (double)strtod(position, &endPtr); //double 
        }

		item = "ZBIN:";
        position = strstr(line, item);
        if(position){  
			position += strlen(item);          
			(*ZBIN) = (double)strtod(position, &endPtr); //double 
        }					

    }
	fclose(fp);

	return 1;
}



int h_initialize_bedload_variables( int nCells, int nWalls,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *DQs, double *lB){
    
	int ic;
	for(ic=0; ic<nCells; ic++){		
		Ab[ic] = 0.0;
		Qs[ic] = 0.0;

		tauB[ic] = 0.0;
		Gcell[ic] = 0.0;	
		qs[ic] = 0.0;

		DQs[ic] = 0.0;
	}

	for(ic=0; ic<nWalls; ic++){		
		lB[ic] = 0.0;
	}	

    return 1;    

}



int h_compute_initial_bedload_variables( int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr){

	double theta;
	double aux1, aux2, aux3;
   
    int ic;
	for(ic=0; ic<nCells; ic++){		
		
		// Bed variables
		if((zb[ic]-zr[ic]) > tol9 ){
			Ab[ic] = M[ic]*(zb[ic]-zr[ic]);	
		}else{
        	Ab[ic] = 0.0;
		}			
		
		// Cell Manning Coefficient
		if((zb[ic]-zr[ic]) > tol9 ){
			n[ic] = nManSed;
		}else{
        	n[ic] = n1;
		}		

		if(h[ic] >= tol9){  //wet cell

			// Shear stress at the bottom
			tauB[ic] = rhow * g * n[ic]*n[ic] * fabs(u[ic]*u[ic]) / cbrt(h[ic]); 

			// Shield parameter
			theta = tauB[ic] / (g*(rhos-rhow)*ds);
			
			// Bedload rate
			if(bedloadCap==0){
				Gcell[ic] = 0.0;
			}else if(bedloadCap==1){	
				Gcell[ic] = Ggrass;
			}else if(bedloadCap==2){
				aux1 = theta-shieldsC;
				if(theta>0.0 && aux1>0.0){
					aux2 = 1.-shieldsC/theta;
					aux3 = n[ic]*n[ic]*n[ic]*sqrt(g) / ((rhos/rhow-1.)*sqrt(h[ic]));
					Gcell[ic] = BetaB * 8.*aux2*sqrt(aux2) * aux3;
				}else{
					Gcell[ic]=0.0;
				}
			}
			if((zb[ic]-zr[ic]) < tol6 ){
				Gcell[ic]=0.0;
			}
			qs[ic] = Gcell[ic] * (u[ic]*u[ic]) * u[ic];  					
			Qs[ic] = M[ic] * qs[ic];
		
		}else{ //dry cell
		
			tauB[ic] = 0.0;
			Gcell[ic] = 0.0;
			qs[ic] = 0.0;  					
			Qs[ic] = 0.0;
		
		}	

	}			

    return 1;

} 



int h_set_bedload_inlet_initial_conditions(int nCells,
	double *h, double *M, double *Qs, double *qs,
	double QSIN){

	int ic=0;	

	if(QSIN > 0.0) {
		if(h[ic] >= tol9){
			qs[ic] = QSIN/M[ic];
		}else{
			qs[ic] = 0.0;
		}
		Qs[ic] = qs[ic]*M[ic];
	}	

    return 1;

}



int h_compute_bedload_mass(int nCells,
	double *Ab, double *dx, double *massBt){

    int ic;

	(*massBt) = 0.0;
	for(ic=0; ic<nCells; ic++){	
		(*massBt) += Ab[ic]*dx[ic];
	}

    return 1;

}



int h_compute_bedload_time_step(int nWalls,
	double *A, double *u, double *zb, double *qs, 
	double *dx, 
	double ds, double poros,
	double *lB, double *dtB){
    
	double uROE;
	double dt=1e6;
	double xi=1./(1.-poros);
	
	int ic;
	for(ic=0; ic<nWalls; ic++){		
		//**** complete code here ********
		//storage the bedlead virtual celerity at each wall in the lB array
		//reduce the minimum bedload time step in the local dt scalar




		//**** end of completed code ********	
	}	

	(*dtB) = dt;		

    return 1;

}   



int h_compute_bedload_wall_fluxes(int nWalls, 
    double *h, double *u, double dt,
    double *Ab, double *zb, double *Qs, double *qs, 
    double *lB, double *DQs,
	double ds, double poros,
    double *dx, int *solidWall){

	double lambdaB;							// Sediment flux wave speed
	double QSx;
	double maxSedVol, eroSedVol;  

    double	xi = 1./(1.-poros);  
	double aux1,aux2,aux3;	

	int ic;
	for(ic=0; ic<nWalls; ic++){	
		
		if( (h[ic] >= ds && h[ic+1] >= ds) &&
			(solidWall[ic]==0) ){	

			// Numerical sediment flux
			//**** complete code here ********
			//compute the numerical bedload rate at the wall in the local Qsx scalar
			//remember that the bedload virtual celerity at each wall was storaged in the lB array




			//**** end of completed code ********			
            
            // Maximum erosion fix for numerical bedload rate
			maxSedVol=0.0;
            if(QSx > 0.0){
				maxSedVol = Ab[ic]*0.5*dx[ic];
			}
            if(QSx < 0.0) {
				maxSedVol = Ab[ic+1]*0.5*dx[ic+1];
			}
            eroSedVol = xi*fabs(QSx)*dt;

            if(eroSedVol > maxSedVol) {
                aux1 = QSx/fabs(QSx);
                QSx = aux1*maxSedVol/(xi*dt);
            }	

			//Cell fluxes
			//**** complete code here ********
			//add the bedload contribution at each cell in the DQs array
			//remember that both left and right contributions must be added to the cells ic and ic+1
			//not that the porosity effect must be included in this step




			//**** end of completed code ********									

            
            // Boundary ghost cells
            if(ic==0){
                QSx = Qs[ic];
                if(QSx < 0.0){
                    maxSedVol = Ab[ic]*0.5*dx[ic];
                    eroSedVol = xi*fabs(QSx)*dt;

                    if(eroSedVol > maxSedVol) {
                        aux1 = QSx/fabs(QSx);
                        QSx = aux1*maxSedVol/(xi*dt);
                    }
                }

                DQs[ic] += -1.*xi*QSx;
            }

            if(ic==nWalls-1){
                QSx = Qs[ic+1];
                if(QSx > 0.0){
                    maxSedVol = Ab[ic+1]*0.5*dx[ic+1];
                    eroSedVol = xi*fabs(QSx)*dt;

                    if(eroSedVol > maxSedVol) {
                        aux1 = QSx/fabs(QSx);
                        QSx = aux1*maxSedVol/(xi*dt);
                    }
                }

                DQs[ic+1] += xi*QSx;
            }				
                

		} // End wet walls

	} // End of fluxes calculation - Wall loop				

    return 1;

}



int h_update_bedload_cells(int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr,
    double *DQs,
	double *dx, double dt){

	double theta;
	double aux1, aux2, aux3;
	
	int ic;
	for(ic=0; ic<nCells; ic++){		
		
		//**** complete code here ********
		//update de bed area Ab for each cell
		//update de bed elevation zb at each cell
		//remember that the bedload contribution at each cell was added in the DQs array




		//**** end of completed code ********			

		
		// Cell Manning Coefficient
		if((zb[ic]-zr[ic]) > tol9 ){
			n[ic] = nManSed;
		}else{
        	n[ic] = n1;
		}	


		if(h[ic] >= tol9){  //wet cell

			// Shear stress at the bottom
			//**** complete code here ********
			//compute the bed shear stress tauB at each cell 



			//**** end of completed code ********				

			// Shield parameter
			theta = tauB[ic] / (g*(rhos-rhow)*ds);
			
			// Bedload rate
			//**** complete code here ********
			//compute the unit bedload rate qs at each cell 
			//remember that the bedload computation depends on the selected bedloadCap model
			//grass-type formulation is recomended by computing the Gcell factor for each cell




			//**** end of completed code ********
			if((zb[idx]-zr[idx]) < tol6 ){
				Gcell[idx]=0.0;
			}		
			qs[idx] = Gcell[idx] * (u[idx]*u[idx]) * u[idx]; 									
								
			Qs[ic] = M[ic] * qs[ic];
		
		}else{ //dry cell
		
			tauB[ic] = 0.0;
			Gcell[ic] = 0.0;
			qs[ic] = 0.0;  					
			Qs[ic] = 0.0;
		
		}	


        // Reset bedload contributions
		//**** complete code here ********
		//re-initialized the bedload contribution DQs array for the next time step 



		//**** end of completed code ********	

	}			

    return 1;

} 



int h_set_bedload_inlet_conditions( int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr,
    double dt, 
	double QSIN, double ZBIN){

	double addBLrate;
	double theta;
    double aux1, aux2, aux3;

	int idx=0;

	addBLrate = 0.0;
	if(ZBIN > 0.0) {		
		aux1=Ab[idx];

		if(ZBIN < zr[idx]) {
			ZBIN < zr[idx];
		}

		zb[idx] = ZBIN;
		Ab[idx] = M[idx]*(zb[idx]-zr[idx]);

		addBLrate = (1.-poros)*(Ab[idx]-aux1); //Bedload solid mass
		addBLrate /= dt; //Total bedload rate
		addBLrate /= M[idx]; //Unit bedload discharge
	}

    // Cell Manning Coefficient
	if((zb[idx]-zr[idx]) > tol9 ){
		n[idx] = nManSed;
	}else{
		n[idx] = n1;
	}	

	if(h[idx] >= tol9){  //wet cell

		// Shear stress at the bottom
		tauB[idx] = rhow * g * n[idx]*n[idx] * fabs(u[idx]*u[idx]) / cbrt(h[idx]); 

		// Shield parameter
		theta = tauB[idx] / (g*(rhos-rhow)*ds);
		
		// Bedload rate
		if(bedloadCap==0){
			Gcell[idx] = 0.0;
		}else if(bedloadCap==1){	
			Gcell[idx] = Ggrass;
		}else if(bedloadCap==2){
			aux1 = theta-shieldsC;
			if(theta>0.0 && aux1>0.0){
				aux2 = 1.-shieldsC/theta;
				aux3 = n[idx]*n[idx]*n[idx]*sqrt(g) / ((rhos/rhow-1.)*sqrt(h[idx]));
				Gcell[idx] = BetaB * 8.*aux2*sqrt(aux2) * aux3;
			}else{
				Gcell[idx]=0.0;
			}
		}
		if((zb[idx]-zr[idx]) < tol6 ){
			Gcell[idx]=0.0;
		}		
		qs[idx] = Gcell[idx] * (u[idx]*u[idx]) * u[idx]; 
		
		if(QSIN > 0.0) {
			qs[idx] = QSIN / M[idx];
		} 
		qs[idx] += addBLrate;
							
		Qs[idx] = M[idx] * qs[idx];

	}else{ //dry cell
	
		tauB[idx] = 0.0;
		Gcell[idx] = 0.0;
		qs[idx] = 0.0;  					
		Qs[idx] = 0.0;
	
	}	

    return 1;
}



int h_set_bedload_outlet_conditions( int nCells,
	int bedloadCap, double Ggrass, double BetaB,
	double rhos, double poros, double ds, 
	double shieldsC, double nManSed, double n1,
	double *Ab, double *Qs, 
	double *tauB, double *Gcell, double *qs, 
	double *h, double *u, 
	double *n, double *M, double *zb, double *zr){

	double theta;
    double aux1, aux2, aux3;

	int idx=nCells-1;;

    // Cell Manning Coefficient
	if((zb[idx]-zr[idx]) > tol9 ){
		n[idx] = nManSed;
	}else{
		n[idx] = n1;
	}	

	if(h[idx] >= tol9){  //wet cell

		// Shear stress at the bottom
		tauB[idx] = rhow * g * n[idx]*n[idx] * fabs(u[idx]*u[idx]) / cbrt(h[idx]); 

		// Shield parameter
		theta = tauB[idx] / (g*(rhos-rhow)*ds);
		
		// Bedload rate
		if(bedloadCap==0){
			Gcell[idx] = 0.0;
		}else if(bedloadCap==1){	
			Gcell[idx] = Ggrass;
		}else if(bedloadCap==2){
			aux1 = theta-shieldsC;
			if(theta>0.0 && aux1>0.0){
				aux2 = 1.-shieldsC/theta;
				aux3 = n[idx]*n[idx]*n[idx]*sqrt(g) / ((rhos/rhow-1.)*sqrt(h[idx]));
				Gcell[idx] = BetaB * 8.*aux2*sqrt(aux2) * aux3;
			}else{
				Gcell[idx]=0.0;
			}
		}
		if((zb[idx]-zr[idx]) < tol6 ){
			Gcell[idx]=0.0;
		}
		qs[idx] = Gcell[idx] * (u[idx]*u[idx]) * u[idx]; 							
		Qs[idx] = M[idx] * qs[idx];

	}else{ //dry cell
	
		tauB[idx] = 0.0;
		Gcell[idx] = 0.0;
		qs[idx] = 0.0;  					
		Qs[idx] = 0.0;
	
	}	

    return 1;
}




