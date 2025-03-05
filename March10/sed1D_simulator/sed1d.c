//   Shallow Water Equations 1D Simulator (SWE1D)
//  ===========================================================================================
//   Rectangular cross section
//  ===========================================================================================
//   Version 1.0 - Jan 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 


#include "define.h"
#include "lib/shallow_water.h"
#include "lib/sediment.h"

// Inputs and storage
char filename[1024];
FILE *fp;
FILE *logFile;
FILE *qFile;




///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


// MAIN CODE FUNTION
main(){


// OPENING DATA STORAGE
	int i; 
	int ic;
	char tempfile[1024], msg[1024];
	FILE *fp;

	//i=system ("rm outputFiles/*.out");	// Remove old files (linux only)

	logFile=fopen("outputFiles/log.out","w");
	fclose(logFile);

	qFile=fopen("outputFiles/discharge.out","w");
	fclose(qFile);
	
	#if DISPLAY
	Gnuplot gp;			// Graphical output pipe (linux only)
	#endif



// SIMULATOR HEADERS
	printf("\n   SHALLOW WATER EQUATIONS  <<<>>> 1D SIMULATOR");
	printf("\n   Version 1.0 - Jan 2025");
	printf("\n   Computational Hydraulics Group - University of Zaragoza");
	printf("\n   --------------------------------------------------");

	printf("\n\n>> Case:	");
	printf("\n");

	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n   SHALLOW WATER EQUATIONS <<<>>> 1D SIMULATOR");
	fprintf(logFile,"\n   Version 1.0 - Jan 2025");
	fprintf(logFile,"\n   Computational Hydraulics Group - University of Zaragoza");
	fprintf(logFile,"\n   --------------------------------------------------");

	fprintf(logFile,"\n\n>> Case:	");
	fprintf(logFile,"\n");
	fclose(logFile);



// GEOMETRY DATA
	// Load geometry data	
	int nCells, nWalls;
	double DeltaX, Xmin, Width, Length;

	sprintf(tempfile,"input/geometry.input");
	read_geometry_file(tempfile, &(nCells), &(nWalls), 
		&(DeltaX), &(Xmin), &(Width), &(Length));

	printf("\n\n>> Geometry: ");
	printf("\n     nCells %d - nWalls %d ", nCells, nWalls);
	printf("\n     DeltaX %lf - Xmin %lf - Width %lf - Length %lf", DeltaX, Xmin, Width, Length);			
	printf("\n");

	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n\n>> Geometry: ");
	fprintf(logFile,"\n     nCells %d - nWalls %d ", nCells, nWalls);
	fprintf(logFile,"\n     DeltaX %lf - Xmin %lf - Width %lf - Length %lf", DeltaX, Xmin, Width, Length);			
	fprintf(logFile,"\n");
	fclose(logFile);	


	// Channel geometry
	double *x, *dx, *dist;
	double *M, *zr; 
	x = (double*) malloc( nCells*sizeof(double) );
	dx = (double*) malloc( nCells*sizeof(double) );
	dist = (double*) malloc( nCells*sizeof(double) );
	M = (double*) malloc( nCells*sizeof(double) );
	zr = (double*) malloc( nCells*sizeof(double) );

	for(ic=0; ic<nCells; ic++) {
		dx[ic] = DeltaX;
		M[ic] = Width;
		x[ic] = Xmin + ic*dx[ic];
		dist[ic] = x[ic] - Xmin;
		zr[ic] = ZRlevel;
	}


	// Bed configuration
	int iniBed;
	double ZBupstream, ZBslope;
	double XBstep, ZBL, ZBR; 
	double XBcoast, ZBbase, slopeCoast;

	sprintf(tempfile,"input/geometry.input");
	read_bed_configuration(tempfile, &(iniBed), 
		&(ZBupstream), &(ZBslope), 
		&(XBstep), &(ZBL), &(ZBR),
		&(XBcoast), &(ZBbase), &(slopeCoast));


	// Flow configuration
	int iniCond;
	double Uconst, Hconst;
	double Xdam, ZSL, ZSR;
	double Xstep, HL, UL, HR, UR;	

	sprintf(tempfile,"input/geometry.input");
	read_flow_configuration(tempfile, &(iniCond), 
		&(Uconst), &(Hconst), 
		&(Xdam), &(ZSL), &(ZSR),
		&(Xstep), &(HL), &(UL), &(HR), &(UR));


	printf("\n\n>> Geometry loaded");
	printf("\n");

	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n\n>> Geometry loaded");
	fprintf(logFile,"\n");
	fclose(logFile);



// MAIN VARIABLE DEFINITION

	// Declare variables
	clock_t CPUtime;

	int iter;										// Iteraction index
	int nout;										// Output index
	
	double *A, *Q;									// Conservative variables
	
	double *h, *q;									// Primitive variables
	double *c, *u;									// Wave speed
	double *Fr;										// Froude number - Se calcula solo al final de cada iteracion	

	double *B;										// Width at the free surface
	double *P;										// Wet perimeter
	double *n;										// 1D averaged Manning Coefficient	
	double *zb;										// Bed surface level	
	
	double *DU1;									// Variation Conserved variable A in Dt
	double *DU2;									// Variation Conserved variable Q in Dt

	int *solidWall;

	double t, dt, dt0;
	double dtmin, dtaux, dtSW;
	double *lSW;
	double condCFL;									// Efective CFL

	// Mass error monitor
	double massWt0;									// Initial mass t = t0
	double massWtn;									// Mass at the time t = tn
	double Qin;										// Discharge at intlet for the time t = t0
	double Qout;									// Discharge at outlet for the time t = t0
	double massWerror;
	double Qbalance;

	#if BEDLOAD
	double *Ab, *Qs;
	double *tauB, *Gcell, *qs;
	double *DQs;

	double dtB;
	double *lB;

	double massBt0;									// Initial bed mass t = t0
	double massBtn;									// Bed mass at the time t = tn
	double QSin;									// Bedload rate at intlet for the time t = t0
	double QSout;									// Bedload rate at outlet for the time t = t0
	double massBerror;
	double QSbalance;
	#endif

	// Auxiliars
	double aux1,aux2,aux3;
	double xcurve, ycurve1, ycurve2, ycurve3;
	std::vector<std::vector<double>> curve1;
	std::vector<std::vector<double>> curve2; 
	std::vector<std::vector<double>> curve3;



	// Allocate memory 
	A = (double*) malloc( nCells* sizeof(double) );
	Q = (double*) malloc( nCells* sizeof(double) );
	h = (double*) malloc( nCells* sizeof(double) );
	q = (double*) malloc( nCells* sizeof(double) );
	u = (double*) malloc( nCells* sizeof(double) );

	c = (double*) malloc( nCells* sizeof(double) );
	Fr = (double*) malloc( nCells* sizeof(double) );
	
	B = (double*) malloc( nCells* sizeof(double) );
	P = (double*) malloc( nCells* sizeof(double) );
	
	n = (double*) malloc( nCells* sizeof(double) );	
	zb = (double*) malloc( nCells* sizeof(double) );

	DU1 = (double*) malloc( nCells* sizeof(double) );
	DU2 = (double*) malloc( nCells* sizeof(double) );

	lSW = (double*) malloc( nWalls* sizeof(double) );	

	solidWall = (int*) malloc( nWalls* sizeof(int) );

	#if BEDLOAD
	Ab = (double*) malloc( nCells* sizeof(double) );
	Qs = (double*) malloc( nCells* sizeof(double) );

	tauB = (double*) malloc( nCells* sizeof(double) );
	Gcell = (double*) malloc( nCells* sizeof(double) );
	qs = (double*) malloc( nCells* sizeof(double) );

	DQs = (double*) malloc( nCells* sizeof(double) );

	lB = (double*) malloc( nWalls* sizeof(double) );
	#endif





// SIMULATION SETUP
	double simTime;
	double Toutput;
	double CFL;
	double n1;

	sprintf(tempfile,"input/simulation.input");
	read_simulation_setup(tempfile, &(simTime), 
		&(Toutput), &(CFL), &(n1));


	double QIN, HIN;
	double HOUT, ZSOUT;
	QIN = -1.0;
	HIN = -1.0;
	HOUT = -1.0;
	ZSOUT = -1.0;

	sprintf(tempfile,"input/simulation.input");
	read_boundary_configuration(tempfile, 
		&(QIN), &(HIN), 
		&(HOUT), &(ZSOUT));	


	#if BEDLOAD
	int bedloadCap;
	double rhos; 
	double poros;
	double ds;
	double shieldsC;
	double nManSed;	
	double Ggrass;
	double BetaB;

	sprintf(tempfile,"input/sediment.input");
	read_bedload_configuration(tempfile, &(bedloadCap), 
		&(rhos), &(poros), &(ds), 
		&(shieldsC), &(nManSed),
		&(Ggrass),&(BetaB));

	double QSIN;
	double ZBIN;	
	QSIN = -1.0;
	ZBIN = -1.0;

	sprintf(tempfile,"input/sediment.input");
	read_bedload_boundary_configuration(tempfile,
		&(QSIN), &(ZBIN));		

	#endif	



	printf("\n   --------------------------------------------------");
	printf("\n>> Simulation setup loaded:");
	printf("\n     Simulation time: %.1lf s",simTime);
	printf("\n     CFL: %.2lf",CFL);
	printf("\n     Friction term activated: %d",friction);
	printf("\n     Friction model: %d",Fmodel);
	printf("\n     Number of cells: %d",nCells);
	printf("\n     dx: %.2lf m", dx[0]);
	printf("\n     Channel: Length %.2lf m - Width %.2lf",Length,Width);
	printf("\n     Data saved each %d iterations",Niter);

	#if BEDLOAD
	printf("\n");
	printf("\n     Sediment density: %04.1lf kg/m3",rhos);
	printf("\n     Sediment diameter: %.3lf mm", 1000.*ds);
	printf("\n     Bed layer porosity: %.3lf",poros);
	#endif

	printf("\n\n>> Simulation setup loaded");
	printf("\n");


	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n\n   --------------------------------------------------");
	fprintf(logFile,"\n>> Simulation setup loaded:");
	fprintf(logFile,"\n     Simulation time: %.1lf s",simTime);
	fprintf(logFile,"\n     CFL: %.2lf",CFL);
	fprintf(logFile,"\n     Friction term activated: %d",friction);
	fprintf(logFile,"\n     Friction model: %d",Fmodel);
	fprintf(logFile,"\n     Number of cells: %d",nCells);
	fprintf(logFile,"\n     dx: %.2lf m",dx[0]);
	fprintf(logFile,"\n     Channel: Length %.2lf m - Width %.2lf",Length,Width);
	fprintf(logFile,"\n     Data saved each %d iterations",Niter);

	#if BEDLOAD
	fprintf(logFile,"\n");
	fprintf(logFile,"\n     Sediment density: %04.1lf kg/m3",rhos);
	fprintf(logFile,"\n     Sediment diameter: %.3lf mm", 1000.*ds);
	fprintf(logFile,"\n     Bed layer porosity: %.3lf",poros);
	#endif	
 
	fprintf(logFile,"\n\n>> Simulation setup loaded");
	fprintf(logFile,"\n");
	fclose(logFile);	



// VARIABLE INITIALIZATION
	iter = 1;					// Iteration number
	t = 0.0;					// Initial time

	h_initialize_variables( nCells, nWalls,
		A,  Q,  h,  q,  u, 
		c,  Fr, 
		B,  P,
		n, zb, 
		DU1,  DU2,
		lSW, solidWall);

	#if BEDLOAD
	h_initialize_bedload_variables( nCells, nWalls,
		Ab,  Qs, 
		tauB,  Gcell, qs, 
		DQs, lB);	
	#endif	

	printf("\n\n>> Flow initialization completed");

	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n\n>> Flow initialization completed");
	fclose(logFile);








/////////////////////////////////////////////////////////////////////////
// INITIAL CONDITIONS 
/////////////////////////////////////////////////////////////////////////

	// Bed initial condition
	h_set_bed_initial_conditions(nCells,
		x, dist, zb, zr,
		iniBed,
		ZBupstream, ZBslope, 
		XBstep, ZBL, ZBR,
		XBcoast, ZBbase, slopeCoast);

	// Flow initial condition
	h_set_flow_initial_conditions( nCells,
		x, zb, h, u,
		iniCond,
		Uconst, Hconst, 
		Xdam, ZSL, ZSR,
		Xstep, HL, UL, HR, UR);

	// Inlet flow condition
	h_set_inlet_initial_conditions( nCells,
		h, u, M,
		QIN, HIN);
	
	// Outlet water depth
	h_set_outlet_initial_conditions( nCells,
		h, u, M, zb,
		HOUT, ZSOUT);

	// Flow variables calculation
	h_compute_initial_flow_variables( nCells,
		A,  Q,  h,  q,  u, 
		c,  Fr, 
		B,  P, 
		M, zb, 
		n, n1);

	#if BEDLOAD
	// Bedload variables calculation
	h_compute_initial_bedload_variables( nCells,
		bedloadCap, Ggrass, BetaB,
		rhos, poros, ds, 
		shieldsC, nManSed, n1,
		Ab, Qs, 
		tauB, Gcell, qs, 
		h, u, 
		n, M, zb, zr);

	// Inlet bedload condition
	h_set_bedload_inlet_initial_conditions( nCells,
		h, M, Qs, qs,
		QSIN);			
	#endif	
		
	printf("\n\n>> Initial conditions loaded");

	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n\n>> Initial conditions loaded");
	fclose(logFile);




	// Initial mass balance
	h_compute_water_mass( nCells,
		A,  dx,  &(massWtn));
	Qin = Q[0];	
	Qout = Q[nCells-1];
	Qbalance = 0.0;

	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n>> INITIAL FLOW MASS %.6lf m3",massWtn);	
	fprintf(logFile,"\n");
	fclose(logFile);

	#if BEDLOAD
	h_compute_bedload_mass( nCells,
		Ab,  dx,  &(massBtn));
	QSin = Qs[0];	
	QSout = Qs[nCells-1];
	QSbalance = 0.0;


	logFile=fopen("outputFiles/log.out","a+");
	fprintf(logFile,"\n>> INITIAL BED MASS %.6lf m3",massWtn);	
	fprintf(logFile,"\n");
	fclose(logFile);
	#endif			




	// Store initial conditions
	//Cell data
	nout = 0;
	sprintf(filename, "outputFiles/celldata%d.out",nout);
	fp = fopen(filename,"w");
	for(ic=0; ic<nCells; ic++){
		fprintf(fp,"%.3lf\t %.6lf\t %.6lf\t %.6lf",
			x[ic],
			zb[ic],
			h[ic],
			u[ic]);
		#if BEDLOAD
		fprintf(fp,"\t %.6lf",
			qs[ic]);
		#endif
		fprintf(fp," \n");
	}
	fclose(fp);

	// Discharge data
	qFile=fopen("outputFiles/discharge.out","a+");	
	fprintf(qFile,"%.3lf\t %.6lf\t %.6lf\t %.6lf\t %.6lf",
		t,
		massWtn,
		Qin,Qout,Qbalance);
	#if BEDLOAD
	fprintf(qFile,"\t %.6lf\t %.6lf\t %.6lf\t %.6lf",
		massBtn,
		QSin,QSout,QSbalance);	
	#endif
	fprintf(qFile," \n");
	fclose(qFile);

	nout += 1;		

	printf("\n\nPress INTRO key to start ...");
	printf("\n\n");
	getchar();







CPUtime = clock();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////     TIME LOOP     /////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
while(t <= simTime) {	

// TIME STEP COMPUTATION
	h_compute_flow_time_step( nWalls,
		A,  u,  B,  dx,
		lSW,  &(dtSW) );
	
	dtmin = dtSW;

	#if BEDLOAD
	h_compute_bedload_time_step( nWalls,
		A,  u,  zb,  qs,
		dx,
		ds, poros,
		lB,  &(dtB) );	

	dtmin = MIN(dtSW, dtB);	
	#endif

	dt = CFL * dtmin;



// DISPLAY TIME AND TIME STEP
	if(iter%Niter == 0){
		printf("\n\n--------------------------------------------------------------------------------------------------------");
		printf("\n>> Time: %.3lf seg  -  Iteration %d  -  Time step %.6lf",t,iter,dt);
	}


// FLUXES CALCULATION
	h_compute_wall_fluxes( nWalls,
		A, Q, h, q, u, 
		c, Fr,
		B, P,
		n, zb,
		DU1, DU2, 
		dx, solidWall);

	#if BEDLOAD
	h_compute_bedload_wall_fluxes(nWalls, 
		h, u, dt,
		Ab, zb, Qs, qs, 
		lB, DQs,
		ds, poros,
		dx, solidWall);
	#endif	

// NEGATIVE WETTED AREA CHECKING
	h_check_depth_posotivity( nCells,
		A,  DU1,  dx,  &(dt));

	  
// SHALLOW WATER CELL ACTUALIZATION
	h_update_cells( nCells,
		A, Q, h, q, u, 
		c, Fr,
		B, P,
		DU1, DU2, 
		M, dx,
		dt);

	h_update_wetdry_fronts(nWalls, 
		h, zb,
		Q, q, u, Fr,
		solidWall);	

	#if BEDLOAD
	h_update_bedload_cells( nCells,
		bedloadCap, Ggrass, BetaB,
		rhos, poros, ds, 
		shieldsC, nManSed, n1,
		Ab, Qs, 
		tauB, Gcell, qs, 
		h, u, 
		n, M, zb, zr,
		DQs,
		dx, dt);
	#endif	
	

// SIMULATION MONITORS
	massWt0 = massWtn;
	h_compute_water_mass( nCells,
		A,  dx,  &(massWtn));

	//Compute mass error
	if(massWt0 > 1e-16) { 
		massWerror = (massWtn - (massWt0-(Qout-Qin)*dt)) / massWt0; 
		if(fabs(massWerror) < 1e-16) massWerror = 1e-16;
	} else { 
		massWerror = 0.0; 
	}

	Qbalance += (Qin-Qout)*dt;

	if(iter%Niter == 0){
		printf("\n\tMass error %.3e\t",massWerror);
		printf("\n\tWater discharge IN %.6lf OUT %.6lf  [m3/s]",Qin,Qout);			
	}


	#if BEDLOAD
	massBt0 = massBtn;
	h_compute_bedload_mass( nCells,
		Ab,  dx,  &(massBtn));

	//Compute bed mass error
	if(massBt0 > 1e-16) { 
		massBerror = (massBtn - (massBt0 - (1./(1.-poros))*(QSout-QSin)*dt)) / massBt0; 
		if(fabs(massBerror) < 1e-16) massBerror = 1e-16;
	} else { 
		massBerror = 0.0; 
	}

	QSbalance += (QSin-QSout)*dt;	

	if(iter%Niter == 0){
		printf("\n\tBedload mass error %.3e\t",massBerror);
		printf("\n\tBedload rate IN %.6lf OUT %.6lf  [m3/s]",QSin,QSout);			
	}
	#endif		



// BOUNDARY CONDITIONS
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Upstream $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	h_set_inlet_conditions( nCells,
		A, Q, h, q, u, 
		c, Fr,
		B, P, 
		M,
		QIN, HIN);
	
	Qin = Q[0];	

	#if BEDLOAD
	h_set_bedload_inlet_conditions( nCells,
		bedloadCap, Ggrass, BetaB,
		rhos, poros, ds, 
		shieldsC, nManSed, n1,
		Ab, Qs, 
		tauB, Gcell, qs, 
		h, u, 
		n, M, zb, zr,
		dt,
		QSIN, ZBIN);

	QSin = Qs[0];
	#endif	
	
		
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Downstream $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	h_set_outlet_conditions( nCells,
		A,  Q,  h,  q, u,
		c,  Fr,
		B,  P, 
		M,  zb,
		HOUT, ZSOUT);

	Qout = Q[nCells-1];
	
	#if BEDLOAD
	h_set_bedload_outlet_conditions( nCells,
		bedloadCap, Ggrass, BetaB,
		rhos, poros, ds, 
		shieldsC, nManSed, n1,
		Ab, Qs, 
		tauB, Gcell, qs, 
		h, u, 
		n, M, zb, zr);

	QSout = Qs[nCells-1];
	#endif



// GRAPHICAL SOLUTION (linux only)
	#if DISPLAY
	if(iter%Niter == 0){
		gp << "\n set grid\n";
		for(ic=0; ic<nCells; ic++) {
			xcurve = x[ic];
			ycurve1 = zb[ic];
			ycurve2 = h[ic]+zb[ic];
			#if BEDLOAD
			ycurve3 = 200.+B[ic]*qs[ic];
			#else
			ycurve3 = u[ic];
			#endif

			curve1.push_back({xcurve,ycurve1});
			curve2.push_back({xcurve,ycurve2});
			curve3.push_back({xcurve,ycurve3});
		}

		gp << "plot ";
		gp << " '-' with lines lt 1 lw 2 lc rgb 'black' title 'z', ";
		gp << " '-' with linespoints lt 1 pt 1 lc rgb 'blue' title 'h+z', ";	
		#if BEDLOAD
		gp << " '-' with linespoints lt 1 pt 2 lc rgb 'red' title 'qs' \n";
		#else
		gp << " '-' with linespoints lt 1 pt 2 lc rgb 'red' title 'u' \n";
		#endif
		gp.send1d(curve1);
		gp.send1d(curve2);
		gp.send1d(curve3);
		gp.flush();

		curve1.clear();
		curve2.clear();
		curve3.clear();
	}
	#endif




// UPDATE TIME
	iter++;
	t = t+dt;

	if(dt <= 1e-9) { break;}
	
	
// DATA OUTPUT
	if(t >= nout*Toutput){
		//Cell data
		sprintf(filename, "outputFiles/celldata%d.out",nout);
		fp = fopen(filename,"w");
		for(ic=0; ic<nCells; ic++){
			fprintf(fp,"%.3lf\t %.6lf\t %.6lf\t %.6lf",
				x[ic],
				zb[ic],
				h[ic],
				u[ic]);
			#if BEDLOAD
			fprintf(fp,"\t %.6lf",
				qs[ic]);
			#endif
			fprintf(fp," \n");
		}
		fclose(fp);

		// Discharge data
		qFile=fopen("outputFiles/discharge.out","a+");
		fprintf(qFile,"%.3lf\t %.6lf\t %.6lf\t %.6lf\t %.6lf",
			t,
			massWtn,
			Qin,Qout,Qbalance);
		#if BEDLOAD
		fprintf(qFile,"\t %.6lf\t %.6lf\t %.6lf\t %.6lf",
			massBtn,
			QSin,QSout,QSbalance);	
		#endif
		fprintf(qFile," \n");
		fclose(qFile);
		
		nout += 1;	
	}	


}
//////////////////////////////////     END TIME LOOP     //////////////////////////////////////////////////////////
CPUtime = clock() - CPUtime ;
double time = ((float)CPUtime)/CLOCKS_PER_SEC;



// DISPLAY FINAL INFORMATION	
printf("\n\n>> Final Time: %.3lf seg",t);
printf("\n\n>> Computation time %.3lf seg",time);
printf("\n\n>> Simulation completed!");
printf("\n\n ");

logFile=fopen("outputFiles/log.out","a+");
fprintf(logFile,"\n>> Final Time: %.3lf seg",t);
fprintf(logFile,"\n\n>> Computation time %.3lf seg",time);
fprintf(logFile,"\n\n>> Simulation completed!");
fprintf(logFile,"\n\n ");
fclose(logFile);


// FREE MEMORY
	#if BEDLOAD
	free( Ab );
	free( Qs );

	free( tauB );
	free( Gcell );
	free( qs );

	free( DQs );

	free( lB );
	#endif

	free( A );
	free( Q );
	free( h );
	free( q );
	free( u );

	free( c );
	free( Fr );
	
	free( B );
	free( P );
	
	free( n );	
	free( zb );

	free( DU1 );
	free( DU2 );

	free( lSW );	

	free( solidWall  );



	free( x );
	free( dx );
	free( dist );
	free( M );
	free( zr );



} // End of main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







