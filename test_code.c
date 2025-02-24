//   Test code
//  ===========================================================================================
//   Version 1.0 - Jan 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza   
//  =========================================================================================== 


#include "define.h"

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
	qFile=fopen("outputFiles/discharge.out","w");
	
	#if DISPLAY
	Gnuplot gp;			// Graphical output pipe (linux only)
	#endif



// SIMULATOR HEADERS
	printf("\n   Test code");
	printf("\n   Version 1.0 - Jan 2025");
	printf("\n   Computational Hydraulics Group - University of Zaragoza");
	printf("\n   --------------------------------------------------");

	fprintf(logFile,"\n   Test code");
	fprintf(logFile,"\n   Version 1.0 - Jan 2025");
	fprintf(logFile,"\n   Computational Hydraulics Group - University of Zaragoza");
	fprintf(logFile,"\n   --------------------------------------------------");



// GEOMETRY DATA
	double x[nCells];
	double zb[nCells];

	for(ic=0; ic<nCells; ic++) {
		x[ic] = ic*dx;
		zb[ic] = ZRlevel;
	}

	int nWalls=nCells-1;


// MAIN VARIABLE DEFINITION

	clock_t CPUtime;

	int iter, nout;	
	double t, dt;
	
	double phi[nCells];	
	double flux;

	// Auxiliars
	double aux1,aux2,aux3;
	double xcurve, ycurve1, ycurve2;
	std::vector<std::vector<double>> curve1;
	std::vector<std::vector<double>> curve2; 



// VARIABLE INITIALIZATION
	iter = 1;					// Iteration number
	t = 0.0;					// Initial time
	nout = 1;

	for(ic=0; ic<nCells; ic++) {
		phi[ic]=0.0;
	}







/////////////////////////////////////////////////////////////////////////
// INITIAL CONDITIONS 
/////////////////////////////////////////////////////////////////////////

for(ic=0; ic<nCells; ic++) {
	phi[ic]=exp(-C1*(ic*dx-C2)*(ic*dx-C2));
}


CPUtime = clock();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////     TIME LOOP     /////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
while(t <= simTime) {	

	dt=0.5*dx/a;

	for(ic=0; ic<nWalls; ic++) {
		flux = a*phi[ic];
		phi[ic] -= flux*dt/dx;
		phi[ic+1] += flux*dt/dx;
	}

	printf("\n>> Time: %.3lf seg  -  Iteration %d  -  Time step %.6lf",t,iter,dt);	

// GRAPHICAL SOLUTION (linux only)
	#if DISPLAY
		gp << "\n set grid\n";
		for(ic=0; ic<nCells; ic++) {
			xcurve = x[ic];
			ycurve1 = zb[ic];
			ycurve2 = zb[ic]+phi[ic];

			curve1.push_back({xcurve,ycurve1});
			curve2.push_back({xcurve,ycurve2});
		}

		gp << "plot ";
		gp << " '-' with lines lt 1 lc rgb 'black' title 'z', ";	
		gp << " '-' with linespoints lt 1 pt 2 lc rgb 'red' title 'phi' \n";
		gp.send1d(curve1);
		gp.send1d(curve2);
		gp.flush();

		curve1.clear();
		curve2.clear();
	#endif




// UPDATE TIME
	iter++;
	t = t+dt;
	
// DATA OUTPUT
	if(t >= nout*Toutput){
		//Cell data
		sprintf(filename, "outputFiles/celldata%d.out",nout);
		fp = fopen(filename,"w");
		for(ic=0; ic<nCells; ic++){
			fprintf(fp,"%.3lf\t %.6lf\t %.6lf \n",
				x[ic],
				zb[ic],
				phi[ic]);
		}
		fclose(fp);

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

fprintf(logFile,"\n>> Final Time: %.3lf seg",t);
fprintf(logFile,"\n\n>> Computation time %.3lf seg",time);
fprintf(logFile,"\n\n>> Simulation completed!");
fprintf(logFile,"\n\n ");






// CLOSING DATA STORAGE
fclose(logFile);
fclose(qFile);



} // End of main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







