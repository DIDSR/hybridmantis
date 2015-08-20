///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        hybridMANTIS v1.0		     //
// 			     //                fastDETECT2 - C code                  //
//			     //		   (optical photons transport)		     //
//  			     //							     //
//                           //               used for Load Balancing                //
//			     //							     //
//			     //////////////////////////////////////////////////////////
//
// 
//
//
// ****Disclaimer****
//  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
//  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
//  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
//  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
//  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
//  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
//  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
//  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are
//  derived from it, and any modified versions bear some notice that they have been modified. 
//
//	Detailed comments are available in "hybridMANTIS_cuda_ver1_0.cu" and "hybridMANTIS_c_ver1_0.c" files.
//
//	Associated publication: Sharma Diksha, Badal Andreu and Badano Aldo, "hybridMANTIS: a CPU-GPU Monte Carlo method for modeling indirect x-ray detectors with
//				columnar scintillators". Physics in Medicine and Biology, 57(8), pp. 2357–2372 (2012)
//
//
//	File:   	hybridMANTIS_c_ver1_0_LB.c 			
//	Author: 	Diksha Sharma (US Food and Drug Administration)
//	Email: 		diksha.sharma@fda.hhs.gov			
//	Last updated:  	Apr 13, 2012
// 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
//
//      Header libraries
//
/////////////////////////////////////////

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/////////////////////////////////////////
//
//      Global variables
//
/////////////////////////////////////////

#define max_photon_per_EDE 900000	// maximum number of optical photons that can be generated per energy deposition event (EDE)

#ifndef USING_CUDA
	#define mybufsizeT 2304000	// CPU buffer size: # of events sent to the CPU
#endif

/////////////////////////////////////////
//
//      Include kernel program
//
/////////////////////////////////////////
#include "kernel_cuda_c_ver1_0_LB.cu"


////////////////////////////////////////////////////////////////////////////
//				MAIN PROGRAM			          //
////////////////////////////////////////////////////////////////////////////

#ifndef USING_CUDA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// cpuoptlb():   Performs optical transport for finding optimal load using load balancing in the CPU 
//	  	 Input arguments: cptime, cpubufsize
//
//		 cptime: 	time taken by GPU to call this routine
//		 cpubufsize:  	CPU buffer size defined by user in PENELOPE tally code. Only to be used for load balancing. This is different from 'mybufsize'.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void cpuoptlb_(double *cptime, int *cpubufsize)
	{    

		float dcos[3]={0}; 		// directional cosines
		float normal[3]={0}; 		// normal to surface in case of TIR
		float pos[3] = {0}; 		// new position
		float old_pos[3] = {0};  	// source coordinates
	
		// command line arguments
		float xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, lbound_x, lbound_y, ubound_x, ubound_y, d_max, yield, sensorRefl;
		int pixelsize, num_primary, min_optphotons, max_optphotons, num_bins;

		int nbytes = (*cpubufsize)*sizeof(struct start_info);
		struct start_info *structa;
		structa = (struct start_info*) malloc(nbytes);
		if( structa == NULL )
			printf("\n Struct start_info array CANNOT BE ALLOCATED - %d !!", (*cpubufsize));

		// get cpu time
		clock_t start, end;
		double num_sec;

		// get current time stamp to initialize seed input for RNG
		time_t seconds;
		seconds = time (NULL);
		struct timeval tv;
		gettimeofday(&tv,NULL);

		float rr=0.0f, theta=0.0f;
		float r=0.0f;		// random number
		float norm=0.0f;
		int jj=0;

		// initialize random number generator (RANECU)
		int seed_input = 271828182 ; // ranecu seed input
		int seed[2];

		// gsl variables
		const gsl_rng_type * Tgsl;
		gsl_rng * rgsl;
		double mu_gsl;	
		int my_index=0;
		int result_algo = 0;
		unsigned long long int *num_rebound;

		// output image variables
		int xdim = 0;
		int ydim = 0;
		int indexi=0, indexj=0;

	      	// create a generator chosen by the environment variable GSL_RNG_TYPE
	       	gsl_rng_env_setup();	     
	       	Tgsl = gsl_rng_default;
	       	rgsl = gsl_rng_alloc (Tgsl);

 		// copy to local variables from PENELOPE buffers
		xdetector = inputargs_.detx;		// x dimension of detector (in um). x in (0,xdetector)
		ydetector = inputargs_.dety;		// y dimension of detector (in um). y in (0,ydetector)
		height = inputargs_.detheight;		// height of column and thickness of detector (in um). z in range (-H/2, H/2)
		radius = inputargs_.detradius;		// radius of column (in um).
		n_C = inputargs_.detnC;			// refractive index of columns
		n_IC = inputargs_.detnIC;		// refractive index of intercolumnar material
		top_absfrac = inputargs_.dettop;	// column's top surface absorption fraction (0.0, 0.5, 0.98)
		bulk_abscoeff = inputargs_.detbulk;	// column's bulk absorption coefficient (in um^-1) (0.001, 0.1 cm^-1) 
		beta = inputargs_.detbeta;		// roughness coefficient of column walls
		d_min = inputargs_.detdmin;		// minimum distance a photon can travel when transmitted from a column
		d_max = inputargs_.detdmax;
		lbound_x = inputargs_.detlboundx;	// x lower bound of region of interest of output image (in um)
		lbound_y = inputargs_.detlboundy;	// y lower bound (in um)
		ubound_x = inputargs_.detuboundx;	// x upper bound (in um) 
		ubound_y = inputargs_.detuboundy;	// y upper bound (in um)
		yield = inputargs_.detyield;		// yield (/eV)
		pixelsize = inputargs_.detpixel;	// 1 pixel = pixelsize microns (in um)
		sensorRefl = inputargs_.detsensorRefl;	// Non-Ideal sensor reflectivity (%)
		num_primary = inputargs_.mynumhist;	// total number of primaries to be simulated
		min_optphotons = inputargs_.minphotons;	// minimum number of optical photons detected to be included in PHS
		max_optphotons = inputargs_.maxphotons;	// maximum number of optical photons detected to be included in PHS
		num_bins = inputargs_.mynumbins;	// number of bins for genrating PHS
		
		// dimensions of PRF image
		xdim = ceil((ubound_x - lbound_x)/pixelsize);
		ydim = ceil((ubound_y - lbound_y)/pixelsize);
		unsigned long long int myimage[xdim][ydim];
		
		// memory for storing histogram of # photons detected/primary
		int *h_num_detected_prim = 0;		
		h_num_detected_prim = (int*)malloc(sizeof(int)*num_primary);
			
		for(indexj=0; indexj < num_primary; indexj++)
		  h_num_detected_prim[indexj] = 0;
		  	 
		// start the clock
		start = clock();

	for(my_index = 0; my_index < (*cpubufsize); my_index++)		// iterate over x-rays
	{

		// reset the global counters
		num_generatedT = 0;
		num_detectT=0;
		num_abs_topT=0;	
		num_abs_bulkT=0;	
		num_lostT=0;
		num_outofcolT=0;
		num_theta1T=0;
		photon_distanceT=0.0f;

		// copying fortran buffer into *structa

		// units in the penelope output file are in cm. Convert to microns.
		structa[my_index].str_x = optical_.xbufopt[my_index] * 10000.0f;	// x-coordinate of interaction event.
		structa[my_index].str_y = optical_.ybufopt[my_index] * 10000.0f;	// y-coordinate
		structa[my_index].str_z = optical_.zbufopt[my_index] * 10000.0f;	// z-coordinate
		structa[my_index].str_E = optical_.debufopt[my_index];			// energy deposited
		structa[my_index].str_histnum = optical_.nbufopt[my_index];		// x-ray history number

		// sample # optical photons based on light yield and energy deposited for this interaction event (using Poisson distribution)
		mu_gsl = (double)structa[my_index].str_E * yield;
		structa[my_index].str_N = gsl_ran_poisson(rgsl,mu_gsl);

		if(structa[my_index].str_N > max_photon_per_EDE)
		{
			printf("\n\n CPU str_n exceeds max photons. program is exiting - %d !! \n\n",structa[my_index].str_N);
			exit(0);
		}

		num_rebound = (unsigned long long int*) malloc(structa[my_index].str_N*sizeof(unsigned long long int));
		if(num_rebound == NULL)
			printf("\n Error allocating num_rebound memory !\n");


		// initialize the RANECU generator in a position far away from the previous history:
		seed_input = (int)(seconds/3600+tv.tv_usec);			// seed input=seconds passed since 1970+current time in micro secs
		init_PRNG(my_index, 50000, seed_input, seed);      		// intialize RNG

		for(jj=0; jj<structa[my_index].str_N; jj++)
			num_rebound[jj] = 0;

		// reset the vectors
		dcos[0]=0.0f; dcos[1]=0.0f; dcos[2]=0.0f;
		normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;

		// set starting location of photon
		pos[0] = structa[my_index].str_x; pos[1] = structa[my_index].str_y; pos[2] = structa[my_index].str_z;	
		old_pos[0] = structa[my_index].str_x; old_pos[1] = structa[my_index].str_y; old_pos[2] = structa[my_index].str_z;

		// initializing the direction cosines for the first particle in each core
		r = (ranecu(seed) * 2.0f) - 1.0f; // random number between (-1,1)
		 	
		while(fabs(r) <= 0.01f)	
		 {
		   	r = (ranecu(seed) * 2.0f) - 1.0f;  	
		 }

		dcos[2] = r;		// random number between (-1,1)
		rr = sqrt(1.0f-r*r);
		theta=ranecu(seed)*twopipen;
		dcos[0]=rr*cos(theta);
		dcos[1]=rr*sin(theta);

		norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))	// normalize
		 {
			dcos[0] = dcos[0]/norm;
			dcos[1] = dcos[1]/norm;
			dcos[2] = dcos[2]/norm;
		 }

		local_counterT=0;
		while(local_counterT < structa[my_index].str_N)
		 { 
			
			absorbedT = 0;
			detectT = 0;
			bulk_absT = 0;

			// set starting location of photon
			pos[0] = structa[my_index].str_x; pos[1] = structa[my_index].str_y; pos[2] = structa[my_index].str_z;	
			old_pos[0] = structa[my_index].str_x; old_pos[1] = structa[my_index].str_y; old_pos[2] = structa[my_index].str_z;
			num_generatedT++;
			result_algo = 0;

			while(result_algo == 0)
			 {
			  	result_algo = algoT(normal, old_pos, pos, dcos, num_rebound, seed, structa[my_index], &myimage[0][0], xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, sensorRefl, d_max, ydim, h_num_detected_prim);      
			 }

		 }	

		// release resources
		free(num_rebound);

	}	// my_index loop ends


	// end the clock		
	end = clock();

	num_sec = (double)((end - start)/CLOCKS_PER_SEC);
        *cptime = num_sec;
        
	// release resources
	free(structa);
	free(h_num_detected_prim);

		return;
	}	// C main() ends
	
#endif
