///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        hybridMANTIS v1.0		     //
// 			     //              fastDETECT2 - CUDA code                 //
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
//	File:   	hybridMANTIS_cuda_ver1_0_LB.cu 			
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

#ifdef USING_CUDA
	#include <cutil_inline.h>
	#include <vector_types.h>
	#include <stdint.h>
#endif

/////////////////////////////////////////
//
//      Global variables
//
/////////////////////////////////////////

#define max_photon_per_EDE 900000	// maximum number of optical photons that can be generated per energy deposition event (EDE)

#ifdef USING_CUDA
	#define gpubufsizeT 2304000	// GPU buffer size: # of events sent to the GPU
#endif

/////////////////////////////////////////
//
//      Include kernel program
//
/////////////////////////////////////////
#include "kernel_cuda_c_ver1_0_LB.cu"

/////////////////////////////////////////
//
//      CUDA parameters
//
/////////////////////////////////////////
#ifdef USING_CUDA
	#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

	#define GRIDSIZE 18000		// number of blocks
	#define BLOCKSIZE 128		// number of threads
#endif




////////////////////////////////////////////////////////////////////////////
//				MAIN PROGRAM			          //
////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// gpuoptlb():   Performs optical transport for finding optimal load using load balancing in the GPU 
//	  	 Input arguments: gptime, lbbuf
//
//		 gptime: time taken by GPU to call this routine
//		 lbbuf:  GPU buffer size defined by user in PENELOPE tally code. Only to be used for load balancing. This is different from 'gpubufsize'.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	extern "C" void gpuoptlb_(double *gptime, int *lbbuf)
	{    

		float xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, lbound_x, lbound_y, ubound_x, ubound_y, d_max, yield, sensorRefl;
		int pixelsize, num_primary, min_optphotons, max_optphotons, num_bins;
		dim3 threads, blocks;
                float gpuTime = 0.0f;
		int devID;
		double mu_gsl;	
		unsigned long long int host_num_generated = 0, host_num_detect = 0;
		unsigned long long int host_num_lost = 0;
		unsigned long long int host_num_abs_top = 0, host_num_abs_bulk = 0;
		unsigned long long int host_num_outofcol = 0;
		unsigned long long int host_num_theta1 = 0;
		const gsl_rng_type * Tgsl;
		gsl_rng * rgsl;
		int xdim = 0;
		int ydim = 0;
		int indexi=0, indexj=0;
		int my_index=0;
		size_t pitch;
		int nbytes = (*lbbuf)*sizeof(struct start_info);

		// allocate memory pointers
		unsigned long long int *myimage = 0;	// device memory for output image
		unsigned long long int *h_myimage = 0; 	// host memory for output image     
		int *num_detected_primary = 0;		// device memory for # detected photons/primary
		int *h_num_detected_primary = 0;		// host memory to get # detected/primary
		struct start_info *h_a = 0;             // pointer to the struct info data in the host memory
		struct start_info *d_a = 0;             // pointers to struct data in the device memory


		// set the device with max GFlops	
		devID = cutGetMaxGflopsDeviceId();
		cudaSetDevice( devID );

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
		
	      	// create a generator chosen by the 
		//  environment variable GSL_RNG_TYPE 
	       	gsl_rng_env_setup();	     
	       	Tgsl = gsl_rng_default;
	       	rgsl = gsl_rng_alloc (Tgsl);
	       	
	       	// dimensions of PRF image
		xdim = ceil((ubound_x - lbound_x)/pixelsize);
		ydim = ceil((ubound_y - lbound_y)/pixelsize);

		// allocate device memory for storing output arrays
		cudaMallocPitch((void**)&myimage, &pitch, xdim*sizeof(unsigned long long int), ydim);		// allocate 2D image array
		cutilSafeCall( cudaMemset2D(myimage, pitch, 0, xdim*sizeof(unsigned long long int), ydim) );	// initialize to 0
		cutilSafeCall( cudaMalloc((void**)&num_detected_primary, sizeof(int)*num_primary) );		// outputting # detected/primary
		cutilSafeCall( cudaMemset(num_detected_primary, 0, sizeof(int)*num_primary) );			// initialize to 0

		// allocate host and device memory for stroing interaction event buffer information
		cutilSafeCall( cudaMallocHost((void**)&h_a, nbytes) ); 
		cutilSafeCall( cudaMalloc((void**)&d_a, nbytes) );

		// reset the host counters
		host_num_generated=0;
		host_num_detect=0;
		host_num_abs_top=0;	
		host_num_abs_bulk=0;	
		host_num_lost=0;
		host_num_outofcol=0;
		host_num_theta1=0;

		// reset device counters to zero
		cutilSafeCall(cudaMemcpyToSymbol("num_detectT",&host_num_detect,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
		cutilSafeCall(cudaMemcpyToSymbol("num_generatedT",&host_num_generated,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
		cutilSafeCall(cudaMemcpyToSymbol("num_abs_topT",&host_num_abs_top,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
		cutilSafeCall(cudaMemcpyToSymbol("num_abs_bulkT",&host_num_abs_bulk,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
		cutilSafeCall(cudaMemcpyToSymbol("num_lostT",&host_num_lost,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMemcpyToSymbol("num_outofcolT",&host_num_outofcol,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMemcpyToSymbol("num_theta1T",&host_num_theta1,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));


		// synchronize threads to ensure that previous kernel has finished
		cudaThreadSynchronize();
    
		cutilSafeCall( cudaMallocHost((void**)&h_myimage, xdim*ydim*sizeof(unsigned long long int)) ); 
		cutilSafeCall( cudaMallocHost((void**)&h_num_detected_primary, sizeof(int)*num_primary) );

		for(indexj=0; indexj < num_primary; indexj++)
		  h_num_detected_primary[indexj] = 0;
		  
		int *h_histogram = 0;		// host memory for storing histogram of # photons detected/primary
		h_histogram = (int*)malloc(sizeof(int)*num_bins);
			
		for(indexj=0; indexj < num_bins; indexj++)
		  h_histogram[indexj] = 0;


		// assign number of threads and blocks
		threads = dim3(BLOCKSIZE,1);
		blocks = dim3(GRIDSIZE,1);

		// reading data from lbbuf
		for(my_index = 0; my_index < (*lbbuf); my_index++)		// iterate over x-rays
		{

			// units in the penelope output file are in cm. Convert to microns.
			h_a[my_index].str_x = optical_.xbufopt[my_index] * 10000.0f;	// x-coordinate of interaction event.
			h_a[my_index].str_y = optical_.ybufopt[my_index] * 10000.0f;	// y-coordinate
			h_a[my_index].str_z = optical_.zbufopt[my_index] * 10000.0f;	// z-coordinate
			h_a[my_index].str_E = optical_.debufopt[my_index];		// energy deposited
			h_a[my_index].str_histnum = optical_.nbufopt[my_index];		// x-ray history number

			// sample # optical photons based on light yield and energy deposited for this interaction event (using Poisson distribution)
			mu_gsl = (double)h_a[my_index].str_E * yield;
			h_a[my_index].str_N = gsl_ran_poisson(rgsl,mu_gsl);

			if(h_a[my_index].str_N > max_photon_per_EDE)
			{
				printf("\n\n GPU str_n exceeds max photons. program is exiting - %d !! \n\n", h_a[my_index].str_N);
				exit(0);
			}

		} // for loop ends

		    // create cuda event handles
		    cudaEvent_t start, stop;
		    cutilSafeCall( cudaEventCreate(&start) );
		    cutilSafeCall( cudaEventCreate(&stop)  );
    

		    // execute the kernel 
			cutilSafeCall( cutilDeviceSynchronize() );
		
			cudaEventRecord(start, 0);

			// asynchronously copy data from host to device	(all to stream 0)
			cutilSafeCall( cudaMemcpyAsync(d_a, h_a, nbytes, cudaMemcpyHostToDevice, 0) );

			// each kernel has BLOCKSIZE threads; each thread transports one event in the buffer (info.str_N optical photons)
			algoT<<<blocks, threads, 0, 0>>>(d_a, myimage, num_detected_primary, pitch, (*lbbuf), xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, d_max, sensorRefl); 
				
			// asynchronously copy image data from device to host
			cutilSafeCall( cudaMemcpy2DAsync((void*)h_myimage,sizeof(unsigned long long int)*xdim,(void*)myimage,pitch,sizeof(unsigned long long int) *xdim,ydim,cudaMemcpyDeviceToHost, 0) );
			cutilSafeCall( cudaMemcpyAsync(h_num_detected_primary, num_detected_primary, sizeof(int)*num_primary, cudaMemcpyDeviceToHost, 0) );

			cudaEventRecord(stop, 0);

			 // have CPU do some work while waiting for stage 1 to finish
			 unsigned long int counter123=0;
			 while( cudaEventQuery(stop) == cudaErrorNotReady )
			 {
				counter123++;
			 }
			 cutilSafeCall( cudaEventElapsedTime(&gpuTime, start, stop) );
			 *gptime = (double)(gpuTime*0.001);	// convert in sec
                                 
		         cutilCheckMsg("algo() execution failed\n");

		         cudaThreadSynchronize();	// to ensure that gpu finished before copying back the final results.


		// copy counters from device to host
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_detect,num_detectT,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));	
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_generated,num_generatedT,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));	
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_abs_top,num_abs_topT,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));	
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_abs_bulk,num_abs_bulkT,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_lost,num_lostT,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_outofcol,num_outofcolT,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));
		cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_theta1,num_theta1T,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));


		// add h_myimage to the new_myimage
		for(indexi = 0; indexi < ydim; indexi++)
		 for(indexj = 0; indexj < xdim; indexj++)
			outputimage_.newimageopt[indexi][indexj] = outputimage_.newimageopt[indexi][indexj] + h_myimage[indexi*xdim + indexj];

		// make histogram of number of detected photons/primary for num_bins
		int binsize=0, newbin=0;
		int bincorr=0;
							
		binsize = floor((max_optphotons-min_optphotons)/num_bins);	// calculate size of each bin. Assuming equally spaced bins.
		bincorr = floor(min_optphotons/binsize);			// correction in bin number if min_optphotons > 0.
			
		for(indexi = 0; indexi < num_primary; indexi++)
		 {
		 	newbin = floor(h_num_detected_primary[indexi]/binsize) - bincorr;	// find bin #
		 	
		 	if(h_num_detected_primary[indexi] > 0)	// store only non-zero bins
		 	{
		 		if(h_num_detected_primary[indexi] <= min_optphotons)	// # detected < minimum photons given by user, add to the 1st bin
			 		h_histogram[0]++;
			 	else if(h_num_detected_primary[indexi] >= max_optphotons)	// # detected > maximum photons given by user, then add to the last bin
			 		h_histogram[num_bins-1]++;
			 	else
				 	h_histogram[newbin]++; 
			}
		 }
			
		// add num_detected_primary to gldetprimary array in PENELOPE
		for(indexi = 0; indexi < num_bins; indexi++)
			outputdetprim_.gldetprimary[indexi] = outputdetprim_.gldetprimary[indexi] + h_histogram[indexi];
				
			   
		// type cast unsigned long long int to double
		double cast_host_num_generated;
		double cast_host_num_detect;
		double cast_host_num_abs_top;
		double cast_host_num_abs_bulk;
		double cast_host_num_lost;
		double cast_host_num_outofcol;
		double cast_host_num_theta1;
		double cast_gputime;

		cast_host_num_generated = (double)host_num_generated;
		cast_host_num_detect    = (double)host_num_detect;
		cast_host_num_abs_top   = (double)host_num_abs_top;
		cast_host_num_abs_bulk  = (double)host_num_abs_bulk;
		cast_host_num_lost      = (double)host_num_lost;
		cast_host_num_outofcol  = (double)host_num_outofcol;
		cast_host_num_theta1    = (double)host_num_theta1;
		cast_gputime		= (double)(gpuTime);

		 // save to global counters
		 optstats_.glgen      = optstats_.glgen      + cast_host_num_generated;
		 optstats_.gldetect   = optstats_.gldetect   + cast_host_num_detect;
		 optstats_.glabstop   = optstats_.glabstop   + cast_host_num_abs_top;
		 optstats_.glabsbulk  = optstats_.glabsbulk  + cast_host_num_abs_bulk;
		 optstats_.gllost     = optstats_.gllost     + cast_host_num_lost;
		 optstats_.gloutofcol = optstats_.gloutofcol + cast_host_num_outofcol;
		 optstats_.gltheta1   = optstats_.gltheta1   + cast_host_num_theta1;
		 optstats_.glgputime  = optstats_.glgputime  + cast_gputime;

	 
		 // release resources
		 cutilSafeCall(cudaFree(d_a));
		 cutilSafeCall(cudaFree(myimage));
		 cutilSafeCall(cudaFree(num_detected_primary));
		 cudaFreeHost(h_a);
		 cudaFreeHost(h_myimage);
		 cudaFreeHost(h_num_detected_primary);
		 
    	         free(h_histogram);
		
		return;
	}	// CUDA main() ends
	
#endif
