///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        hybridMANTIS v1.0		     //
// 			     //   	  fastDETECT2 kernel - CUDA + C  	     //
//			     //		   (optical photons transport)		     //
//  			     //							     //
//                           //               used for Load Balancing                //
//			     //							     //
//			     //////////////////////////////////////////////////////////
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
//	Detailed comments included in "kernel_cuda_c_ver1_0.cu" file.
//
//	Associated publication: Sharma Diksha, Badal Andreu and Badano Aldo, "hybridMANTIS: a CPU-GPU Monte Carlo method for modeling indirect x-ray detectors with
//				columnar scintillators". Physics in Medicine and Biology, 57(8), pp. 2357–2372 (2012)
//
//
//
//	File:   	kernel_cuda_c_ver1_0_LB.cu 			
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

	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <sys/time.h>
	#include <time.h>

/////////////////////////////////////////
//
//       Constants
//
/////////////////////////////////////////

	#define twopipen 6.283185308	// 2*PI
	#define pi 3.14159265		// PI
	#define epsilon 8.1929093e-6	// a very small number for float comparisons


/////////////////////////////////////////////////////////////////////////////////////
//
//     Data structure for storing a scintillation event location and deposited energy
//
/////////////////////////////////////////////////////////////////////////////////////

	struct start_info
	{
		double str_x;		
		double str_y;		
		double str_z;		
		double str_E;		
		int str_histnum;	
		int str_N;		
	};



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Fortran structure declarations - using PENELOPE 2006 (coded in Fortran)
//	A 'common' block in Fortran needs to be declared here to allow calling function interexchangebly.	
// 	
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Similar structure to 'start_info' - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			double xbufopt[gpubufsizeT];	
			double ybufopt[gpubufsizeT];	
			double zbufopt[gpubufsizeT];	
			double debufopt[gpubufsizeT];	
			int nbufopt[gpubufsizeT];	
			int myctropt;			
		        int cpu_num_real;		
		} optical_;
	#else
		extern struct
		{
			double xbufopt[mybufsizeT];
			double ybufopt[mybufsizeT];
			double zbufopt[mybufsizeT];
			double debufopt[mybufsizeT];	
			int nbufopt[mybufsizeT];		
			int myctropt;			
		        int cpu_num_real;			
		} optical_;
	#endif

// Storing optical output statistics - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			double glgen;			
			double gldetect;		
			double glabstop;		
			double glabsbulk;		
			double gllost;			
			double gloutofcol;		
			double gltheta1;		
			double glgputime;		
		} optstats_;
	#else
		extern struct
		{
			double glgen;
			double gldetect;
			double glabstop;
			double glabsbulk;
			double gllost;
			double gloutofcol;
			double gltheta1;
			double glgputime;		
		} optstats_;
	#endif

// Storing deposited energy and # optical photons detected - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			int gldetprimary[1000];		
		} outputdetprim_;
	#else
		extern struct
		{
			int gldetprimary[1000];
		} outputdetprim_;
	#endif

// Structure for storing point response functions - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			unsigned long long int newimageopt[501][501];	
			unsigned long long int tempimageopt[501][501];
		} outputimage_;
	#else
		extern struct
		{
			unsigned long long int newimageopt[501][501];
			unsigned long long int tempimageopt[501][501];
		} outputimage_;
	#endif

// Storing the memory addresses of arrays, in order to call fastDETECT2 in the GPU asynchronously - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			unsigned long long int gpuimage;		
			unsigned long long int gpudetect;		
			unsigned long long int hosta;			
			unsigned long long int deva;			
			unsigned long long int devpitch;			
		} gpumemaddr_;
	#endif

// Storing the input arguments - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			double detx;		
			double dety;		
			double detheight;	
			double detradius;	
			double detnC;		
			double detnIC;		
			double dettop;		
			double detbulk;		
			double detbeta;		
			double detdmin;		
			double detdmax;		
			double detlboundx;	
			double detlboundy;	
			double detuboundx;	
			double detuboundy;	
			double detyield;	
			double detsensorRefl;	
			int detpixel;		
			int rungpu;
			int machinenum;
			int mynumhist;	
			int minphotons;	
			int maxphotons;	
			int mynumbins;	
		} inputargs_;
	#else
		extern struct
		{
			double detx;
			double dety;
			double detheight;
			double detradius;
			double detnC;
			double detnIC;
			double dettop;
			double detbulk;
			double detbeta;
			double detdmin;
			double detdmax;
			double detlboundx;
			double detlboundy;
			double detuboundx;
			double detuboundy;
			double detyield;
			double detsensorRefl;
			int detpixel;
			int rungpu;
			int machinenum;
			int mynumhist;	
			int minphotons;	
			int maxphotons;	
			int mynumbins;	
		} inputargs_;
	#endif

/////////////////////////////////////////
//
//       Function declarations
//
/////////////////////////////////////////

// transports optical photon from its generation until it ends (detected/absorbed/lost).
	#ifdef USING_CUDA
		__global__ void algoT(struct start_info *info, unsigned long long int *myimage, int *num_detected_primary, size_t pitch, int rowsread, float xdetector, 
		float ydetector, float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y,
		float ubound_x, float ubound_y, float d_max, float sensorRefl);	
	#else
		int algoT(float *normal, float *old_pos, float *pos, float *dcos, unsigned long long int *num_rebound, int* seed, struct start_info info, 
		unsigned long long int *myimage, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, 
		float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, float sensorRefl, float d_max, int ydim, int *h_num_detected_prim);
	#endif

// photon within a column. calculate if it gets absorbed or moves inside the column.
	#ifdef USING_CUDA
		__device__ inline int isotropicT(float3 *pos, float3 *dcos, int2* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector,
	 	struct start_info *info, unsigned long long int mynum_rebound, float *XcT, float *YcT, int mytid);
	#else
		int isotropicT(float *pos, float *dcos, int* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector, 
		struct start_info info, unsigned long long int mynum_rebound);
	#endif

// photon within a column. calculate distance to next position in the same column and move it.
	#ifdef USING_CUDA
		__device__ float dist_to_surfaceT(float3 *pos, float3 *dcos, float R, float H, float xdetector, float ydetector, struct start_info *info, 
		unsigned long long int mynum_rebound, float *XcT, float *YcT, int mytid);	
	#else
		float dist_to_surfaceT(float *pos, float *dcos, float R, float H, float xdetector, float ydetector, struct start_info info, 
		unsigned long long int mynum_rebound);
	#endif

// photon within/between columns. calculate if it gets reflected or transmitted.
	#ifdef USING_CUDA
		__device__ int boundary_analysisT(float3 *normal, float3 *pos, float3 *dcos, int2* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, 
		float top_absfrac, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, float *XcT,
		float *YcT, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float d_max, float sensorRefl);	
	#else
		int boundary_analysisT(float *normal, float *pos, float *dcos, int* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, 			float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, 
		float d_max, float sensorRefl, int ydim, int *h_num_detected_prim);	
	#endif

// transmit photon to another column. calculates the new position when it transmits, build new column and move photon here.
	#ifdef USING_CUDA
		__device__ int transmitT(float3 *pos, float3 *dcos, float3 *normal, int2* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min,
		int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, size_t pitch, struct start_info *info, int mytid, 
		int *num_detected_primary, float d_max, float sensorRefl, int flagCCT);	
	#else
		int transmitT(float *pos, float *dcos, float *normal, int* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min, int pixelsize, 			float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, float d_max, float sensorRefl, int ydim, 
		int flagCCT, int *h_num_detected_prim);	
	#endif

// called when photon reflects from sensor plane (bottom surface) of the detector, outside of any column.
	#ifdef USING_CUDA
		__device__ int refl_bottomT(float3 *pos, float3 *dcos, float3 *normal, float xdetector, float ydetector, int2* seed, float beta, float d_min, float H, float d_max);	
	#else
		int refl_bottomT(float *pos, float *dcos, float *normal, float xdetector, float ydetector, int* seed, float beta, float d_min, float H, float d_max);	
	#endif

// calculate dot product of two vectors to give cosine of angle between them.
	#ifdef USING_CUDA
		__device__ inline float dot_productT(float3 *aa, float3 *b);	
	#else
		float dot_productT(float *aa, float *b);	
	#endif

// determine if photon got detected at sensor plane or is reflected back within the column
	#ifdef USING_CUDA
		__device__ inline int detectionT(float3 *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage,
		size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float sensorRefl, float d_min, int2* seed, float3 *dcos, float3 *normal, 
		float bulk_abscoeff, float R, float xdetector, float ydetector, unsigned long long int mynum_rebound, float *XcT, float *YcT);  
	#else
		int detectionT(float *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, 
		struct start_info info, float sensorRefl, float d_min, int* seed, float *dcos, float *normal, float bulk_abscoeff, float R, float xdetector, float ydetector, 
		unsigned long long int mynum_rebound, int ydim, int *h_num_detected_prim); 
	#endif

// calculate directional cosines of reflected/refracted vector.
	#ifdef USING_CUDA
		__device__ inline void trans_dir_cosT(float3 *dcos, float3 *normal, float refl_theta, float trans_theta, int flag_ref, int mytid, 
		struct start_info *info);  
	#else
		void trans_dir_cosT(float *dcos, float *normal, float refl_theta, float trans_theta, int flag_ref, struct start_info info);
	#endif

// calculate new rough normal vector depending on value of 'beta' (roughness coefficient).
	#ifdef USING_CUDA
		__device__ inline void RoughSurfaceT(float3 *normal, int2* seed, float beta);  
	#else
		void RoughSurfaceT(float *normal, int* seed, float beta); 
	#endif



// RANECU pseudo random number generator
	#ifdef USING_CUDA
		__device__ inline void init_PRNGT(int history_batch, int histories_per_thread, int seed_input, int2* seed);
		__device__ inline int abMODmT(int m, int a, int s);
		__device__ inline float ranecuT(int2* seed);
	#else
		void init_PRNGT(int history_batch, int histories_per_thread, int seed_input, int* seed);
		int abMODmT(int m, int a, int s);
		float ranecuT(int* seed);
	#endif


/////////////////////////////////////////
//
//       Global variables
//
/////////////////////////////////////////

	#ifdef USING_CUDA
		// counters
		__device__ unsigned long long int num_generatedT; 	
		__device__ unsigned long long int num_detectT;	 	
		__device__ unsigned long long int num_abs_topT;	 	
		__device__ unsigned long long int num_abs_bulkT;	
		__device__ unsigned long long int num_lostT;	 	
		__device__ unsigned long long int num_outofcolT;	
		__device__ unsigned long long int num_theta1T;	 
		__device__ float photon_distanceT;     		 
	#else
		// counters
		unsigned long long int num_generatedT=0;	
		unsigned long long int num_detectT=0;	
		unsigned long long int num_abs_topT=0;	
		unsigned long long int num_abs_bulkT=0;	
		unsigned long long int num_lostT=0;	
		unsigned long long int num_outofcolT=0;	 
		unsigned long long int local_counterT=0;	 
		unsigned long long int num_theta1T=0;	

		//flags
		int absorbedT=0;	
		int detectT=0;		
		int bulk_absT=0;	

		float XcT=0.0f;		
		float YcT=0.0f;
		float photon_distanceT=0.0f; 

		FILE *fp1;
	#endif

/////////////////////////////////////////
//
//    Functions definition
//
/////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// transports optical photon from its generation until it ends (detected/absorbed/lost).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__global__ void algoT(struct start_info *info, unsigned long long int *myimage, int *num_detected_primary, size_t pitch, int rowsread, float xdetector, float ydetector, 
	float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, 
	float ubound_y, float d_max, float sensorRefl)
	{

		unsigned long long int local_counter = 0; // total number of photons terminated (either detected at bottom, absorbed at top or in the bulk)
		unsigned long long int local_num_generated = 0;
		unsigned long long int num_rebound=0;
		float3 dcos, normal, pos; 
		float rr=0.0f, theta=0.0f;
		float XcT=0.0f;				// center x,y of the current cylinder
		float YcT=0.0f;
		int tid = threadIdx.x + blockIdx.x * blockDim.x;	// thread Id

		// flags
		int absorbedT=0;			// flag for photons absorbed at the top surface of the detector
		int detectT=0;				// flag for photons detected at bottom of the detector
		int bulk_absT=0;			// flag for photons absorbed in the material of a column


	if(tid < rowsread)	// number of threads = number of rows read from pen output file, should simulate photons, rest do nothing
	{
		int NUM_EACH_THREAD = info[tid].str_N;			// number of photons to be simulated by each thread

		// Initialize variables
		dcos.x = 0.0f; dcos.y = 0.0f; dcos.z = 0.0f;
		normal.x = 0.0f; normal.y = 0.0f; normal.z = 0.0f;

		pos.x = info[tid].str_x; pos.y = info[tid].str_y; pos.z = info[tid].str_z;	// starting location given by the host

		int seed_input = 271828182+tid; // ranecu seed input
		int2 seed;

		// Initialize the RANECU generator in a position far away from the previous history:
		init_PRNGT(tid, 50000, seed_input, &seed);     

		// Initalize the device memory - dcos
		dcos.z = (ranecuT(&seed) * 2.0f) - 1.0f;	// generate random number between -1.0 and 1.0

		rr = sqrt(1.0f - dcos.z*dcos.z);
		theta = ranecuT(&seed) * twopipen;	// generate random number between 0 and 2pi
	
		dcos.x = rr*cos(theta);
		dcos.y = rr*sin(theta);

	
		if (((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) < (1.0f - epsilon)) || ((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) > (1.0f + epsilon)))
		 {
			dcos.x = dcos.x/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
			dcos.y = dcos.y/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
			dcos.z = dcos.z/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
		 }

	   local_num_generated++;		// first particle generated for this thread

	   while(local_num_generated < (NUM_EACH_THREAD+1))	//run until NUM_EACH_THREAD particles generated
	     {

		if(absorbedT == 0)        
		 {
			bulk_absT = isotropicT(&pos, &dcos, &seed, bulk_abscoeff, R, H, xdetector, ydetector, &info[tid], num_rebound, &XcT, &YcT, tid);

			if(bulk_absT == 0)
			{
				detectT = detectionT(&pos, H, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, tid, num_detected_primary, sensorRefl, d_min, &seed, &dcos, &normal, bulk_abscoeff, R, xdetector, ydetector, num_rebound, &XcT, &YcT);
			}
		 }

	 
		if( ((detectT == 1) || (absorbedT == 1) || (bulk_absT == 1)) && (local_counter < (NUM_EACH_THREAD-1)) ) // particle terminated
		 {

			local_counter++;

			// re-initialize all the arrays
			dcos.z = (ranecuT(&seed) * 2.0f) - 1.0f;	// generate random number between -1.0 and 1.0

			rr = sqrt(1.0f - dcos.z*dcos.z);
			theta = ranecuT(&seed) * twopipen;	// generate random number between 0 and 2pi
	
			dcos.x = rr*cos(theta);
			dcos.y = rr*sin(theta);

			if (((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) < (1.0f - epsilon)) || ((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) > (1.0f + epsilon)))
			 {
				dcos.x = dcos.x/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
				dcos.y = dcos.y/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
				dcos.z = dcos.z/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
			 }


			// set starting location of photon
			pos.x = info[tid].str_x; pos.y = info[tid].str_y; pos.z = info[tid].str_z;
			normal.x = 0.0f; normal.y = 0.0f; normal.z = 0.0f;
	
			if(beta > 0.0f)
				RoughSurfaceT(&normal, &seed, beta);	// new normal for rough surface

			local_num_generated++;
			absorbedT = 0;
			detectT = 0;
			bulk_absT = 0;
			num_rebound = 0;
			XcT = 0.0f;
			YcT = 0.0f;

		 }
		else if( ((detectT == 1) || (absorbedT == 1) || (bulk_absT == 1)) && (local_counter == (NUM_EACH_THREAD-1)) )
		 {
			local_counter++;
			break;
	
		 }
		else if( (detectT == 0) && (absorbedT == 0) && (bulk_absT == 0) && (fabs(dcos.z - 0.0f) < epsilon) )  // checking for trapped particle going back and forth with dcos(z)=0
		 {
			// kill the particle and generate a new one instead - do not increment the counter

			// re-initialize all the arrays
	 		dcos.z = (ranecuT(&seed) * 2.0f) - 1.0f;	// generate random number between -1.0 and 1.0

			rr = sqrt(1.0f - dcos.z*dcos.z);
			theta = ranecuT(&seed) * twopipen;	// generate random number between 0 and 2pi
	
			dcos.x = rr*cos(theta);
			dcos.y = rr*sin(theta);

			if (((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) < (1.0f - epsilon)) || ((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) > (1.0f + epsilon)))
			 {
				dcos.x = dcos.x/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
				dcos.y = dcos.y/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
				dcos.z = dcos.z/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
			 }

			// set starting location of photon
			pos.x = info[tid].str_x; pos.y = info[tid].str_y; pos.z = info[tid].str_z;

			normal.x = 0.0f; normal.y = 0.0f; normal.z = 0.0f;
	
			if(beta > 0.0f)
				RoughSurfaceT(&normal, &seed, beta);	// new normal for rough surface

			absorbedT = 0;
			detectT = 0;
			bulk_absT = 0;
			num_rebound = 0;
			XcT = 0.0f;
			YcT = 0.0f;
		 }
		else
		 {
			num_rebound++;
		    	absorbedT = boundary_analysisT(&normal, &pos, &dcos, &seed, xdetector, ydetector, R, H, n1, n2, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, &XcT, &YcT, pitch, info, tid, num_detected_primary, d_max, sensorRefl);

		 }

	  } // while loop ends

	atomicAdd(&num_generatedT, local_num_generated);

	}	// if tid=rowsread ends

	 return;
	}
#else
	int algoT(float *normal, float *old_pos, float *pos, float *dcos, unsigned long long int *num_rebound, int* seed, struct start_info info,
	unsigned long long int *myimage, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float d_min, 
	int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, float sensorRefl, float d_max, int ydim, int *h_num_detected_prim)
	{


		float rr=0.0f, theta=0.0f;
		float norm = 0.0f;
		float rnd_num = 0.0f;
		int myresult = 0;


		if(absorbedT == 0)        
		 {
			bulk_absT = isotropicT(pos, dcos, seed, bulk_abscoeff, R, H, xdetector, ydetector, info, num_rebound[local_counterT]);

			if(bulk_absT == 0)
			{
				detectT = detectionT(pos, H, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, sensorRefl, d_min, seed, dcos, normal, bulk_abscoeff, R, xdetector, ydetector, num_rebound[local_counterT], ydim, h_num_detected_prim);
			}
		 }

	 
		if( (detectT == 1) || (absorbedT == 1) || (bulk_absT == 1) )
		 {
			local_counterT++;

			// re-initialize all the arrays
			rnd_num = (ranecuT(seed) * 2.0f) - 1.0f; 
	 	
			while(fabs(rnd_num) <= 0.01)	
			 {
		   		rnd_num = (ranecuT(seed) * 2.0f) - 1.0f;  	
		 	 }

			dcos[2] = rnd_num;		// random number between (-1,1)
			rr = sqrt(1.0-rnd_num*rnd_num);
			theta=ranecuT(seed)*twopipen;
			dcos[0]=rr*cos(theta);
			dcos[1]=rr*sin(theta);

			norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

			if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
			 {
				dcos[0] = dcos[0]/norm;
				dcos[1] = dcos[1]/norm;
				dcos[2] = dcos[2]/norm;
			 }


			// set starting location of photon
			pos[0] = info.str_x; pos[1] = info.str_y; pos[2] = info.str_z;
			old_pos[0] = info.str_x; old_pos[1] = info.str_y; old_pos[2] = info.str_z;

			normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 0.0f;
	
			if(beta > 0.0f)
			  RoughSurfaceT(normal, seed, beta);	// new normal for rough surface

			absorbedT = 0;
			detectT = 0;
			bulk_absT = 0;

			myresult = 1;
	
		 }
		else if( (detectT == 0) && (absorbedT == 0) && (bulk_absT == 0) && (fabs(dcos[2] - 0.0f) < epsilon) )  // checking for trapped particle going back and forth with dcos(z)=0
		 {
			// kill the particle and generate a new one instead - do not increment the counter

			// re-initialize all the arrays
			rnd_num = (ranecuT(seed) * 2.0f) - 1.0f; 
	 	
			while(fabs(rnd_num) <= 0.01)	
			 {
		   		rnd_num = (ranecuT(seed) * 2.0f) - 1.0f;  	
		 	 }

			dcos[2] = rnd_num;		// random number between (-1,1)
			rr = sqrt(1.0-rnd_num*rnd_num);
			theta=ranecuT(seed)*twopipen;
			dcos[0]=rr*cos(theta);
			dcos[1]=rr*sin(theta);

			norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

			if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
			 {
				dcos[0] = dcos[0]/norm;
				dcos[1] = dcos[1]/norm;
				dcos[2] = dcos[2]/norm;
			 }

			// set starting location of photon
			pos[0] = info.str_x; pos[1] = info.str_y; pos[2] = info.str_z;

			normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 0.0f;
	
			if(beta > 0.0f)
				RoughSurfaceT(normal, seed, beta);	// new normal for rough surface

			absorbedT = 0;
			detectT = 0;
			bulk_absT = 0;

			myresult = 0;
		 }
		else
		 {
			num_rebound[local_counterT]++;
		    	absorbedT = boundary_analysisT(normal, pos, dcos, seed, xdetector, ydetector, R, H, n1, n2, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, h_num_detected_prim);

	 		myresult = 0;

		 }

	 return myresult;

	}
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// determine where photon hits next within the column or if it gets absorbed in the material
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline int isotropicT(float3 *pos, float3 *dcos, int2* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector,
	struct start_info *info, unsigned long long int mynum_rebound, float *XcT, float *YcT, int mytid)
	{
		float dsurf = 999.0f;
		float dabs = 999.0f;
		int flag_bulkabs = 0;

		dsurf = dist_to_surfaceT(pos, dcos, R, H, xdetector, ydetector, info, mynum_rebound, XcT, YcT, mytid);	// distance to surface

		if (bulk_abscoeff > 0.0f)	
			dabs = (-1.0f/bulk_abscoeff) * log(ranecuT(seed));					// distance to absorption
		else
			dabs = 999999.0f;

		if (fabs(dsurf-(-99.0f)) < epsilon)				// particle lost because it went out of limit in dist_to_surface()
		{
			flag_bulkabs = 1;
		}
		else if ( (dsurf < dabs) && (dsurf >= 0.0f) )
		 {		
			flag_bulkabs = 0;
		 }
		else if ( (dsurf >= dabs) && (dabs >= 0.0f) )
		 {
			flag_bulkabs = 1;

			atomicAdd(&num_abs_bulkT,1);
		 }

	   return flag_bulkabs;
	}
#else
	int isotropicT(float *pos, float *dcos, int* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector,struct start_info info,
	unsigned long long int mynum_rebound)
	{
		float dsurf = 999.0f;
		float dabs = 999.0f;
		int flag_bulkabs = 0;

		dsurf = dist_to_surfaceT(pos, dcos, R, H, xdetector, ydetector, info, mynum_rebound);	// distance to surface

		if (bulk_abscoeff > 0.0f)	
			dabs = (-1.0f/bulk_abscoeff) * log(ranecuT(seed));					// distance to absorption
		else
			dabs = 999999.0f;

		if (fabs(dsurf-(-99.0f)) < epsilon)				// particle lost because it went out of limit in dist_to_surface()
		{
			flag_bulkabs = 1;
		}
		else if ( (dsurf < dabs) && (dsurf >= 0.0f) )
		 {		
			flag_bulkabs = 0;
		 }
		else if ( (dsurf >= dabs) && (dsurf >= 0.0f) )
		 {
			flag_bulkabs = 1;

			num_abs_bulkT++;
		 }

	   return flag_bulkabs;
	}
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// calculate distance to surface (within same column)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ float dist_to_surfaceT(float3 *pos, float3 *dcos, float R, float H, float xdetector, float ydetector, struct start_info *info, 
	unsigned long long int mynum_rebound, float *XcT, float *YcT, int mytid)
	{
		float d=999.0f;
		float d1=999.0f, d2=999.0f;
		float d_plane=999.0f, d_cyl=999.0f;
		float3 temp_pos = {0.0f};
		float3 my1 = {0.0f};
		float R1 = 999.0f, R2 = 999.0f;
		float stepsize = 0.5f;
		int repeat = 0, ctr1 = 0; 	// number of times photon should be moved in steps towards the column before killing. 
						// Valid only when goes out of column.

		temp_pos.x = pos->x;
		temp_pos.y = pos->y;
		temp_pos.z = pos->z;

		// center of first column (assumed as x,y position of the energy deposition from Penelope)
		if(mynum_rebound == 0)				
		{
			*XcT = info->str_x;
			*YcT = info->str_y;
		}

		// solving quadratic equation for distance from a point to the surface of cylinder
		my1.x = dcos->x*dcos->x + dcos->y*dcos->y;
		my1.y = 2.0f*( ((pos->x-(*XcT))*dcos->x) + ((pos->y-(*YcT))*dcos->y) );
		my1.z = (pos->x-(*XcT))*(pos->x-(*XcT)) + (pos->y-(*YcT))*(pos->y-(*YcT)) - (R*R);
	
		// actual distance d = (d1 or d2)/sin_theta2	
		d1 = (-my1.y + (sqrt( (my1.y*my1.y) - (4.0f*my1.x*my1.z) )))/(2.0f * my1.x);
		d2 = (-my1.y - (sqrt( (my1.y*my1.y) - (4.0f*my1.x*my1.z) )))/(2.0f * my1.x);


		// might hit either upper half surface or top
		if(dcos->z > 0.0f) 
		 {

		  if((fabs(dcos->x - 0.0f) < epsilon) && (fabs(dcos->y - 0.0f) < epsilon) && (fabs(dcos->z - 1.0f) < epsilon))  
		  // if particle travel straight in +z axis direction
		   {
			d = (H/2.0f - pos->z)/dcos->z;

			pos->z = H/2.0f;	

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;
		   }
		  else
		   {		

			// calculate distance to infinite plane at z=H/2
			d_plane = (H/2.0f - pos->z)/(dcos->z);

			// calculate the distance to the upper half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos->z = temp_pos.z + d*dcos->z;
			 }
			else
			 {
				d = d_plane;		
				pos->z = H/2.0f;
			 }

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;

		   }	// else loop ends
	
		
		 } // if loop for dcos.z > 0 ends

		else if(dcos->z < 0.0f) // might hit either lower half or bottom of cylinder
		 {
		  // if particle travels in -Z direction staright, then it should get detected
		  if ((fabs(dcos->x-0.0f) < epsilon) && (fabs(dcos->y-0.0f) < epsilon) && (fabs(dcos->z - (-1.0f)) < epsilon))  
		   {
			d = (-H/2.0f - pos->z)/dcos->z;  
			pos->z = -H/2.0f;
				
			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;
		   }
		  else
		   {

			// calculate distance to infinite plane at z=-H/2
			d_plane = (-H/2.0f - pos->z)/(dcos->z);

			// calculate the distance to the lower half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos->z = temp_pos.z + d*dcos->z;
			 }
			else
			 {
				d = d_plane;		
				pos->z = -H/2.0f;
			 }

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;

		   }	// else loop ends

		 }	// else if loop for dcos.z < 0 ends

		else	// when dcos.z=0.0 (will hit only the side of the cylinder)
		 {
			// calculate the distance to the side of cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		
		
			d = d_cyl;
			pos->z = temp_pos.z + d*dcos->z;

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;

		 }

		// condition to check that pos is within detector boundaries - if true, particle LOST
		if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) || (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
			{
				d = -99.0f;
				atomicAdd(&num_lostT,1);
				goto distexit;
			}
		else
			atomicAdd(&photon_distanceT, d);		// add distance travelled to global variable


		// CHECK IF PHOTON IS OUTSIDE THE CURRENT COLUMN
		R1 = sqrt((pos->x - (*XcT))*(pos->x - (*XcT)) + (pos->y - (*YcT))*(pos->y - (*YcT)));

		// check if photon is out of current column
		repeat = 0;
		ctr1 = 0;

		while( (R1 > (R-1e-5)) && (repeat < 10) && (ctr1 < 10) ) // R1 > R1-some small value..because of single precision errors that comparison with R may generate
		{

			// store current position
			temp_pos.x = pos->x;
			temp_pos.y = pos->y;

			// move particle by 0.5 um in the incident direction
			pos->x = pos->x + stepsize*(-dcos->x);
			pos->y = pos->y + stepsize*(-dcos->y);

			R2 = sqrt((pos->x - (*XcT))*(pos->x - (*XcT)) + (pos->y - (*YcT))*(pos->y - (*YcT)));

			if(R2 > R1) // means the photon is moving farther away from the column 
			{	    // this can happen if the stepsize if too big and the photon passes through the column and gets out on other side.

				// move it back to previous position, reduce the stepsize and try moving it again.
				pos->x = temp_pos.x;
				pos->y = temp_pos.y;

				stepsize = stepsize/2.0f;
				ctr1++;
			}
			else
			{
				R1 = R2;
				repeat++;
			}

		}

		// kill the particle if still outside the column
		if(R1 > (R-1e-5))
		 {
			d = -99.0f;
			atomicAdd(&num_outofcolT,1);
			goto distexit;
		 }

		// condition to check that pos is within detector boundaries - if true, particle LOST
		if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) || (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
			{
				d = -99.0f;
				atomicAdd(&num_lostT,1);
				goto distexit;
			}

	distexit:
	 return d;
	}	// CUDA dist_to_surface function ends
#else

	float dist_to_surfaceT(float *pos, float *dcos, float R, float H, float xdetector, float ydetector, struct start_info info, unsigned long long int mynum_rebound)
	{

		float d=999.0f;
		float d1=999.0f, d2=999.0f;
		float d_plane=999.0f, d_cyl=999.0f;
		float temp_pos[3] = {0.0f};
		float my1[3] = {0.0f};
		float R1 = 999.0f, R2 = 999.0f;
		float stepsize = 0.5f;
		int repeat = 0, ctr1 = 0; 	// number of times photon should be moved in steps towards the column before killing. 
					// Valid only when goes out of column.

		temp_pos[0] = pos[0];
		temp_pos[1] = pos[1];
		temp_pos[2] = pos[2];

		// center of first column (assumed as x,y position of the energy deposition from Penelope)
		if(mynum_rebound == 0)				
		{
			XcT = info.str_x;
			YcT = info.str_y;
		}

		// solving quadratic equation for distance from a point to the surface of cylinder
		my1[0] = dcos[0]*dcos[0] + dcos[1]*dcos[1];
		my1[1] = 2.0f*( ((pos[0]-(XcT))*dcos[0]) + ((pos[1]-(YcT))*dcos[1]) );
		my1[2] = (pos[0]-(XcT))*(pos[0]-(XcT)) + (pos[1]-(YcT))*(pos[1]-(YcT)) - (R*R);
	
		// actual distance d = (d1 or d2)/sin_theta2	
		d1 = (-my1[1] + (sqrt( (my1[1]*my1[1]) - (4.0f*my1[0]*my1[2]) )))/(2.0f * my1[0]);
		d2 = (-my1[1] - (sqrt( (my1[1]*my1[1]) - (4.0f*my1[0]*my1[2]) )))/(2.0f * my1[0]);

	
		// might hit either upper half surface or top
		if(dcos[2] > 0.0f) 
		 {

		  if((fabs(dcos[0] - 0.0f) < epsilon) && (fabs(dcos[1] - 0.0f) < epsilon) && (fabs(dcos[2] - 1.0f) < epsilon))  
		  // if particle travel straight in +z axis direction
		   {
			d = (H/2.0f - pos[2])/dcos[2];

			pos[2] = H/2.0f;	

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];
		   }
		  else
		   {		

			// calculate distance to infinite plane at z=H/2
			d_plane = (H/2.0f - pos[2])/(dcos[2]);

			// calculate the distance to the upper half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos[2] = temp_pos[2] + d*dcos[2];
			 }
			else
			 {
				d = d_plane;		
				pos[2] = H/2.0f;
			 }

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		   }	// else loop ends
	
		
		 } // if loop for dcos[2] > 0 ends

		else if(dcos[2] < 0.0f) // might hit either lower half or bottom of cylinder
		 {
		  // if particle travels in -Z direction staright, then it should get detected
		  if ((fabs(dcos[0]-0.0f) < epsilon) && (fabs(dcos[1]-0.0f) < epsilon) && (fabs(dcos[2] - (-1.0f)) < epsilon))  
		   {
			d = (-H/2.0f - pos[2])/dcos[2];  
			pos[2] = -H/2.0f;
				
			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];
		   }
		  else
		   {

			// calculate distance to infinite plane at z=-H/2
			d_plane = (-H/2.0f - pos[2])/(dcos[2]);

			// calculate the distance to the lower half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos[2] = temp_pos[2] + d*dcos[2];
			 }
			else
			 {
				d = d_plane;		
				pos[2] = -H/2.0f;
			 }

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		   }	// else loop ends

		 }	// else if loop for dcos[2] < 0 ends

		else	// when dcos[2]=0.0 (will hit only the side of the cylinder)
		 {
			// calculate the distance to the side of cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		
		
			d = d_cyl;
			pos[2] = temp_pos[2] + d*dcos[2];

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		 }

		// condition to check that pos is within detector boundaries - if true, particle LOST
		if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) || (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
			{
				d = -99.0f;
				num_lostT++;
				goto distexit;
			}
		else
			photon_distanceT = photon_distanceT + d;		// add distance travelled to global variable


		// CHECK IF PHOTON IS OUTSIDE THE CURRENT COLUMN
		R1 = sqrt((pos[0] - XcT)*(pos[0] - XcT) + (pos[1] - YcT)*(pos[1] - YcT));

		// check if photon is out of current column
		repeat = 0;
		ctr1 = 0;

		while( (R1 > (R-1e-5)) && (repeat < 10) && (ctr1 < 10) ) // R1 > R1-some small value..because of single precision errors that comparison with R may generate
		{

			// store current position
			temp_pos[0] = pos[0];
			temp_pos[1] = pos[1];

			// move particle by 0.5 um in the incident direction
			pos[0] = pos[0] + stepsize*(-dcos[0]);
			pos[1] = pos[1] + stepsize*(-dcos[1]);

			R2 = sqrt((pos[0] - XcT)*(pos[0] - XcT) + (pos[1] - YcT)*(pos[1] - YcT));

			if(R2 > R1) // means the photon is moving farther away from the column 
			{	    // this can happen if the stepsize if too big and the photon passes through the column and gets out on other side.

				// move it back to previous position, reduce the stepsize and try moving it again.
				pos[0] = temp_pos[0];
				pos[1] = temp_pos[1];

				stepsize = stepsize/2.0f;
				ctr1++;
			}
			else
			{
				R1 = R2;
				repeat++;
			}


		}

		// kill the particle if still outside the column
		if(R1 > (R-1e-5))
		 {

			d = -99.0f;
			num_outofcolT++;
			goto distexit;
		 }

		// condition to check that pos is within detector boundaries - if true, particle LOST
		if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) || (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
			{
				d = -99.0f;
				num_lostT++;
				goto distexit;
			}

	distexit:
	 return d;
	}	// C dist_to_surface function ends
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// calculate the directional cosines of the reflected vector
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ int boundary_analysisT(float3 *normal, float3 *pos, float3 *dcos, int2* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, 
	float top_absfrac, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, float *XcT, 
	float *YcT, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float d_max, float sensorRefl)
	{
		float3 dcos_temp = {0.0f};
		float3 old_normal = {0.0f};
		float3 old_dcos = {0.0f};
		float Pr = 0.0f, Pt = 0.0f;	// Prob. of reflection and transmission
		float theta1 = 0.0f, theta2 = 0.0f;
		float temp_norm = 0.0f;
		float mag = 0.0f;
		float rr_rnd = 0.0f;
		float theta_rnd = 0.0f;
		float angle_oldN_R = 0.0f;
		float newdepthT = 0.0f;
		float cos_newangle = 0.0f;
		float newangle = 0.0f;
		float cct1 = 0.0f;	// columnar crosstalk
		int trans_flag = 0.0f;
		int flag_abs = 0;		// flag - particle got absorbed at top surface or exited during the transmission to another boundary
		int flag_call_transmit = 1;	// flag - particle is going to move within a column (flag = 0) [call isotropic()] or between columns (flag = 1) [call transmit()] 
		int reperturb_ctr = 0;
		int flagCCT = 0;	// flag to indicate in transmit() that the photon needs to cross over
		int theta1ctr=0;	// counter for theta1 > 90 degrees (max resampling 100 times)
		int oldN_Rctr=0;	// counter for angle_oldN_R > 90 degrees (max resampling 25 times)
		int newnormalctr=0;

		// determine the coordinates of normal
		if ( (fabs(pos->z - (float)(H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	// reached top surface and dir. cosine in z-direction is positive
		{
	
			// top surface absorption - using absorption coefficient 'top_absfrac'
			if ( (top_absfrac > 0.0f) && (ranecuT(seed) < top_absfrac) )	// particle gets absorbed
			{
				flag_abs = 1;
				atomicAdd(&num_abs_topT, 1);
			}
			else
			{

				// specular reflector		
				normal->x = 0.0f;
				normal->y = 0.0f;
				normal->z = -1.0f;

				// assign new directional cosines
				dcos->z = -fabs((ranecuT(seed) * 2.0f) - 1.0f);
				rr_rnd = sqrt(1.0f - dcos->z*dcos->z);
				theta_rnd = ranecuT(seed)*twopipen;	

				dcos->x=rr_rnd*cos(theta_rnd);
				dcos->y=rr_rnd*sin(theta_rnd);
	
				flag_abs = 0;
			}

		}	
		else 	// compute the normal and check if gets reflected or transmitted
		{	
			// Columnar crosstalk
			newdepthT = H*0.2f;	// top 20% depth CCT=1. considering CsI layer only. NO organic polymer coating.
			
			if( (pos->z <= H/2.0f) && (pos->z >= (H/2.0f - newdepthT)) )	// top 20% - 100% cct
			{
				cct1 = 1.0f;
			}
			else if( (pos->z < (H/2.0f - newdepthT)) && (pos->z >= 0.0f) )  // from 20% depth to 50% - linear 100% to 50% 
			{
				cct1 = (pos->z/(2.0f*(H/2.0f - newdepthT))) + 0.5;	
			}
			else if( (pos->z < 0.0f) && (pos->z >= (-H/2.0f)) ) // bottom 50% to -H/2 - 50% to 100% CCT
			{
				cct1 = ( (pos->z - (-H/2.0f))/(2.0f * (-H/2.0f)) ) + 1.0 ;
			}
	
			if(ranecuT(seed) < cct1)		// columnar cross talk occurs
			{

				// photon crosses over to new column with random orientation. dcos do not change.
				flagCCT = 1;

				trans_flag = transmitT(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, mytid, num_detected_primary, d_max, sensorRefl, flagCCT);

				if (trans_flag == 1)		// photon exited
					flag_abs = 1;
				else if (trans_flag == 0)
				{
					// calculate new column's center coordinates
					*XcT = (float)( pos->x + R*(-normal->x) );
					*YcT = (float)( pos->y + R*(-normal->y) );

					flag_abs = 0;
				}
			}
			else
			{
			prpt:

				// within the column
				if(flag_call_transmit == 1)			// photon is currrently within a column with center Xc,Yc
				{
					mag = sqrt( (((*XcT)-pos->x) * ((*XcT)-pos->x)) + (((*YcT)-pos->y) * ((*YcT)-pos->y)) );
					normal->x = ((*XcT)-pos->x)/mag;
					normal->y = ((*YcT)-pos->y)/mag;
					normal->z = 0.0f;
		
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);	// new normal for rough surface

					flag_abs = 0;
				}
				// outside the column
				else if (flag_call_transmit == 0)		// photon is currently between columns and has not entered any column yet. New normal is sampled in the transmit(), so do not calculate normal here.
				{
					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					*XcT = (float)( pos->x + R*(-normal->x) );
					*YcT = (float)( pos->y + R*(-normal->y) );

				        flag_abs = 0;
				}

				// -dcos -> inverted the incident vector to get the smaller angle, else would have to do angle = 180-angle
				dcos_temp.x = -dcos->x;
				dcos_temp.y = -dcos->y;
				dcos_temp.z = -dcos->z;

				old_normal.x = normal->x;
				old_normal.y = normal->y;
				old_normal.z = normal->z;
	
				old_dcos.x = dcos->x;
				old_dcos.y = dcos->y;
				old_dcos.z = dcos->z;

			reperturb:
				normal->x = old_normal.x;
				normal->y = old_normal.y;
				normal->z = old_normal.z;

				dcos->x = old_dcos.x;
				dcos->y = old_dcos.y;
				dcos->z = old_dcos.z;

				dcos_temp.x = -dcos->x;
				dcos_temp.y = -dcos->y;
				dcos_temp.z = -dcos->z;

				if( (flag_call_transmit == 1) && (reperturb_ctr != 0) )
				 {
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);	// new normal for rough surface
				 }
				if( (flag_call_transmit == 0) && (reperturb_ctr != 0) )
				 {
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);	// the sampled normal for transmitted photon needs to be perturbed

					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					*XcT = (float)( pos->x + R*(-normal->x) );
					*YcT = (float)( pos->y + R*(-normal->y) );

				 }
			

				// Using Snell's law, calculate theta1 (angle between normal and reflected) and theta2 (angle between normal and transmitted)
			no_perturbation:
				theta1 = dot_productT(&dcos_temp,normal);		// cosine of angle between incident in opposite direction and normal (in radians)


				if ( (theta1 > 1.0f) || (theta1 < 0.0f) )	// if incidence angle > 1.57 radian or < 0 radian, then recalculate normal
				{
					// if particle was transmitted, then new normal has to be sampled again
					if(flag_call_transmit == 0)
					{
					mynewnormal:
						normal->x = dcos_temp.x;		// invert dcos of incident vector
						normal->y = dcos_temp.y;
						normal->z = dcos_temp.z;

						RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

						mag = sqrt(normal->x*normal->x + normal->y*normal->y);

						// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
						normal->z = 0.0f;			// normal_z of a cylinder is always zero
						normal->x = normal->x/mag;		// re-normalize
						normal->y = normal->y/mag;

						// perturb the normal according to Beta
						if(beta > 0.0f)
							RoughSurfaceT(normal, seed, beta);

						// find the angle between Normal and -Dcos
						cos_newangle = dot_productT(&dcos_temp, normal);
						newangle = acosf(cos_newangle);

						if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
						{						// keep looping until get a theta within 90 degrees
							if(newnormalctr < 100)
							{
								newnormalctr++;
								goto mynewnormal;			
							}
							else // kill it
							{
								atomicAdd(&num_theta1T,1);
								flag_abs = 1;
								newnormalctr = 0;
								goto baexit;
							}
						}

					}
		
					if(theta1ctr < 100)	// resample max 100 times	
					{
						theta1ctr++;
						goto prpt;
					}
					else	// kill it
					{
						atomicAdd(&num_theta1T,1);
						flag_abs = 1;
						theta1ctr = 0;
						goto baexit;
					}
				}
				else
					theta1 = acosf(theta1);
		

				// check for conditions where photon can only reflect
				if (flag_call_transmit == 1)	// only valid when photon within the column and can transmit outside the column. asin(n1/n2) -> nan
				{
					if (theta1 > asin(n2/n1))	// critical angle condition for TIR
					{
						Pr = 1.0f;		// TIR occurs
						Pt = 0.0f;
					}
			       		else if ( theta1 < epsilon ) 	// theta1 ~= 0, then always reflect
					{
				        	theta1 = 0.00042;       // make theta1 a very smal number, to avoid getting nan probabilities
					        theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians

					        // Using Fresnel's law, compute probability of reflection and transmission 
					        Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
				        	Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else    // the ray will transmit
					{
				        	theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians
					        
					        // Using Fresnel's law, compute probability of reflection and transmission 
				        	Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
					        Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}
				else if (flag_call_transmit == 0)
				{		
					if((n1/n2) < 1.57f)	// TIR can occur
					{
						if (theta1 > asin(n1/n2))	// critical angle condition for TIR
						{
							Pr = 1.0f;		// TIR occurs
							Pt = 0.0f;
						}
					}
					else if ( theta1 < epsilon )	// theta1 ~= 0, then always reflect
					{
						theta1 = 0.00042;	// make theta1 a very smal number, to avoid getting nan probabilities
						theta2 = asinf((float)(n1/n2)*sin(theta1)); 	// refracted/transmitted angle in radians
	
						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else	// the ray will transmit
					{
						theta2 = asinf((float)(n2/n1)*sin(theta1));

						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}


				// normalize Pr and Pt
				temp_norm = Pr + Pt;
				Pr = Pr/temp_norm;
				Pt = Pt/temp_norm;


				if(ranecuT(seed) < Pr)				// reflection
				{
					trans_dir_cosT(dcos, normal, theta1, theta2, 0, mytid, info);


					// condition to check that reflected vector is within 90 degrees from original normal
					angle_oldN_R = dot_productT(&old_normal, dcos);
					angle_oldN_R = acosf(angle_oldN_R);


					if (angle_oldN_R > 1.57f) // > 90 degrees, reperturb the normal
					{
						reperturb_ctr++;

						if(reperturb_ctr < 4)	// maximum 3 times reperturb, else calculate using smooth surface normal (old_normal)
							goto reperturb;
						else
						{
							normal->x = old_normal.x;
							normal->y = old_normal.y;
							normal->z = old_normal.z;

							dcos->x = old_dcos.x;
							dcos->y = old_dcos.y;
							dcos->z = old_dcos.z;

							dcos_temp.x = -dcos->x;
							dcos_temp.y = -dcos->y;
							dcos_temp.z = -dcos->z;

							reperturb_ctr = 0;
				
							if(oldN_Rctr < 25)	// resample max 100 times (25 * reperturb 4 times)	
							{
								oldN_Rctr++;
								goto no_perturbation;
							}
							else	// kill it
							{
								atomicAdd(&num_theta1T,1);
								flag_abs = 1;
								oldN_Rctr = 0;
								goto baexit;
							}

						}

					}

					if (flag_call_transmit == 0)		// it is reflecting between columns, so need to calculate distance using transmit()
					{
						flagCCT = 0;	

						trans_flag = transmitT(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, mytid, num_detected_primary, d_max, sensorRefl, flagCCT);

						if (trans_flag == 1)		// photon exited
							flag_abs = 1;
						else if (trans_flag == 0)
							goto prpt;				
					}
				}
				else						// transmission
				{
					trans_dir_cosT(dcos, normal, theta1, theta2, 1, mytid, info);

					if (flag_call_transmit == 1)		// photon travels between columns
					{
						flag_call_transmit = 0;
						flagCCT = 0;

						trans_flag = transmitT(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, mytid, num_detected_primary, d_max, sensorRefl, flagCCT);

						if (trans_flag == 1)		// particle exited
							flag_abs = 1;
						else if (trans_flag == 0)	// hits a column
							goto prpt;		// check again to see if it gets reflected or transmitted
					}			
				}
		
			} // else 'prpt ends
		
		} // main else ends

	baexit:
	   return flag_abs;

	}	// CUDA boundary analysis function ends
#else

	int boundary_analysisT(float *normal, float *pos, float *dcos, int* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, 
	float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, float d_max, 		float sensorRefl, int ydim, int *h_num_detected_prim)
	{

		float dcos_temp[3] = {0.0f};
		float old_normal[3] = {0.0f};
		float old_dcos[3] = {0.0f};
		float Pr = 0.0f, Pt = 0.0f;	// Prob. of reflection and transmission
		float theta1 = 0.0f;
		float theta2 = 0.0f;
		float temp_norm = 0.0f;
		float mag = 0.0f;
		float rr_rnd = 0.0f;
		float theta_rnd = 0.0f;
		float newdepthT = 0.0f;
		float angle_oldN_R = 0.0f;
		float cos_newangle = 0.0f;
		float newangle = 0.0f;
		float cct1 = 0.0f;	// columnar cross talk
		int reperturb_ctr = 0;
		int trans_flag = 0.0f;
		int flag_abs = 0;		// flag - particle got absorbed at top surface or exited during the transmission to another boundary
		int flag_call_transmit = 1;	// flag - particle is going to move within a column (flag = 1) [call isotropic()] or between columns (flag = 0) [call transmit()]
		int theta1ctr=0;	// counter for theta1 > 90 degrees (max resampling 100 times)
		int oldN_Rctr=0;	// counter for angle_oldN_R > 90 degrees (max resampling 25 times)
		int newnormalctr=0;
		int flagCCT = 0;	// flag to indicate in transmit() that the photon needs to cross over

		// determine the coordinates of normal
		if ( (fabs(pos[2] - (float)(H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	// reached top surface and dir. cosine in z-direction is positive
		{
	
			// top surface absorption - using absorption coefficient 'top_absfrac'
			if ( (top_absfrac > 0.0f) && (ranecuT(seed) < top_absfrac) )	// particle gets absorbed
			{
				flag_abs = 1;
				num_abs_topT++;
			}
			else
			{

				// specular reflector		
				normal[0] = 0.0f;
				normal[1] = 0.0f;
				normal[2] = -1.0f;

				// assign new directional cosines
				dcos[2] = -fabs((ranecuT(seed) * 2.0f) - 1.0f);
				rr_rnd = sqrt(1.0f - dcos[2]*dcos[2]);
				theta_rnd = ranecuT(seed)*twopipen;	

				dcos[0]=rr_rnd*cos(theta_rnd);
				dcos[1]=rr_rnd*sin(theta_rnd);
	
				flag_abs = 0;
			}

		}	
		else 	// compute the normal and check if gets reflected or transmitted
		{	

			// Columnar crosstalk
			newdepthT = H*0.2f;	// top 20% depth CCT=1. considering CsI layer only. NO organic polymer coating.
		
			if( (pos[2] <= H/2.0f) && (pos[2] >= (H/2.0f - newdepthT)) )	// top 20% - 100% cct
			{
				cct1 = 1.0f;
			}
			else if( (pos[2] < (H/2.0f - newdepthT)) && (pos[2] >= 0.0f) )  // from 20% depth to 50% - linear 100% to 50% 
			{
				cct1 = (pos[2]/(2.0f*(H/2.0f - newdepthT))) + 0.5;	
			}
			else if( (pos[2] < 0.0f) && (pos[2] >= (-H/2.0f)) ) // bottom 50% to (-H/2 - 4 um polymer) - 50% to 100% CCT
			{
				cct1 = ( (pos[2] - (-H/2.0f))/(2.0f * (-H/2.0f)) ) + 1.0 ;
			}
	
			if(ranecuT(seed) < cct1)		// columnar cross talk occurs
			{

				// photon crosses over to new column with random orientation. dcos do not change.
				flagCCT = 1;

				trans_flag = transmitT(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, flagCCT, h_num_detected_prim);

				if (trans_flag == 1)		// photon exited
					flag_abs = 1;
				else if (trans_flag == 0)
				{
					// calculate new column's center coordinates
					XcT = (float)( pos[0] + R*(-normal[0]) );
					YcT = (float)( pos[1] + R*(-normal[1]) );

					flag_abs = 0;
				}
			}
			else
			{

			prpt:

				// within the column
				if(flag_call_transmit == 1)			// photon is currrently within a column with center Xc,Yc
				{
					mag = sqrt( (((XcT)-pos[0]) * ((XcT)-pos[0])) + (((YcT)-pos[1]) * ((YcT)-pos[1])) );
					normal[0] = ((XcT)-pos[0])/mag;
					normal[1] = ((YcT)-pos[1])/mag;
					normal[2] = 0.0f;
		
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);	// new normal for rough surface

					flag_abs = 0;
				}
				// outside the column
				else if (flag_call_transmit == 0)		// photon is currently between columns and has not entered any column yet. New normal is sampled in the transmit(), so do not calculate normal here.
				{
					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					XcT = (float)( pos[0] + R*(-normal[0]) );
					YcT = (float)( pos[1] + R*(-normal[1]) );

				        flag_abs = 0;
				}

				// -dcos -> inverted the incident vector to get the smaller angle, else would have to do angle = 180-angle
				dcos_temp[0] = -dcos[0];
				dcos_temp[1] = -dcos[1];
				dcos_temp[2] = -dcos[2];

				old_normal[0] = normal[0];
				old_normal[1] = normal[1];
				old_normal[2] = normal[2];
	
				old_dcos[0] = dcos[0];
				old_dcos[1] = dcos[1];
				old_dcos[2] = dcos[2];


			reperturb:
				normal[0] = old_normal[0];
				normal[1] = old_normal[1];
				normal[2] = old_normal[2];

				dcos[0] = old_dcos[0];
				dcos[1] = old_dcos[1];
				dcos[2] = old_dcos[2];

				dcos_temp[0] = -dcos[0];
				dcos_temp[1] = -dcos[1];
				dcos_temp[2] = -dcos[2];

				if( (flag_call_transmit == 1) && (reperturb_ctr != 0) )
				 {
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);	// new normal for rough surface
				 }
				if( (flag_call_transmit == 0) && (reperturb_ctr != 0) )
				 {
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);	// the sampled normal for transmitted photon needs to be perturbed

					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					XcT = (float)( pos[0] + R*(-normal[0]) );
					YcT = (float)( pos[1] + R*(-normal[1]) );

				 }
			

				// Using Snell's law, calculate theta1 (angle between normal and reflected) and theta2 (angle between normal and transmitted)
			no_perturbation:
				theta1 = dot_productT(dcos_temp,normal);		// cosine of angle between incident in opposite direction and normal (in radians)

				if ( (theta1 > 1.0f) || (theta1 < 0.0f) )	// if incidence angle > 1.57 radian or < 0 radian, then recalculate normal
				{
					// if particle was transmitted, then new normal has to be sampled again
					if(flag_call_transmit == 0)
					{
					mynewnormal:
						normal[0] = dcos_temp[0];		// invert dcos of incident vector
						normal[1] = dcos_temp[1];
						normal[2] = dcos_temp[2];

						RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

						mag = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

						// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
						normal[2] = 0.0f;			// normal_z of a cylinder is always zero
						normal[0] = normal[0]/mag;		// re-normalize
						normal[1] = normal[1]/mag;

						// perturb the normal according to Beta
						if(beta > 0.0f)
							RoughSurfaceT(normal, seed, beta);

						// find the angle between Normal and -Dcos
						cos_newangle = dot_productT(dcos_temp, normal);
						newangle = acosf(cos_newangle);

						if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
						{
							if(newnormalctr < 100)
							{
								newnormalctr++;
								goto mynewnormal;			
							}
							else // kill it
							{
								num_theta1T++;
								flag_abs = 1;
								newnormalctr = 0;
								goto baexit;
							}				// keep looping until get a theta within 90 degrees
						}

					}

					if(theta1ctr < 100)	
					{
						theta1ctr++;
						goto prpt;
					}
					else	// kill it
					{
						num_theta1T++;
						flag_abs = 1;
						theta1ctr = 0;
						goto baexit;
					}
				}
				else
					theta1 = acosf(theta1);
		

				// check for conditions where photon can only reflect
				if (flag_call_transmit == 1)	// only valid when photon within the column and can transmit outside the column. asin(n1/n2) -> nan
				{
					if (theta1 > asin(n2/n1))	// critical angle condition for TIR
					{
						Pr = 1.0f;		// TIR occurs
						Pt = 0.0f;
					}
			       		else if ( theta1 < epsilon ) 	// theta1 ~= 0, then always reflect
					{
				        	theta1 = 0.00042;       // make theta1 a very smal number, to avoid getting nan probabilities
					        theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians

					        // Using Fresnel's law, compute probability of reflection and transmission 
					        Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
				        	Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else    // the ray will transmit
					{
				        	theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians
					        
					        // Using Fresnel's law, compute probability of reflection and transmission 
				        	Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
					        Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}
				else if (flag_call_transmit == 0)
				{		
					if((n1/n2) < 1.57f)	// TIR can occur
					{
						if (theta1 > asin(n1/n2))	// critical angle condition for TIR
						{
							Pr = 1.0f;		// TIR occurs
							Pt = 0.0f;
						}
					}
					else if ( theta1 < epsilon )	// theta1 ~= 0, then always reflect
					{
						theta1 = 0.00042;	// make theta1 a very smal number, to avoid getting nan probabilities
						theta2 = asinf((float)(n1/n2)*sin(theta1)); 	// refracted/transmitted angle in radians
	
						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else	// the ray will transmit
					{
						theta2 = asinf((float)(n2/n1)*sin(theta1));

						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}


				// normalize Pr and Pt
				temp_norm = Pr + Pt;
				Pr = Pr/temp_norm;
				Pt = Pt/temp_norm;


				if(ranecuT(seed) < Pr)				// reflection
				{
					trans_dir_cosT(dcos, normal, theta1, theta2, 0, info);


					// condition to check that reflected vector is within 90 degrees from original normal
					angle_oldN_R = dot_productT(old_normal, dcos);
					angle_oldN_R = acosf(angle_oldN_R);


					if (angle_oldN_R > 1.57f) // > 90 degrees, reperturb the normal
					{
						reperturb_ctr++;

						if(reperturb_ctr < 4)	// maximum 3 times reperturb, else calculate using smooth surface normal (old_normal)
							goto reperturb;
						else
						{
							normal[0] = old_normal[0];
							normal[1] = old_normal[1];
							normal[2] = old_normal[2];

							dcos[0] = old_dcos[0];
							dcos[1] = old_dcos[1];
							dcos[2] = old_dcos[2];

							dcos_temp[0] = -dcos[0];
							dcos_temp[1] = -dcos[1];
							dcos_temp[2] = -dcos[2];

							reperturb_ctr = 0;


							if(oldN_Rctr < 25)	// max resample 25 times (25*4reperturb = 100 times)
							{
								oldN_Rctr++;
								goto no_perturbation;
							}
							else	// kill it
							{
								num_theta1T++;
								flag_abs = 1;
								oldN_Rctr = 0;
								goto baexit;
							}
						}
					}

					if (flag_call_transmit == 0)		// it is reflecting between columns, so need to calculate distance using transmit()
					{
						flagCCT = 0;

						trans_flag = transmitT(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, flagCCT, h_num_detected_prim);

						if (trans_flag == 1)		// photon exited
							flag_abs = 1;
						else if (trans_flag == 0)
							goto prpt;				
					}
				}
				else						// transmission
				{
					trans_dir_cosT(dcos, normal, theta1, theta2, 1, info);

					if (flag_call_transmit == 1)		// photon exits current column
					{
						flag_call_transmit = 0;
						flagCCT = 0;

						trans_flag = transmitT(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, flagCCT, h_num_detected_prim);

						if (trans_flag == 1)		// particle exited
							flag_abs = 1;
						else if (trans_flag == 0)	// hits a column
							goto prpt;		// check again to see if it gets reflected or transmitted
					}			
				}

			} //else prpt ends

		}	// main else ends

	baexit:
	   return flag_abs;

	}	// C boundary analysis function ends
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Photon gets transmitted, calculate the new position where it hits next column or boundary
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ int transmitT(float3 *pos, float3 *dcos, float3 *normal, int2* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min, 
	int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, size_t pitch, struct start_info *info, int mytid, 
	int *num_detected_primary, float d_max, float sensorRefl, int flagCCT)
	{
		float3 temp_pos = {0.0f};
		float d_nextCol = 0.0f;		// distance to next column	
		float d_top = 0.0f;		// distance to top surface
		float d_bottom = 0.0f;		// distance to bottom surface
		int particle_exit = 0;		// flag to indicate if photon enters another column or gets lost/detected/absorbed
		float newangle = 0.0f;
		float cos_newangle = 0.0f;
		float3 temp_dcos = {0.0f};
		float rr_rnd = 0.0f, theta_rnd = 0.0f;
		float tmp_deno = 0.0f;
		int iii = 0, jjj = 0;
		int reflbtm = 0;

		int newnormalctr=0;	// counter for resampling rough normal to new column (max 100 times, else kill it)
		int newnormalctr2=0;

		temp_pos.x = pos->x;
		temp_pos.y = pos->y;
		temp_pos.z = pos->z;

		temp_dcos.x = -dcos->x;
		temp_dcos.y = -dcos->y;
		temp_dcos.z = -dcos->z;


		if(flagCCT == 1)	// CCT occurs
		{
			// no change in dcos. d_nextcol = 0. new column has random orientation.
			newnormal1:
				normal->x = temp_dcos.x;		// invert dcos of incident vector
				normal->y = temp_dcos.y;
				normal->z = temp_dcos.z;

				RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

				tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

				// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
				normal->z = 0.0f;			// normal_z of a cylinder is always zero
				normal->x = normal->x/tmp_deno;		// re-normalize
				normal->y = normal->y/tmp_deno;

				// perturb the normal according to Beta
				if(beta > 0.0f)
					RoughSurfaceT(normal, seed, beta);

				// find the angle between Normal and -Dcos
				cos_newangle = dot_productT(&temp_dcos, normal);
				newangle = acosf(cos_newangle);

				if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
				 {						// keep looping until get a theta within 90 degrees
					if(newnormalctr < 100)
					{
						newnormalctr++;
						goto newnormal1;			
					}
					else // kill it
					{
						atomicAdd(&num_theta1T,1);
						particle_exit = 1;
						newnormalctr = 0;
						goto exitnow;
					}
						
				 }

				particle_exit = 0;
		}
		else
		{

			// sample distance uniformly between d_min and d_max to next column
			d_nextCol = ranecuT(seed) * (d_max - d_min) + d_min;

			// compute the new position of the photon. 
			pos->x = temp_pos.x + dcos->x * d_nextCol;
			pos->y = temp_pos.y + dcos->y * d_nextCol;
			pos->z = temp_pos.z + dcos->z * d_nextCol;

			// calculate distance to top and bottom surface: if d_top < d_nextCol then photon should reflect from the top surface; else if d_bottom < d_nextCol, photon should get detected.
			d_top = ((H/2.0f) - temp_pos.z)/dcos->z;
			d_bottom  = ((-H/2.0f) - temp_pos.z)/dcos->z;

			// condition to check that pos is within detector boundaries - if true, particle LOST
			if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) )
			{
				atomicAdd(&num_lostT, 1);
				particle_exit = 1;
				goto exitnow;
			}

			if ( (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
				{
					if( (d_top < d_nextCol) && (d_top > epsilon) )
					{
						pos->x = temp_pos.x + dcos->x * d_top;
						pos->y = temp_pos.y + dcos->y * d_top;
						pos->z = H/2.0f;
				
						atomicAdd(&photon_distanceT, d_top);
						particle_exit = 0;			
					}
					else if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )
					{
						pos->x = temp_pos.x + dcos->x * d_bottom;
						pos->y = temp_pos.y + dcos->y * d_bottom;
						pos->z = -H/2.0f;

						atomicAdd(&photon_distanceT, d_bottom);

						// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
						if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
						{
		
							particle_exit = 0;

							// normal pointing (0,0,1)
							normal->x = 0.0f; normal->y = 0.0f; normal->z = 1.0f;

							// obtain reflected dcos from the bottom (specular reflection; 
							// bottom surface is smooth, so no need to perturb the normal)
							// this condition is called only when photon hits the bottom surface OUTSIDE any column
							trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, mytid, info);	// reflection only so refl_theta,trans_theta = 0

							// sample new distance and place new column
							reflbtm = refl_bottomT(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

							if(reflbtm == 1)
							{
								particle_exit = 1;
								goto exitnow;
							}

							// if it hits top surface after reflecting back
							if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	
							{
								goto mytopsurface;
							}

					

						}
						else
						{
							particle_exit = 1;	
							atomicAdd(&num_detectT, 1);


							iii = floor((pos->x-lbound_x)/pixelsize);	// determine pixel number in x and y direction
							jjj = floor((pos->y-lbound_y)/pixelsize);

							// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
							if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
							 {	
								unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + iii * pitch);
								atomicAdd(&current_img[jjj],1);
							 }
		
							atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);// start array from 0.str_histnum starts from 1
							
							goto exitnow;
						}	
					}
					else
					{
						atomicAdd(&num_lostT, 1);
						particle_exit = 1;
						goto exitnow;
					}
				}
			else
				atomicAdd(&photon_distanceT, d_nextCol);		// add distance travelled to global variable


			// sample new normal to determine orientation of new column.
		newnormal:
			normal->x = temp_dcos.x;		// invert dcos of incident vector
			normal->y = temp_dcos.y;
			normal->z = temp_dcos.z;

			RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

			tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

			// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
			normal->z = 0.0f;			// normal_z of a cylinder is always zero
			normal->x = normal->x/tmp_deno;		// re-normalize
			normal->y = normal->y/tmp_deno;

			// perturb the normal according to Beta
			if(beta > 0.0f)
				RoughSurfaceT(normal, seed, beta);

			// find the angle between Normal and -Dcos
			cos_newangle = dot_productT(&temp_dcos, normal);
			newangle = acosf(cos_newangle);

			if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
			 {						// keep looping until get a theta within 90 degrees
				if(newnormalctr < 100)
				{
					newnormalctr++;
					goto newnormal;			
				}
				else // kill it
				{
					atomicAdd(&num_theta1T,1);
					particle_exit = 1;
					newnormalctr = 0;
					goto exitnow;
				}
						
			 }
	
			// check if the photon enters another column or got lost (hit detector side)/ reflected (detector top)/ detected (detector bottom)
	
			// hit side of detector?
			if ( (fabs(pos->x-0.0f) < epsilon) || (fabs(pos->x-xdetector) < epsilon) || (fabs(pos->y-0.0f) < epsilon) || (fabs(pos->y-ydetector) < epsilon) )
			{
				atomicAdd(&num_lostT, 1);
				particle_exit = 1;
				goto exitnow;
			}
			

			// hit top?
		     mytopsurface:

			if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	// gets specularly reflected or absorbed
			{
				normal->x = 0.0f;
				normal->y = 0.0f;
				normal->z = -1.0f;

				// top surface absorption - using absorption coefficient 'top_absfrac'
				if ( (top_absfrac > 0.0f) && (ranecuT(seed) < top_absfrac) )	// particle gets absorbed
				{
					atomicAdd(&num_abs_topT, 1);
					particle_exit = 1;
					goto exitnow;
				}
				else
				{

					// assign new directional cosines
					dcos->z = -fabs((ranecuT(seed) * 2.0f) - 1.0f);
					rr_rnd = sqrt(1.0f - dcos->z*dcos->z);
					theta_rnd = ranecuT(seed)*twopipen;	
	
					dcos->x=rr_rnd*cos(theta_rnd);
					dcos->y=rr_rnd*sin(theta_rnd);

					temp_pos.x = pos->x;
					temp_pos.y = pos->y;
					temp_pos.z = pos->z;

					temp_dcos.x = -dcos->x;
					temp_dcos.y = -dcos->y;
					temp_dcos.z = -dcos->z;

					// sample distance uniformly between d_min and d_max to next column
					d_nextCol = ranecuT(seed) * (d_max - d_min) + d_min;

					// calculate distance to bottom surface: if d_bottom < d_nextCol, photon should get detected.
					d_bottom  = ((-H/2.0f) - temp_pos.z)/dcos->z;

					// compute the new position of the photon. 
					pos->x = temp_pos.x + dcos->x * d_nextCol;
					pos->y = temp_pos.y + dcos->y * d_nextCol;
					pos->z = temp_pos.z + dcos->z * d_nextCol;

					// condition to check that pos is within detector boundaries - if true, particle LOST
					if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) )
					{
						atomicAdd(&num_lostT, 1);
						particle_exit = 1;
						goto exitnow;
					}

					if ( (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
						{
							if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )
							{
								pos->x = temp_pos.x + dcos->x * d_bottom;
								pos->y = temp_pos.y + dcos->y * d_bottom;
								pos->z = -H/2.0f;

								atomicAdd(&photon_distanceT, d_bottom);

								// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
								if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
								{
									particle_exit = 0;

									// normal pointing (0,0,1)
									normal->x = 0.0f; normal->y = 0.0f; normal->z = 1.0f;

									// obtain reflected dcos from the bottom (specular reflection; 
									// bottom surface is smooth, so no need to perturb the normal)
									// this condition is called only when photon hits the bottom surface OUTSIDE any column
									trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, mytid, info);	// reflection only so refl_theta,trans_theta = 0

									// sample new distance and place new column
									reflbtm = refl_bottomT(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

									if(reflbtm == 1)
									{
										particle_exit = 1;
										goto exitnow;
									}

									// if it hits top surface after reflecting back
									if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	
									{
										goto mytopsurface;
									}

								}
								else
								{
										particle_exit = 1;
										atomicAdd(&num_detectT, 1);

										iii = floor((pos->x-lbound_x)/pixelsize);// determine pixel number in x and y direction
										jjj = floor((pos->y-lbound_y)/pixelsize);

										// if photon gets detected within lower and upper bounds: accumulate signal contribution
										if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
										 {	
											unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + iii * pitch);
											atomicAdd(&current_img[jjj],1);
										 }

										 atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);

										goto exitnow;	
								}	
							}
							else
							{
								atomicAdd(&num_lostT, 1);
								particle_exit = 1;
								goto exitnow;
							}
						}
					else
						atomicAdd(&photon_distanceT, d_nextCol);		// add distance travelled to global variable

					// sample new normal to determine orientation of new column.
			  newnormal_TOP:
					normal->x = temp_dcos.x;		// invert dcos of incident vector
					normal->y = temp_dcos.y;
					normal->z = temp_dcos.z;

					RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

					tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

					// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
					normal->z = 0.0f;			// normal_z of a cylinder is always zero
					normal->x = normal->x/tmp_deno;		// re-normalize
					normal->y = normal->y/tmp_deno;

					// perturb the normal according to Beta
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);

					// find the angle between Normal and -Dcos
					cos_newangle = dot_productT(&temp_dcos, normal);
					newangle = acosf(cos_newangle);

					if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
					 {
						if(newnormalctr2 < 100)	// resample max 100 times
						{
							newnormalctr2++;
							goto newnormal_TOP;			// keep looping until get a theta within 90 degrees		
						}
						else // kill it
						{
							atomicAdd(&num_theta1T,1);
							particle_exit = 1;
							newnormalctr2 = 0;
							goto exitnow;
						}
				
					 }
	
					particle_exit = 0;
				}
			}	// hit top ends
	

			// hit bottom? z of detector can be in the range (-H/2, H/2).
			if ( fabs(pos->z - (-H/2.0f)) < epsilon )	// gets detected
			{
				// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
				if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
				{
					particle_exit = 0;

					// normal pointing (0,0,1)
					normal->x = 0.0f; normal->y = 0.0f; normal->z = 1.0f;

					// obtain reflected dcos from the bottom (specular reflection; 
					// bottom surface is smooth, so no need to perturb the normal)
					// this condition is called only when photon hits the bottom surface OUTSIDE any column
					trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, mytid, info);	// reflection only so refl_theta,trans_theta = 0

					// sample new distance and place new column
					reflbtm = refl_bottomT(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

					if(reflbtm == 1)
					{
						particle_exit = 1;
						goto exitnow;
					}

					// if it hits top surface after reflecting back
					if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	
					{
						goto mytopsurface;
					}

				}
				else
				{
					atomicAdd(&num_detectT, 1);
					particle_exit = 1;

					iii = floor((pos->x-lbound_x)/pixelsize);	// determine pixel number in x and y direction
					jjj = floor((pos->y-lbound_y)/pixelsize);

					// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
					if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
					 {	
						unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + iii * pitch);
						atomicAdd(&current_img[jjj],1);
					 }

					atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);	// start array from 0.str_histnum starts from 1
					
					goto exitnow;
				}
			}

		} // else CCT ends
	
	exitnow:
	 return particle_exit;	
	}	// CUDA transmit function ends
#else
	int transmitT(float *pos, float *dcos, float *normal, int* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min, int pixelsize, 
	float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, float d_max, float sensorRefl, int ydim, int flagCCT,
	int *h_num_detected_prim)
	{

		float temp_pos[3] = {0.0f};
		float d_nextCol = 0.0f;		// distance to next column	
		float d_top = 0.0f;		// distance to top surface
		float d_bottom = 0.0f;		// distance to bottom surface
		int particle_exit = 0;		// flag to indicate if photon enters another column or gets lost/detected/absorbed
		float newangle = 0.0f;
		float cos_newangle = 0.0f;
		float temp_dcos[3] = {0.0f};
		float rr_rnd = 0.0f, theta_rnd = 0.0f;
		float tmp_deno = 0.0f;
		int iii = 0, jjj = 0;
		int reflbtm = 0;

		int newnormalctr = 0;
		int newnormalctr2 = 0;

		temp_pos[0] = pos[0];
		temp_pos[1] = pos[1];
		temp_pos[2] = pos[2];

		temp_dcos[0] = -dcos[0];
		temp_dcos[1] = -dcos[1];
		temp_dcos[2] = -dcos[2];

		if(flagCCT == 1)	// CCT occurs
		{
			// no change in dcos. d_nextcol = 0. new column has random orientation.
			newnormal1:
				normal[0] = temp_dcos[0];		// invert dcos of incident vector
				normal[1] = temp_dcos[1];
				normal[2] = temp_dcos[2];

				RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

				tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

				// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
				normal[2] = 0.0f;			// normal_z of a cylinder is always zero
				normal[0] = normal[0]/tmp_deno;		// re-normalize
				normal[1] = normal[1]/tmp_deno;

				// perturb the normal according to Beta
				if(beta > 0.0f)
					RoughSurfaceT(normal, seed, beta);

				// find the angle between Normal and -Dcos
				cos_newangle = dot_productT(temp_dcos, normal);
				newangle = acosf(cos_newangle);

				if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
				 {						// keep looping until get a theta within 90 degrees
					if(newnormalctr < 100)
					{
						newnormalctr++;
						goto newnormal1;			
					}
					else // kill it
					{
						num_theta1T++;
						particle_exit = 1;
						newnormalctr = 0;
						goto exitnow;
					}
						
				 }

				particle_exit = 0;
		}
		else
		{

			// sample distance uniformly between d_min and d_max to next column
			d_nextCol = ranecuT(seed) * (d_max - d_min) + d_min;

			// compute the new position of the photon. 
			pos[0] = temp_pos[0] + dcos[0] * d_nextCol;
			pos[1] = temp_pos[1] + dcos[1] * d_nextCol;
			pos[2] = temp_pos[2] + dcos[2] * d_nextCol;

			// calculate distance to top and bottom surface: if d_top < d_nextCol then photon should reflect from the top surface; else if d_bottom < d_nextCol, photon should get detected.
			d_top = ((H/2.0f) - temp_pos[2])/dcos[2];
			d_bottom  = ((-H/2.0f) - temp_pos[2])/dcos[2];

			// condition to check that pos is within detector boundaries - if true, particle LOST
			if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) )
			{
				num_lostT++;
				particle_exit = 1;
				goto exitnow;
			}

			if ( (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
				{
					if( (d_top < d_nextCol) && (d_top > epsilon) )
					{
						pos[0] = temp_pos[0] + dcos[0] * d_top;
						pos[1] = temp_pos[1] + dcos[1] * d_top;
						pos[2] = H/2.0f;
				
						photon_distanceT = photon_distanceT + d_top;
						particle_exit = 0;			
					}
					else if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )
					{
						pos[0] = temp_pos[0] + dcos[0] * d_bottom;
						pos[1] = temp_pos[1] + dcos[1] * d_bottom;
						pos[2] = -H/2.0f;

						photon_distanceT = photon_distanceT + d_bottom;

						// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
						if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
						{
		
							particle_exit = 0;

							// normal pointing (0,0,1)
							normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 1.0f;

							// obtain reflected dcos from the bottom (specular reflection; 
							// bottom surface is smooth, so no need to perturb the normal)
							// this condition is called only when photon hits the bottom surface OUTSIDE any column
							trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, info);	// reflection only so refl_theta,trans_theta = 0

							// sample new distance and place new column
							reflbtm = refl_bottomT(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

							if(reflbtm == 1)
							{
								particle_exit = 1;
								goto exitnow;
							}


							// if it hits top surface after reflecting back
							if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	
							{
								goto mytopsurface;
							}

						}
						else
						{
							particle_exit = 1;	
							num_detectT++;


							iii = floor((pos[0]-lbound_x)/pixelsize);	// determine pixel number in x and y direction
							jjj = floor((pos[1]-lbound_y)/pixelsize);

							// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
							if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
							 {	
								outputimage_.tempimageopt[iii][jjj]++;
							 }

							goto exitnow;
						}	
					}
					else
					{
						num_lostT++;
						particle_exit = 1;
						goto exitnow;
					}
				}
			else
				photon_distanceT = photon_distanceT + d_nextCol;		// add distance travelled to global variable


			// sample new normal to determine orientation of new column.
		newnormal:
			normal[0] = temp_dcos[0];		// invert dcos of incident vector
			normal[1] = temp_dcos[1];
			normal[2] = temp_dcos[2];

			RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

			tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

			// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
			normal[2] = 0.0f;			// normal_z of a cylinder is always zero
			normal[0] = normal[0]/tmp_deno;		// re-normalize
			normal[1] = normal[1]/tmp_deno;

			// perturb the normal according to Beta
			if(beta > 0.0f)
				RoughSurfaceT(normal, seed, beta);

			// find the angle between Normal and -Dcos
			cos_newangle = dot_productT(temp_dcos, normal);
			newangle = acosf(cos_newangle);

			if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
			 {						// keep looping until get a theta within 90 degrees
				if(newnormalctr < 100)
				{
					newnormalctr++;
					goto newnormal;			
				}
				else // kill it
				{
					num_theta1T++;
					particle_exit = 1;
					newnormalctr = 0;
					goto exitnow;
				}
						
			 }
	
			// check if the photon enters another column or got lost (hit detector side)/ reflected (detector top)/ detected (detector bottom)
	
			// hit side of detector?
			if ( (fabs(pos[0]-0.0f) < epsilon) || (fabs(pos[0]-xdetector) < epsilon) || (fabs(pos[1]-0.0f) < epsilon) || (fabs(pos[1]-ydetector) < epsilon) )
			{
				num_lostT++;
				particle_exit = 1;
				goto exitnow;
			}
		

			// hit top?
		     mytopsurface:

			if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	// gets specularly reflected or absorbed
			{
				normal[0] = 0.0f;
				normal[1] = 0.0f;
				normal[2] = -1.0f;

				// top surface absorption - using absorption coefficient 'top_absfrac'
				if ( (top_absfrac > 0.0f) && (ranecuT(seed) < top_absfrac) )	// particle gets absorbed
				{
					num_abs_topT++;
					particle_exit = 1;
					goto exitnow;
				}
				else
				{

					// assign new directional cosines
					dcos[2] = -fabs((ranecuT(seed) * 2.0f) - 1.0f);
					rr_rnd = sqrt(1.0f - dcos[2]*dcos[2]);
					theta_rnd = ranecuT(seed)*twopipen;	
	
					dcos[0]=rr_rnd*cos(theta_rnd);
					dcos[1]=rr_rnd*sin(theta_rnd);

					temp_pos[0] = pos[0];
					temp_pos[1] = pos[1];
					temp_pos[2] = pos[2];

					temp_dcos[0] = -dcos[0];
					temp_dcos[1] = -dcos[1];
					temp_dcos[2] = -dcos[2];

					// sample distance uniformly between d_min and d_max to next column
					d_nextCol = ranecuT(seed) * (d_max - d_min) + d_min;

					// calculate distance to bottom surface: if d_bottom < d_nextCol, photon should get detected.
					d_bottom  = ((-H/2.0f) - temp_pos[2])/dcos[2];

					// compute the new position of the photon. 
					pos[0] = temp_pos[0] + dcos[0] * d_nextCol;
					pos[1] = temp_pos[1] + dcos[1] * d_nextCol;
					pos[2] = temp_pos[2] + dcos[2] * d_nextCol;

					// condition to check that pos is within detector boundaries - if true, particle LOST
					if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) )
					{
						num_lostT++;
						particle_exit = 1;
						goto exitnow;
					}

					if ( (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
						{
							if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )
							{
								pos[0] = temp_pos[0] + dcos[0] * d_bottom;
								pos[1] = temp_pos[1] + dcos[1] * d_bottom;
								pos[2] = -H/2.0f;

								photon_distanceT = photon_distanceT + d_bottom;

								// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
								if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
								{

									particle_exit = 0;		
			
									// normal pointing (0,0,1)
									normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 1.0f;

									// obtain reflected dcos from the bottom (specular reflection; 
									// bottom surface is smooth, so no need to perturb the normal)
									// this condition is called only when photon hits the bottom surface OUTSIDE any column
									trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, info);	// reflection only so refl_theta,trans_theta = 0

									// sample new distance and place new column
									reflbtm = refl_bottomT(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

									if(reflbtm == 1)
									{
										particle_exit = 1;
										goto exitnow;
									}

									// if it hits top surface after reflecting back
									if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	
									{
										goto mytopsurface;
									}

								}
								else
								{
										particle_exit = 1;
										num_detectT++;

										iii = floor((pos[0]-lbound_x)/pixelsize);
										jjj = floor((pos[1]-lbound_y)/pixelsize);

										// if photon gets detected within lower and upper bounds: accumulate signal contribution
										if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
										 {	
											outputimage_.tempimageopt[iii][jjj]++;
										 }

										goto exitnow;	
								}	
							}
							else
							{
								num_lostT++;
								particle_exit = 1;
								goto exitnow;
							}
						}
					else
						photon_distanceT = photon_distanceT + d_nextCol;		// add distance travelled to global variable

					// sample new normal to determine orientation of new column.
			  newnormal_TOP:
					normal[0] = temp_dcos[0];		// invert dcos of incident vector
					normal[1] = temp_dcos[1];
					normal[2] = temp_dcos[2];

					RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

					tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

					// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
					normal[2] = 0.0f;			// normal_z of a cylinder is always zero
					normal[0] = normal[0]/tmp_deno;		// re-normalize
					normal[1] = normal[1]/tmp_deno;

					// perturb the normal according to Beta
					if(beta > 0.0f)
						RoughSurfaceT(normal, seed, beta);

					// find the angle between Normal and -Dcos
					cos_newangle = dot_productT(temp_dcos, normal);
					newangle = acosf(cos_newangle);

					if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
					 {
						if(newnormalctr2 < 100)	// resample max 100 times
						{
							newnormalctr2++;
							goto newnormal_TOP;			// keep looping until get a theta within 90 degrees		
						}
						else // kill it
						{
							num_theta1T++;
							particle_exit = 1;
							newnormalctr2 = 0;
							goto exitnow;
						}
				
					 }
	
					particle_exit = 0;
				}
			}	// hit top ends
	

			// hit bottom? z of detector can be in the range (-H/2, H/2).
			if ( fabs(pos[2] - (-H/2.0f)) < epsilon )	// gets detected
			{
				// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
				if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
				{
					particle_exit = 0;		
	
					// normal pointing (0,0,1)
					normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 1.0f;

					// obtain reflected dcos from the bottom (specular reflection; 
					// bottom surface is smooth, so no need to perturb the normal)
					// this condition is called only when photon hits the bottom surface OUTSIDE any column
					trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, info);	// reflection only so refl_theta,trans_theta = 0

					// sample new distance and place new column
					reflbtm = refl_bottomT(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

					if(reflbtm == 1)
					{
						particle_exit = 1;
						goto exitnow;
					}

					// if it hits top surface after reflecting back
					if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	
					{
						goto mytopsurface;
					}

				}
				else
				{
					num_detectT++;
					particle_exit = 1;

					iii = floor((pos[0]-lbound_x)/pixelsize);	// determine pixel number in x and y direction
					jjj = floor((pos[1]-lbound_y)/pixelsize);

					// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
					if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
					 {	
						outputimage_.tempimageopt[iii][jjj]++;
					 }

					goto exitnow;
				}
			}

		} // else CCT ends

	exitnow:
	 return particle_exit;	
	}	// C transmit function ends
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// photon reflects from sensor_plane or bottom surface, when in between columns. 
// Obtains the next column where it hits.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ int refl_bottomT(float3 *pos, float3 *dcos, float3 *normal, float xdetector, float ydetector, int2* seed, float beta, float d_min, float H, float d_max)
	{

		float3 temp_pos, temp_dcos;
		float d_nextCol=0.0f, d_top=0.0f;
		float tmp_deno=0.0f, cos_newangle=0.0f, newangle=0.0f;
		int pexit=0;

		int newnormalctr=0;

		temp_pos.x = pos->x;
		temp_pos.y = pos->y;
		temp_pos.z = pos->z;

		temp_dcos.x = -dcos->x;
		temp_dcos.y = -dcos->y;
		temp_dcos.z = -dcos->z;

		// sample distance uniformly between d_min and d_max to next column
		d_nextCol = ranecuT(seed) * (d_max - d_min) + d_min;

		// calculate distance to bottom surface: if d_bottom < d_nextCol, photon should get detected.
		d_top  = ((H/2.0f) - temp_pos.z)/dcos->z;

		// compute the new position of the photon. 
		pos->x = temp_pos.x + dcos->x * d_nextCol;
		pos->y = temp_pos.y + dcos->y * d_nextCol;
		pos->z = temp_pos.z + dcos->z * d_nextCol;

		// condition to check that pos is within detector boundaries - if true, particle LOST
		if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) )
		{
			atomicAdd(&num_lostT, 1);
			pexit = 1;
			goto myexit;
		}

		if ( (pos->z > H/2.0f)  )		// check if photon's new z position is above top surface
		{
				if( (d_top < d_nextCol) && (d_top > epsilon) )		// photon will hit top surface
				{
					pos->x = temp_pos.x + dcos->x * d_top;
					pos->y = temp_pos.y + dcos->y * d_top;
					pos->z = H/2.0f;
				
					atomicAdd(&photon_distanceT, d_top);
					pexit = 0;			
				}
				else
				{
					atomicAdd(&num_lostT, 1);
					pexit = 1;
					goto myexit;
				}
		}
		else
		{
			atomicAdd(&photon_distanceT, d_nextCol);		// add distance travelled to global variable

			// sample new normal to determine orientation of new column.
		  	mynewnormal:
				normal->x = temp_dcos.x;		// invert dcos of incident vector
				normal->y = temp_dcos.y;
				normal->z = temp_dcos.z;

				RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

				tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

				// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
				normal->z = 0.0f;			// normal_z of a cylinder is always zero
				normal->x = normal->x/tmp_deno;		// re-normalize
				normal->y = normal->y/tmp_deno;

				// perturb the normal according to Beta
				RoughSurfaceT(normal, seed, beta);

				// find the angle between Normal and -Dcos
				cos_newangle = dot_productT(&temp_dcos, normal);
				newangle = acosf(cos_newangle);

				if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
				 {

					if(newnormalctr < 100)	// resample max 100 times
					{
						newnormalctr++;
						goto mynewnormal;			// keep looping until get a theta within 90 degrees	
					}
					else // kill it
					{
						atomicAdd(&num_theta1T,1);
						pexit = 1;
						newnormalctr = 0;
						goto myexit;
					}
				
				 }
	
				pexit = 0;
		}

	myexit:

	 return pexit;
	}	// CUDA refl_bottom function ends
#else
	int refl_bottomT(float *pos, float *dcos, float *normal, float xdetector, float ydetector, int* seed, float beta, float d_min, float H, float d_max)
	{
		float temp_pos[3], temp_dcos[3];
		float d_nextCol=0.0f, d_top=0.0f;
		float tmp_deno=0.0f, cos_newangle=0.0f, newangle=0.0f;
		int pexit=0;
	
		int newnormalctr = 0;

		temp_pos[0] = pos[0];
		temp_pos[1] = pos[1];
		temp_pos[2] = pos[2];

		temp_dcos[0] = -dcos[0];
		temp_dcos[1] = -dcos[1];
		temp_dcos[2] = -dcos[2];

		// sample distance uniformly between d_min and d_max to next column
		d_nextCol = ranecuT(seed) * (d_max - d_min) + d_min;

		// calculate distance to bottom surface: if d_bottom < d_nextCol, photon should get detected.
		d_top  = ((H/2.0f) - temp_pos[2])/dcos[2];

		// compute the new position of the photon. 
		pos[0] = temp_pos[0] + dcos[0] * d_nextCol;
		pos[1] = temp_pos[1] + dcos[1] * d_nextCol;
		pos[2] = temp_pos[2] + dcos[2] * d_nextCol;

		// condition to check that pos is within detector boundaries - if true, particle LOST
		if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) )
		{
			num_lostT++;
			pexit = 1;
			goto myexit;
		}

		if ( (pos[2] > H/2.0f)  )		// check if photon's new z position is above top surface
		{
				if( (d_top < d_nextCol) && (d_top > epsilon) )		// photon will hit top surface
				{
					pos[0] = temp_pos[0] + dcos[0] * d_top;
					pos[1] = temp_pos[1] + dcos[1] * d_top;
					pos[2] = H/2.0f;
				
					photon_distanceT = photon_distanceT + d_top;
					pexit = 0;			
				}
				else
				{
					num_lostT++;
					pexit = 1;
					goto myexit;
				}
		}
		else
		{
			photon_distanceT = photon_distanceT + d_nextCol;		// add distance travelled to global variable

			// sample new normal to determine orientation of new column.
		  	mynewnormal:

				normal[0] = temp_dcos[0];		// invert dcos of incident vector
				normal[1] = temp_dcos[1];
				normal[2] = temp_dcos[2];

				RoughSurfaceT(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-90 degrees of inverted dcos.

				tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

				// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
				normal[2] = 0.0f;			// normal_z of a cylinder is always zero
				normal[0] = normal[0]/tmp_deno;		// re-normalize
				normal[1] = normal[1]/tmp_deno;

				// perturb the normal according to Beta
				RoughSurfaceT(normal, seed, beta);

				// find the angle between Normal and -Dcos
				cos_newangle = dot_productT(temp_dcos, normal);
				newangle = acosf(cos_newangle);

				if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 90 degrees from inverted dcos
				 {

					if(newnormalctr < 100)	// resample max 100 times
					{
						newnormalctr++;
						goto mynewnormal;			// keep looping until get a theta within 90 degrees	
					}
					else // kill it
					{
						num_theta1T++;
						pexit = 1;
						newnormalctr = 0;
						goto myexit;
					}
				
				 }
	
				pexit = 0;
		}

	myexit:

	 return pexit;
	}	// C refl_bottom function ends
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// calculate dot product of two vectors
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline float dot_productT(float3 *aa, float3 *b)
	{
		float result = 0.0f;

		result = aa->x*b->x + aa->y*b->y + aa->z*b->z;

	  return result;
	}
#else
	float dot_productT(float *aa, float *b)
	{
		float result = 0.0f;

		result = aa[0]*b[0] + aa[1]*b[1] + aa[2]*b[2];

	  return result;
	}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// compute directional cosines of transmitted vector
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline void trans_dir_cosT(float3 *dcos, float3 *normal, float refl_theta, float trans_theta, int flag_ref, int mytid, struct start_info *info)
	{
		float cos_angle = 0.0f;
		float norm = 0.0f;
		float3 dcos_temp = {0.0f};

		dcos_temp.x = -dcos->x;
		dcos_temp.y = -dcos->y;
		dcos_temp.z = -dcos->z;
	
		cos_angle = dot_productT(&dcos_temp,normal);	// cosine of angle between incident in opposite direction and normal

		if (flag_ref == 0)				// reflection
		{
				dcos->x = 2.0f*cos_angle*normal->x + dcos->x;  // specular ray
				dcos->y = 2.0f*cos_angle*normal->y + dcos->y;
				dcos->z = 2.0f*cos_angle*normal->z + dcos->z;
		}
		else if (flag_ref == 1)				// transmission	
		{
			 dcos->x= -normal->x*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos->x+(cos_angle*normal->x));
			 dcos->y= -normal->y*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos->y+(cos_angle*normal->y));
			 dcos->z= -normal->z*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos->z+(cos_angle*normal->z));
		}

		norm = sqrt(dcos->x*dcos->x + dcos->y*dcos->y + dcos->z*dcos->z);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
		 {
			dcos->x = dcos->x/norm;
			dcos->y = dcos->y/norm;
			dcos->z = dcos->z/norm;
		 } 

	return;	
	}
#else
	void trans_dir_cosT(float *dcos, float *normal, float refl_theta, float trans_theta, int flag_ref, struct start_info info)
	{
		float cos_angle = 0.0f;
		float norm = 0.0f;
		float dcos_temp[3] = {0.0f};

		dcos_temp[0] = -dcos[0];
		dcos_temp[1] = -dcos[1];
		dcos_temp[2] = -dcos[2];
	
		cos_angle = dot_productT(dcos_temp,normal);	// cosine of angle between incident in opposite direction and normal

		if (flag_ref == 0)				// reflection
		{
				dcos[0] = 2.0f*cos_angle*normal[0] + dcos[0];  // specular ray
				dcos[1] = 2.0f*cos_angle*normal[1] + dcos[1];
				dcos[2] = 2.0f*cos_angle*normal[2] + dcos[2];
		}
		else if (flag_ref == 1)				// transmission	
		{
			 dcos[0]= -normal[0]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[0]+(cos_angle*normal[0]));
			 dcos[1]= -normal[1]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[1]+(cos_angle*normal[1]));
			 dcos[2]= -normal[2]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[2]+(cos_angle*normal[2]));
		}

		norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
		 {
			dcos[0] = dcos[0]/norm;
			dcos[1] = dcos[1]/norm;
			dcos[2] = dcos[2]/norm;
		 } 

	return;	
	}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// add roughness to the surface of the column according to roughness coefficient 'beta'
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline void RoughSurfaceT(float3 *normal, int2* seed, float beta)
	{

		float theta = 0.0f;
		float status = 0.0f;
		float rr = 0.0f;
		float3 normalpert = {0.0f};
		float3 rough_normal = {0.0f};
		float normalize_base = 0.0f;

		// generate the perturbation vector
		status = ranecuT(seed);
		normalpert.z = 2.0f*status - 1.0f;
		rr = sqrt(1.0f - status*status);
		status = ranecuT(seed);
		theta = status * 2.0f * pi;

		normalpert.x = rr * cos(theta);
		normalpert.y = rr * sin(theta);

		// normalize the perturbed vector
		normalize_base = sqrt( pow(normalpert.x,2) + pow(normalpert.y,2) + pow(normalpert.z,2) );
	
		normalpert.x = normalpert.x/normalize_base;
		normalpert.y = normalpert.y/normalize_base;
		normalpert.z = normalpert.z/normalize_base;


		// rough normal = beta*perturbed + original normal
		rough_normal.x = beta * normalpert.x + normal->x;	
		rough_normal.y = beta * normalpert.y + normal->y;
		rough_normal.z = beta * normalpert.z + normal->z;

		// normalize new normal
		normalize_base = sqrt( pow(rough_normal.x,2) + pow(rough_normal.y,2) + pow(rough_normal.z,2) );

		normal->x = rough_normal.x/normalize_base; 
		normal->y = rough_normal.y/normalize_base;
		normal->z = rough_normal.z/normalize_base;


	return;
	}	// CUDA RoughSurface function ends
#else
	void RoughSurfaceT(float *normal, int* seed, float beta)
	{

		float theta = 0.0f;
		float status = 0.0f;
		float rr = 0.0f;
		float normalpert[3] = {0.0f};
		float rough_normal[3] = {0.0f};
		float normalize_base = 0.0f;

		// generate the perturbation vector
		status = ranecuT(seed);
		normalpert[2] = 2.0f*status - 1.0f;
		rr = sqrt(1.0f - status*status);
		status = ranecuT(seed);
		theta = status * 2.0f * pi;

		normalpert[0] = rr * cos(theta);
		normalpert[1] = rr * sin(theta);

		// normalize the perturbed vector
		normalize_base = sqrt( pow(normalpert[0],2) + pow(normalpert[1],2) + pow(normalpert[2],2) );
	
		normalpert[0] = normalpert[0]/normalize_base;
		normalpert[1] = normalpert[1]/normalize_base;
		normalpert[2] = normalpert[2]/normalize_base;


		// rough normal = beta*perturbed + original normal
		rough_normal[0] = beta * normalpert[0] + normal[0];	
		rough_normal[1] = beta * normalpert[1] + normal[1];
		rough_normal[2] = beta * normalpert[2] + normal[2];

		// normalize new normal
		normalize_base = sqrt( pow(rough_normal[0],2) + pow(rough_normal[1],2) + pow(rough_normal[2],2) );

		normal[0] = rough_normal[0]/normalize_base; 
		normal[1] = rough_normal[1]/normalize_base;
		normal[2] = rough_normal[2]/normalize_base;


	return;
	}	// C RoughSurface function ends
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// determine if the photon gets detected
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline int detectionT(float3 *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, 
	size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float sensorRefl, float d_min, int2* seed, float3 *dcos, float3 *normal, float bulk_abscoeff, 
	float R, float xdetector, float ydetector, unsigned long long int mynum_rebound, float *XcT, float *YcT)
	{
		int result = 0;
		int ii = 0, jj = 0;
		int absflag = 0;

		// equation of plane is z = -H/2
		// if a point satisfies above equation, it is detected

		if (fabs(pos->z - (float)(-H/2.0f)) < epsilon) 
		 {
			// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
			if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
			{

				// normal pointing (0,0,1)
				normal->x = 0.0f; 
				normal->y = 0.0f; 
				normal->z = 1.0f;

				// obtain reflected dcos from the bottom (specular reflection; bottom surface is smooth, no need to perturb the normal)
				// this function is called only when photon hits bottom of a column
				trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, mytid, info); // as their is only reflection so refl_theta,trans_theta = 0

				// using above calculated dcos, move the particle within the column
				absflag = isotropicT(pos, dcos, seed, bulk_abscoeff, R, H, xdetector, ydetector, &info[mytid], mynum_rebound, XcT, YcT, mytid);

				if(absflag == 1)	// if gets absorbed in bulk, photon exits
					result = 1;
				else
					result = 0;

			}
			else	// photon absorbed or detected  
			{
				result = 1;
				atomicAdd(&num_detectT, 1);
		
				ii = floor((pos->x-lbound_x)/pixelsize);	// determine pixel number in x and y direction
				jj = floor((pos->y-lbound_y)/pixelsize);

				// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
				if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
				 {	
					unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + ii * pitch);
					atomicAdd(&current_img[jj],1);
				 }

				atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);	// start array from 0.str_histnum starts from 1
			}
			 
		 }
		else
		    	result = 0;
	  
	//exit1:
	 return result;
	}	// CUDA detection function ends
#else
	int detectionT(float *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, 		float sensorRefl, float d_min, int* seed, float *dcos, float *normal, float bulk_abscoeff, float R, float xdetector, float ydetector, unsigned long long int mynum_rebound, 
	int ydim, int *h_num_detected_prim)
	{

		int result = 0;
		int ii = 0, jj = 0;
		int absflag = 0;


		// equation of plane is z = -H/2
		// if a point satisfies above equation, it is detected

		if (fabs(pos[2] - (float)(-H/2.0f)) < epsilon) 
		 {
			// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
			if(ranecuT(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
			{

				// normal pointing (0,0,1)
				normal[0] = 0.0f; 
				normal[1] = 0.0f; 
				normal[2] = 1.0f;

				// obtain reflected dcos from the bottom (specular reflection; bottom surface is smooth, so no need to perturb the normal)
				// this function is called only when photon hits bottom of a column
				trans_dir_cosT(dcos, normal, 0.0f, 0.0f, 0, info); // as their is only reflection so made refl_theta, trans_theta = 0

				// using above calculated dcos, move the particle within the column
				absflag = isotropicT(pos, dcos, seed, bulk_abscoeff, R, H, xdetector, ydetector, info, mynum_rebound);

				if(absflag == 1)	// if gets absorbed in bulk, photon exits
					result = 1;
				else
					result = 0;

			}
			else	// photon absorbed or detected  
			{
				result = 1;
				num_detectT++;
		
				ii = floor((pos[0]-lbound_x)/pixelsize);	// determine pixel number in x and y direction
				jj = floor((pos[1]-lbound_y)/pixelsize);

				// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
				if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
				 {	
					outputimage_.tempimageopt[ii][jj]++;
				 }

			}
			 
		 }
		else
		    	result = 0;
	  
	 return result;
	}
#endif


////////////////////////////////////////////////////////////////////////////////
//! Initialize the pseudo-random number generator (PRNG) RANECU to a position
//! far away from the previous history (leap frog technique).
//!
//! Each calculated seed initiates a consecutive and disjoint sequence of
//! pseudo-random numbers with length LEAP_DISTANCE, that can be used to
//! in a parallel simulation (Sequence Splitting parallelization method).
//! The basic equation behind the algorithm is:
//!    S(i+j) = (a**j * S(i)) MOD m = [(a**j MOD m)*S(i)] MOD m  ,
//! which is described in:
//!   P L'Ecuyer, Commun. ACM 31 (1988) p.742
//!
//! This function has been adapted from "seedsMLCG.f", see:
//!   A Badal and J Sempau, Computer Physics Communications 175 (2006) p. 440-450
//!
//!       @param[in] history   Particle bach number.
//!       @param[in] seed_input   Initial PRNG seed input (used to initiate both MLCGs in RANECU).
//!       @param[out] seed   Initial PRNG seeds for the present history.
//!
////////////////////////////////////////////////////////////////////////////////
// -- Upper limit of the number of random values sampled in a single track:
#define  LEAP_DISTANCE    1000
// -- Multipliers and moduli for the two MLCG in RANECU:
#define  a1_RANECU       40014
#define  m1_RANECU  2147483563
#define  a2_RANECU       40692
#define  m2_RANECU  2147483399

#ifdef USING_CUDA
	__device__ inline void init_PRNGT(int history_batch, int histories_per_thread, int seed_input, int2* seed)
	{
	  // -- Move the RANECU generator to a unique position for the current batch of histories:
	  //    I have to use an "unsigned long long int" value to represent all the simulated histories in all previous batches
	  //    The maximum unsigned long long int value is ~1.8e19: if history >1.8e16 and LEAP_DISTANCE==1000, 'leap' will overflow.
	  // **** 1st MLCG:
	  unsigned long long int leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  int y = 1;
	  int z = a1_RANECU;
	  // -- Calculate the modulo power '(a^leap)MOD(m)' using a divide-and-conquer algorithm adapted to modulo arithmetic
	  for(;;)
	  {
	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  Equivalent: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2);  
	      y = abMODmT(m1_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODmT(m1_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm1 = y;                 // Exponentiation finished:  AjMODm = expMOD = y = a^j

	  // -- Compute and display the seeds S(i+j), from the present seed S(i), using the previously calculated value of (a^j)MOD(m):
	  //         S(i+j) = [(a**j MOD m)*S(i)] MOD m
	  //         S_i = abMODmT(m,S_i,AjMODm)
	  seed->x = abMODmT(m1_RANECU, seed_input, y);     // Using the input seed as the starting seed

	  // **** 2nd MLCG (repeating the previous calculation for the 2nd MLCG parameters):
	  leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  y = 1;
	  z = a2_RANECU;
	  for(;;)
	  {
	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  Equivalent: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	      y = abMODmT(m2_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODmT(m2_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm2 = y;
	  seed->y = abMODmT(m2_RANECU, seed_input, y);     // Using the input seed as the starting seed
	}
#else
	void init_PRNGT(int history_batch, int histories_per_thread, int seed_input, int* seed)
	{
	  // -- Move the RANECU generator to a unique position for the current batch of histories:
	  //    I have to use an "unsigned long long int" value to represent all the simulated histories in all previous batches
	  //    The maximum unsigned long long int value is ~1.8e19: if history >1.8e16 and LEAP_DISTANCE==1000, 'leap' will overflow.
	  // **** 1st MLCG:
	  unsigned long long int leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  int y = 1;
	  int z = a1_RANECU;
	  // -- Calculate the modulo power '(a^leap)MOD(m)' using a divide-and-conquer algorithm adapted to modulo arithmetic
	  for(;;)
	  {

	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  !!DeBuG!! OLD: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2); 
	      y = abMODmT(m1_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODmT(m1_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm1 = y;                 // Exponentiation finished:  AjMODm = expMOD = y = a^j

	  // -- Compute and display the seeds S(i+j), from the present seed S(i), using the previously calculated value of (a^j)MOD(m):
	  //         S(i+j) = [(a**j MOD m)*S(i)] MOD m
	  //         S_i = abMODmT(m,S_i,AjMODm)
	  seed[0] = abMODmT(m1_RANECU, seed_input, y);     // Using the input seed as the starting seed

	  // **** 2nd MLCG (repeating the previous calculation for the 2nd MLCG parameters):
	  leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  y = 1;
	  z = a2_RANECU;
	  for(;;)
	  {
	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  !!DeBuG!! OLD: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2); 
	      y = abMODmT(m2_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODmT(m2_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm2 = y;
	  seed[1] = abMODmT(m2_RANECU, seed_input, y);     // Using the input seed as the starting seed

	}
#endif


/////////////////////////////////////////////////////////////////////
//!  Calculate "(a1*a2) MOD m" with 32-bit integers and avoiding   **
//!  the possible overflow, using the Russian Peasant approach     **
//!  modulo m and the approximate factoring method, as described   **
//!  in:  L'Ecuyer and Cote, ACM Trans. Math. Soft. 17 (1991)      **
//!                                                                **
//!  This function has been adapted from "seedsMLCG.f", see:       **
//!  Badal and Sempau, Computer Physics Communications 175 (2006)  **
//!                                                                **
//!    Input:          0 < a1 < m                                  **
//!                    0 < a2 < m                                  **
//!                                                                **
//!    Return value:  (a1*a2) MOD m                                **
//!                                                                **
/////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline int abMODmT(int m, int a, int s)
	{
	  // CAUTION: the input parameters are modified in the function but should not be returned to the calling function! (pass by value!)
	  int q, k;
	  int p = -m;            // p is always negative to avoid overflow when adding

	  // ** Apply the Russian peasant method until "a =< 32768":
	  while (a>32768)        // We assume '32' bit integers (4 bytes): 2^(('32'-2)/2) = 32768
	  {
	    if (0!=(a&1))        // Store 's' when 'a' is odd     Equivalent code:   if (1==(a%2))
	    {
	      p += s;
	      if (p>0) p -= m;
	    }
	    a >>= 1;             // Half a (move bits 1 position right)   Equivalent code: a = a/2;
	    s = (s-m) + s;       // Double s (MOD m)
	    if (s<0) s += m;     // (s is always positive)
	  }

	  // ** Employ the approximate factoring method (a is small enough to avoid overflow):
	  q = (int) m / a;
	  k = (int) s / q;
	  s = a*(s-k*q)-k*(m-q*a);
	  while (s<0)
	    s += m;

	  // ** Compute the final result:
	  p += s;
	  if (p<0) p += m;

	  return p;
	}
#else
	int abMODmT(int m_par, int a_par, int s_par)
	{
	  // CAUTION: the input parameters are modified in the function but should not be returned to the calling function! (pass by value!)   !!DeBuG!!
	  int mval,aval,sval;
	  mval=m_par; aval=a_par; sval=s_par;
	  
	  int qval, kval;
	  int pval = -mval;            // p is always negative to avoid overflow when adding

	  // ** Apply the Russian peasant method until "a =< 32768":
	  while (aval>32768)        // We assume '32' bit integers (4 bytes): 2^(('32'-2)/2) = 32768
	  {
	    if (0!=(aval&1))        // Store 's' when 'a' is odd    !!DeBuG!! OLD code:   if (1==(a%2))
	    {
	      pval += sval;
	      if (pval>0) pval -= mval;
	    }
	    aval >>= 1;             // Half a (move bits 1 position right)        
	    sval = (sval-mval) + sval;       // float s (MOD m)
	    if (sval<0) sval += mval;     // (s is always positive)
	  }

	  // ** Employ the approximate factoring method (a is small enough to avoid overflow):
	  qval = (int) mval / aval;
	  kval = (int) sval / qval;
	  sval = aval*(sval-kval*qval)-kval*(mval-qval*aval);
	  while (sval<0)
	    sval += mval;

	  // ** Compute the final result:
	  pval += sval;
	  if (pval<0) pval += mval;

	  return pval;
	}
#endif


////////////////////////////////////////////////////////////////////////////////
//! Pseudo-random number generator (PRNG) RANECU returning a float value
//! (single precision version).
//!
//!       @param[in,out] seed   PRNG seed (seed kept in the calling function and updated here).
//!       @return   PRN double value in the open interval (0,1)
//!
////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline float ranecuT(int2* seed)
	{
	//return (float(seed->x%100)*0.01f+0.005f)  ;

	  int i1 = (int)(seed->x/53668);
	  seed->x = 40014*(seed->x-i1*53668)-i1*12211;

	  int i2 = (int)(seed->y/52774);
	  seed->y = 40692*(seed->y-i2*52774)-i2*3791;

	  if (seed->x < 0) seed->x += 2147483563;
	  if (seed->y < 0) seed->y += 2147483399;

	  i2 = seed->x-seed->y;
	  if (i2 < 1) i2 += 2147483562;

	  return (__int2float_rn(i2)*4.65661305739e-10f);        // 4.65661305739e-10 == 1/2147483563

	}
#else
	float ranecuT(int* seed)
	{
	  int i1 = (int)(seed[0]/53668);
	  seed[0] = 40014*(seed[0]-i1*53668)-i1*12211;

	  int i2 = (int)(seed[1]/52774);
	  seed[1] = 40692*(seed[1]-i2*52774)-i2*3791;

	  if (seed[0] < 0) seed[0] += 2147483563;
	  if (seed[1] < 0) seed[1] += 2147483399;

	  i2 = seed[0]-seed[1];
	  if (i2 < 1) i2 += 2147483562;

	  const float USCALE = 1.0/2147483563.0;       
	  return ((float)(i2*USCALE));

	}
#endif

