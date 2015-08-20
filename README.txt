##########################################################################################################################################################
#
# ****Disclaimer****
#  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
#  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
#  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
#  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
#  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
#  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
#  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
#  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are
#
#	@file    README.txt
#       @author  Diksha Sharma (Diksha.Sharma@fda.hhs.gov)
#       @date    Apr 19, 2012
#
##########################################################################################################################################################


*************
INTRODUCTION
*************

hybridMANTIS is a Monte Carlo tool for modeling x-ray detectors with columnar scintillators. It uses a novel hybrid approach to maximize the utilization of computing resources like CPUs and GPUs in modern workstations. hybridMANTIS was developed using the CUDA programming model from NVIDIA to achieve maximum performance on NVIDIA GPUs. The code can also be compiled with a standard C compiler to be executed in a regular CPU.

hybridMANTIS uses PENELOPE 2006 for the x-ray and electron transport with the penEasy main program and tallies, and fastDETECT2 for the optical transport. fastDETECT2 is a new and improved version of the DETECT2 (optical transport used in MANTIS). The penEasy program was modified to output the energy and location of individual energy deposition events which are then used as input to fastDETECT2 which then spawns independent optical transport kernel calls.

The source code is free and open software in the public domain, as explained in the Disclaimer section below. 
The software distribution website is: http://code.google.com/p/hybridmantis/. 


*******************************
CODE COMPILATION AND EXECUTION
*******************************

hybridMANTIS has been tested only on the Linux operating system. The CUDA libraries, GNU gcc and gfortran compiler and GNU scientific library needs to be pre installed before running hybridMANTIS. A bash script is included to compile the CUDA and C codes with PENELOPE and penEasy files. This file may have to be edited to modify the library paths. A pre-compiled executable is also attached. It was compiled using CUDA version 4.0, gcc version 4.4.5 and gsl version 1.14.

To run hybridMANTIS:
	hybridMANTIS ver1 0.x < penEasy CsI input.in
	
Sample inputs and outputs for simulating 100,000 x-ray histories are included under the 'example' folder. For more details, read 'MANUAL_hybridMANTIS.pdf'.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																//
//    NOTE: PENELOPE 2006 SOURCE CODE FILES PENELOPE.F AND PENGEOM.F ARE NOT DISTRIBUTED WITH hybridMANTIS PACKAGE, 		//
//    BUT ARE NEEDED FOR COMPILING. IF THESE FILES ARE NEEDED, PLEASE CONTACT DIKSHA SHARMA AT Diksha.Sharma(at)fda.hhs.gov.	//
//																//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


**********************************
hybridMANTIS v1.0 PACKAGE CONTENTS
**********************************

	* fastDETECT2 SOURCE CODES: 	'hybridMANTIS_****' and 'kernel_****' CUDA and C files. 
	
	* penEasy SOURCE CODES: 	All fortran files. Files 'hybridMANTIS_peneasy.f' and 'hybridMANTIS_tallyEnergyDepositionEvents.f' were derived from original penEasy files. These files have been modified to output energy deposition events information and to include the load balancing algorithm to be able to run with fastDETECT2 codes via hybrid technique.
	
	* penEASY IMAGING:	 	This folder contains penEasy documentation, copyright notice and gnuplot scripts.
	
	* example:			This folder contains sample input and output files for simulating 100,000 x-ray histories using hybridMANTIS.
	
	* GNUPLOT_scripts:		hybridMANTIS gnuplot scripts to plot pulse-height spectrum and point response function using output files 'detected_', 'myimage_'.
	
	* hybridMANTIS_input.in:	hybridMANTIS input file required for running hybridMANTIS.
	
	* penEasy_CsI_input.in:		penEasy input file required for running hybridMANTIS.
	
	* MANUAL_hybridMANTIS.pdf: 	Reference manual for hybridMANTIS v1.0.
	
	* hybridMANTIS_ver1_0.x:	hybridMANTIS v1.0 executable.
	
	* compile_ver1_0.sh:		Script for compiling hybridMANTIS. The user will require two more files PENELOPE.F, PENGEOM.F to compile this package successfully. These files are not distributed with this package. If these files are needed, please contact Diksha Sharma at Diksha.Sharma (at) fda.hhs.gov.
	
	* README.txt:			this file.
	

For more details, read 'MANUAL_hybridMANTIS.pdf'.
	
