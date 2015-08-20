#!/bin/bash

# CUDA codes
nvcc -O3 -use_fast_math -DUSING_CUDA -c hybridMANTIS_cuda_ver1_0_LB.cu -I/usr/local/cuda/include -I/opt/cuda_sdk_4.0/C/common/inc/ -L/opt/cuda_sdk_4.0/C/lib/ -L/usr/local/cuda/lib64/ -lcutil -lcudart -lgsl -lgslcblas -lm -arch sm_20 --ptxas-options=-v

nvcc -O3 -use_fast_math -DUSING_CUDA -c hybridMANTIS_cuda_ver1_0.cu -I/usr/local/cuda/include -I/opt/cuda_sdk_4.0/C/common/inc/ -L/opt/cuda_sdk_4.0/C/lib/ -L/usr/local/cuda/lib64/ -lcutil -lcudart -lgsl -lgslcblas -lm -arch sm_20 --ptxas-options=-v

# C codes
gcc -I/usr/include -c hybridMANTIS_c_ver1_0_LB.c -lgsl -lgslcblas -lm

gcc -I/usr/include -c hybridMANTIS_c_ver1_0.c -lgsl -lgslcblas -lm

# link with PENELOPE
gfortran -O3 -g hybridMANTIS_penEasy.f hybridMANTIS_tallyEnergyDepositionEvents.f penelope.f pengeom.f penvared.f hybridMANTIS_cuda_ver1_0_LB.o hybridMANTIS_cuda_ver1_0.o hybridMANTIS_c_ver1_0_LB.o hybridMANTIS_c_ver1_0.o -o hybridMANTIS_ver1_0.x -I/usr/local/cuda/include -I/opt/cuda_sdk_4.0/C/common/inc/  -L/opt/cuda_sdk_4.0/C/lib/ -L/usr/local/cuda/lib64/ -lcutil -lcudart -lgsl -lgslcblas -lm
