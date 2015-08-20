# hybridMANTIS
## Project Summary

hybridMANTIS is a Monte Carlo package for modeling indirect x-ray detectors with columnar scintillators using a novel hybrid CPU-GPU technique. This hybrid approach maximizes the utilization of CPUs and GPUs in modern workstations.

hybridMANTIS uses PENELOPE for the x-ray/electron transport and fastDETECT2 for the optical transport. fastDETECT2 is a complete rewrite and improved version of DETECT2 (optical transport code used in MANTIS). It includes several new features like on-the-fly column geometry and columnar crosstalk to model the columnar arrays more realistically as compared to MANTIS. A load balancer is implemented to dynamically allocate optical transport showers to the GPU and CPU computing cores.

The load balancing algorithm and the use of GPUs in hybrid with the CPU makes hybridMANTIS significantly faster than MANTIS. In MANTIS, the optical transport takes most of the time as compared to the x-ray and electron transport. Using hybridMANTIS, we were able to successfully hide hours of optical transport time by running it in parallel with the x-ray and electron transport, thus shifting the computational bottleneck from optical to x-ray transport. The new code requires much less memory than MANTIS and, as a result, allows us to efficiently simulate clinical-size, large-area imaging detectors.

hybridMANTIS is being developed at the U.S. Food and Drug Administration, Center for Devices and Radiological Health, Office of Science and Engineering Laboratories, Division of Imaging and Applied Mathematics. As explained in the Disclaimer below, this package is under the public domain and is free to be downloaded and distributed. This code is still under development, please report to the authors any issue/bug you may encounter.

hybridMANTIS has been mentioned in several publications, but the below paper should be used as a reference by researchers using this code.

    Diksha Sharma, Andreu Badal and Aldo Badano, "hybridMANTIS: a CPU-GPU Monte Carlo method for modeling indirect x-ray detectors with columnar scintillators", Physics in Medicine and Biology, 57 (8), p. 2357-72 (2012). 


## Comparison with experimental and MANTIS results

hybridMANTIS data for modeling a CsI detector was compared with the experimental and MANTIS results (obtained from Freed et al., Medical Physics, 36(11), 4944–56, 2009). Metrics like point response and modulation transfer functions were used for comparison. The results were analyzed quantitatively based on the Swank factor and root mean square values. Our results suggests that hybridMANTIS matches the experimental data as good as or better than MANTIS, being significantly more computationally efficient. This work was presented at the 11th International Workshop on Breast Imaging (IWDM 2012). 

## Disclaimer

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 
