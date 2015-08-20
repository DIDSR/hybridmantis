                      penEasy README FILE
                      version 2008-06-15
                       by Josep Sempau
                Universitat Politecnica de Catalunya
                  e-mail: josep.sempau@upc.es

CONTENTS:

0) COPYRIGHT AND DISCLAIMER
1) PURPOSE OF THIS SOFTWARE
2) WHAT YOU GET
3) HOW TO USE IT
4) EXAMPLE
5) VOXELIZED GEOMETRIES
6) PARALLEL EXECUTION
7) KNOWN ISSUES AND LIMITATIONS
8) WHERE DO I GET THE LATEST UPDATE?
9) WHO DO I COMPLAIN TO?



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
0) COPYRIGHT AND DISCLAIMER

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c  penEasy                                                                     c
c  Copyright (c) 2004-2008                                                     c
c  Universitat Politecnica de Catalunya                                        c
c                                                                              c
c  Permission to use, copy, modify and re-distribute copies of this software   c
c  or parts of it and its documentation for any purpose is hereby granted      c
c  without fee, provided that this copyright notice appears in all copies.     c
c  The Universitat Politecnica de Catalunya makes no representations about     c
c  the suitability of this software for any purpose. It is provided "as is"    c
c  without express or implied warranty.                                        c
c                                                                              c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
1) PURPOSE OF THIS SOFTWARE

PenEasy is a structured general-purpose main program for PENELOPE. PENELOPE
is a Monte Carlo simulation package capable of transporting photons,
electrons and positrons in complex geometries and arbitrary materials. It
is freely distributed by the OECD Nuclear Energy Agency Data Bank
(http://www.nea.fr). Detailed information about PENELOPE can be found in
its accompanying documentation (see [1]). Hereafter it is assumed that the
reader is familiar with the basic concepts of Monte Carlo simulation of
radiation transport and with the operation of PENELOPE.

Although most of the intricate features of the transport process are taken
care of by PENELOPE's subroutines, users of this simulation system are
required to write a steering main program which, among other things, should
define the initial particle states (i.e. the radiation source) and the
tallies. PenEasy is intended to accomplish two goals. First, it provides a
set of source models and tally schemes that are directly applicable to a
large variety of practical situations. Thus, users are only asked to
provide an input data file and no programming is required. In some cases
the problem in hand does not fit into one of the pre-defined models. The
second goal is then to furnish the programmer with a structured and modular
code that facilitates the adaptation of the existing routines and the
creation of new ones, reducing the programming effort to a minimum.

PenEasy, like PENELOPE, is free and open software mostly written in FORTRAN
77, although it has recourse to some features included in the Fortran 95
standard. In particular, the time routines comply with this standard.

The current version of penEasy is compatible with PENELOPE 2006.


Reference:

[1] F. Salvat, J.M. Fernández-Varea and J. Sempau, PENELOPE-2006: A Code
  System for Monte Carlo Simulation of Electron and Photon Transport (3rd
  edition), OECD-NEA 2006, Issy-les-Moulineaux, France. Freely available
  from http://www.nea.fr/html/science/pubs/2006/nea6222-penelope.pdf



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
2) WHAT YOU GET

The penEasy package consists of an all-purpose main program, called
penEasy.f, and a set of subroutine libraries. After decompressing the
distribution ZIP file your working directory (represented below as ~/)
should contain the following files and subdirectories (in alphabetical
order):

- ~/documentation/changeHistory.txt
  List of changes and bugs corrected with respect to previous versions.

- ~/documentation/sourceXXX.txt
  ~/documentation/tallyXXX.txt
  Text files containing a description of the input data for each source and
  tally.

- ~/fortranCode/penaux.f
  Auxiliary routines that initialize penEasy/PENELOPE.

- ~/fortranCode/penEasy.f
  The main program source code.

- ~/fortranCode/penpatch.f
  Subroutines that supersede their counterparts in the PENELOPE package.
  The physics in these routines is identical to the original, but provide
  more information on the outcome of some interaction processes.

- ~/fortranCode/penvox.f
  Geometry package for the simulation of voxelized geometries.

- ~/fortranCode/sourceXXX.f
  Each source model is coded in a different file, e.g.
  sourcePhaseSpaceFile.f. They are described in the accompanying
  ~/documentation/sourceXXX.txt files.

- ~/fortranCode/tallyXXX.f
  Each tally is coded in a different file, e.g. tallySpatialDoseDistrib.f.
  They are described in the accompanying ~/documentation/tallyXXX.txt
  files.

- ~/fortranCode/timing.f
  Time routines compatible with the Fortran 95 and subsequent standards.

- ~/gnuplotScripts/*.gpl
  Gnuplot (4.2) scripts to represent graphically the simulation results.
  Each script has a name associated to a tally. See the example provided
  below to learn how to execute them and refer to the gnuplot manual for
  details on how to modify these scripts if needed.

- ~/run/command.in
  Users may change some simulation settings by editing this file while the
  program is running--e.g. to stop the execution.

- ~/run/penEasy.exe
  The executable code for Windows systems obtained with Compaq Visual
  Fortran 6.1. Note that in previous releases the main program was called
  penmain.

- ~/run/penEasy.in
  Sample input file for penEasy. It also serves to run the example case
  described below.

- ~/run/phantom.geo
  A simple quadrics (PENGEOM model) geometry file that is used in the
  example case described below. It also contains, for your convenience, the
  general layout of a quadric geometry file (taken from the PENELOPE
  distribution).

- ~/run/water.mat
  This is the PENELOPE material data file for water. It is included for
  your convenience, so that the example in this README (see below) can be
  run without further preparations.

- ~/voxGeoSample/sample.vox
  A simple voxels geometry file that complies with the penVox syntax. It
  contains a detailed description of the data format.

- ~/README.txt
  This file.

All files in the package (except the Windows executable) are in plain text
(i.e. ASCII) format and, except for this readme, use the Unix new line
convention, which differs from that of Windows. As a result, some text
editors that run on Windows, notably Microsoft Notepad, may not interpret
correctly the new line mark and display scrambled text. This problem can be
solved by using a better editor such as the Programmer's File Editor (PFE)
which can read and save files in both Unix and Windows formats. PFE is
freely available from the Internet. Alternatively, the MS-DOS Editor
available in most Windows systems (execute "edit") can also read both
formats. Most compilers that run on Windows accept both new line
conventions quietly.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
3) HOW TO USE IT

NOTE: In this section it is assumed that the installation has been carried
      out on a GNU/Linux system. For concreteness, the g95 Fortran compiler
      (freely available from http://www.g95.org) is also assumed. For
      Windows systems, an executable binary file is included in the
      distribution (see above) and, therefore, step c) is unnecessary.

The "installation" is accomplished following these steps:

a) Copy all the penEasy files into a directory on your disk. We shall
   assume it has been named ~/penEasy/.

b) Copy the files penelope.f, pengeom.f and penvared.f from the PENELOPE
   package into ~/penEasy/fortranCode/.

c) Compile and link penEasy.f, as in
   $ g95 penEasy.f -o penEasy.x -O

   This generates the executable file (penEasy.x), which should be moved to
   ~/penEasy/run/. This completes the "installation" procedure.

   Note that, although you must not include all the Fortran source files in
   the compilation command, those files should be present in the same
   directory where penEasy is located for the compilation to be performed
   successfully.

d) In order to be able to visualize the results in graphic format, the
   penEasy package includes a set of gnuplot scripts. (gnuplot is a
   function and data plotting software freely available from
   http://www.gnuplot.info ). If you plan to use these scripts, you will
   also need to have gnuplot installed on your machine.


To perform a simulation, take these steps:

a') It is advisable to create a new directory in your working area where
    all the files for your job will be stored. We shall suppose that this
    directory is named ~/mySimul/. Always keep the original, unaltered
    files in ~/penEasy/.

b') Copy the following files to ~/mySimul/.
        ~/penEasy/fortranCode/penEasy.x
        ~/penEasy/run/command.in
        ~/penEasy/run/penEasy.in
        ~/penEasy/gnuplotScripts/*
    The geometry and material files required by PENELOPE need also be in
    ~/mySimul/.

c') Edit the file penEasy.in, where the input data for penEasy is
    introduced. It is structured in sections, each one starting with a
    string of the type '[SECTION ...]'. Each source or tally has its own
    section. Refer to the file itself for instructions on the configuration
    of the various sections.

d') To run penEasy, execute this

    $ penEasy.x < penEasy.in > penEasy.out

    The main program reads from the standard input (the keyboard) and
    writes to the standard output (the screen), so it must be run
    redirecting these devices to the appropiate external files with the
    symbols '<' and '>', as shown in the previous command line.

e') The course of the simulation can be controlled by sending commands to
    the program while it is running. This is accomplished by modifying the
    file command.in. The available commands and their coding are briefly
    explained inside the file itself. For instance, the simulation can be
    terminated and the final results printed by setting the number of
    requested histories to zero at any time. The file command.in is scanned
    periodically, as determined by the update interval defined in
    penEasy.in in terms of either the simulation time or the number of
    histories.

f') Simulation results are written in a separate file for each tally. These
    filed are named after the corresponding tally and have extension *.dat.
    The gnuplot scripts included in the distribution use these DAT files to
    display simulation results graphically. A summary of the execution is
    also printed to the output file (penEasy.out in the example above).

    Note that the DAT files are regularly written at each update interval
    so that users can track the progress of the simulation.

    All uncertainties reported by penEasy, and indicated by '+-', are two
    standard deviations (2 sigma).


Note for **advanced users** on the use of TALLY:

This technical note is useful in case you plan to modify an existing tally
or wish to create your own. Notice that TALLY, which takes the arguments
MODE (integer) and ARG (real*8), is called at various points of the penEasy
main program where it is expected to perform different actions. The
possible values of MODE identify the point at which TALLY is called so that
the appropiate action, if any, can be programmed. The following table
specifies the situation that corresponds to each mode and the value passed
as ARG in each case.

In this table 'E' represents the particle's kinetic energy. By a "history"
we mean the simulation of the primary particle (or particles) produced by a
single call to the SOURCE routine and the simulation of all the secondaries
generated by it (or them).


MODE  situation/required action(if any)       ARG
----  --------------------------------------  ----------------------------
-99   A new particle (either primary or       -E, the kinetic energy of the
      secondary) has been retrieved from      retrieved particle with reversed
      the stack and its simulation is about   sign; ARG should be *ADDED* to
      to begin. Its kinetic energy should     the deposited energy counters.
      be subtracted from the deposited
      energy counters.

-98   The particle has been absorbed. Its     E. Usually, KNOCK sets E=0 after
      remaining kinetic energy should be      absorption, but a particle can also
      added to all deposited energy           be absorbed BEFORE any interaction
      counters.                               when it enters a medium with an Eabs
                                              larger than E.

-97   A positron has been absorbed with       2mc^2 = 1022 keV aprox.
      E>0 (i.e. not as a result of a
      KNOCK but because it entered a region
      with lower Eabs) and two photons have
      been pushed to the stack.

from  KNOCK informed that the particle had    DE, the energy lost.
  -1  an interaction of type ICOL=-MODE (see
  to  PENELOPE manual for a description of
  -8  what interaction corresponds to each
      value of ICOL returned by KNOCK). The
      energy DE lost by the particle should
      be added to all deposited energy
      counters.

  0   Routine SOURCE calls TALLY with         E.
      MODE=0 each time it adds a new
      born particle to the stack. Its
      kinetic energy E should be added
      to all deposited energy counters.

  =>  Notice that in all MODE.le.0 cases the value of ARG should be added
      to the deposited energy counters.

  1   A new history is about to begin.        N, the current history number.

  2   The history number N has been changed.  N, the actual history number.
      This is used by the PSF source
      to adjust N to the value dictated by
      the contents of the PSF.

  3   The particle has been moved a distance  DS, the distance travelled.
      DS by STEP without any interface
      crossing and it is about to interact.
      This information is useful e.g. for
      track length estimators (to obtain
      the fluence).

  4   The particle has been moved a distance  DSEF, the distance travelled.
      DSEF by STEP and an interface crossing
      has ocurred at the end of the track.
      Again, useful e.g. for track length
      estimators.

  6   A history has been finished and it is   dble(N), the history number
      necessary to perform some bookkeeping   just finished
      to store the average and variance of
      the quantities being calculated.


The setting MODE=5, present in previous versions of penEasy, is no longer
in use.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
4) EXAMPLE

NOTE: In this section it is assumed that the installation has been
      performed on a GNU/Linux system.

Consider a 10 MeV electron beam impinging on a semi-infinite (z>0) water
phantom. The beam is modelled as emerging from a point source located at
the Cartesian coordinates (0,0,-100), i.e., at 100 cm from the water
surface. Electrons are emitted with an initial direction limited to a cone
with its axis along the vector (0,0,1) and with an angular semiaperture
alpha of 2.86241 deg. Since tg(alpha)=0.05, the field on the water surface
is a circle with radius equal to 100*tg(alpha)=5 cm.

We wish to calculate the energy deposited per unit depth interval, that is,
a depth-dose curve. More precisely, penEasy will report dE/(rho*dz) in
keV*cm^2/g, where dE is the energy deposited in a bin, rho is the mass
density of water and dz is the bin width along the z direction. The dose
distribution in cylindrical coordinates (as a function of z and the radial
distance from the beam axis) will also be computed, although the radiation
field does not have an exact cylindrical symmetry due to its divergence,
determined by alpha.

The depth interval from z=0 up to z=7 cm will be partitioned in 40 bins.
For the cylindrical distribution, the radial interval [0,8] cm will be
divided into 80 bins. Results will be reported per history, that is, per
unit emitted electron.

Electrons, positrons and photons will be followed down to 100, 100 and 10
keV, respectively. The cutoffs for the production of secondary electrons
and bremsstrahlung photons will be set to 100 and 10 keV, respectively. The
simulation will be stopped after 100,000 histories, or when the average
uncertainty in the depth-dose is 1% (2 sigma), whatever comes first. No
interaction forcing (a variance reduction technique available from
PENELOPE) will be applied.

Now, follow these steps:

a) Create the directory ~/mySimul/ in your working area.

b) Copy these files to ~/mySimul/ :
        ~/penEasy/run/command.in
        ~/penEasy/run/penEasy.x
             (See installation procedure above; on a Windows system
              you may use the ready-made penEasy.exe provided with
              the distribution.)
        ~/penEasy/run/penEasy.in  (input file prepared for this example)
        ~/penEasy/run/phantom.geo (semiinfinite phantom, PENGEOM syntax)
        ~/penEasy/run/water.mat   (PENELOPE material file)
        ~/gnuplotScripts/*

   As in any other PENELOPE simulation, a material data file, water in our
   case, generated with the program MATERIAL (included in PENELOPE) is
   required. For convenience, this file is provided with penEasy for this
   example.

c) You may want to edit penEasy.in and check that the information
   introduced in this file reflects the definition of our depth-dose
   problem. Except for the tally ParticleTrackStructure, tallies other than
   those producing the requested distributions have been turned OFF to save
   CPU time.

   A detailed description of the simulation parameters EABS,C1,C2,WCC,WCR
   and DSMAX can be found in the PENELOPE documentation.

d) Run the program

   $ penEasy.x < penEasy.in > penEasy.out &
   (on Unix systems, the final '&' puts the process in background)

   When the simulation ends--it takes about one minute--the sought depth-
   dose data can be found in tallySpatialDoseDistrib.dat. This file is
   updated file every 50 s ('UPDATE' parameter in penEasy.in) whilst the
   program is running. Recall that the quoted uncertainties are at the 2
   sigma level.

e) In case you installed gnuplot, you can visualize the depth dose curve by
   executing the corresponding script from the ~/mySimul/ directory:

   $ gnuplot tallySpatialDoseDistrib-1D.gpl

   To visualize the cylindrical dose distribution use:

   $ gnuplot tallyCylindricalDoseDistrib-rz.gpl

   You can also visualize a few particles tracks executing:

   $ gnuplot tallyParticleTrackStructure.gpl

f) You may send commands to penEasy while it is running by editing the file
   command.in. Notice that, in the present example, command.in is scanned
   every 50 s (UPDATE INTERVAL).



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
5) VOXELIZED GOMETRIES

Three possible geometry models can be employed in penEasy: (i) quadrics;
(ii) voxels, that is, homogeneous volume elements shaped as rectangular
prisms, which can be obtained, e.g., from a CT scan; and (iii) a mixture of
quadrics and voxels. Detailed instructions on how to select one of these
models are given in the file penEasy.in.

To use a voxelized geometry users must provide a valid voxels file. The
format of the voxels file is described in detail in the sample.vox file
provided under ~/voxGeoSample/. Keep this file for future reference.

A tally that reports the absorbed dose in each voxel is also provided. As
with any other tally, a detailed description of its operation is given in
the documentation found in ~/documentation/.

Not all tallies are compatible with voxelized geometries. More precisely,
spatial, spherical and cylindrical dose distributions (obtained with the
tallies bearing those same names) and fluence spectra are not reliable when
reporting data that refers to the voxelized region. See the documentation
of these tallies under ~/documentation/ for more information. To compute
the dose distribution inside the voxelized region, the tally specifically
tailored for voxels should be used.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
6) PARALLEL EXECUTION

Monte Carlo (MC) simulations can be parallelized relatively easily. Various
strategies and tools have been developed to facilitate the implementation
of parallel algorithms, e.g., MPI, openMP, PVM and openMOSIX. All these
have been successfully used with MC, but they suffer from the drawback of
requiring the modification of the sequential code and/or the installation
of additional software, which may be inconvenient for some users.

An alternative solution for penEasy that does not have these limitations is
provided by the package clonEasy, which can be freely downloaded from
http://www.upc.es/inte/downloads/clonEasy.htm . Its principles and usage
are described in [1].

In brief: the same MC job is distributed to several CPUs (the clones) but
different random seeds are provided for each clone. After all the
executions are done, the partial results are collected and averaged
appropriately. To distribute the jobs, the secure shell (ssh) protocol is
used. This protocol is usually embedded in Unix-like systems, including all
GNU/Linux distributions that we are aware of. This means that any
accessible Unix-like computer--e.g. connected to the Internet with a
permanent IP address and with a valid user account--can be a clone for our
parallel computation. This is achieved without the installation of any
additional software.

Random seeds for the PENELOPE pseudo-random number generator (named RANECU)
that initiate disjoint sequences, intended to prevent correlations between
clones, can be easily obtained with the help of the seedsMLCG code. This
code can also be freely downloaded from
http://www.upc.es/inte/downloads/seedsMLCG.htm.

To feed different clones with different seeds it is necessary to exploit
the feature of penEasy that allows the introduction of the initial seed
values through an external file. See the instructions in penEasy.in for
more details.


Reference:

[1] Andreu Badal and Josep Sempau, A package of Linux scripts for the
    parallelization of Monte Carlo simulations, Comput. Phys. Commun. 175
    (2006) pp 440-450.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
7) KNOWN ISSUES AND LIMITATIONS

- All tallies support variance reduction, i.e., they take the statistical
  weight into account when scoring some quantity. However, care must be
  exercised for the tallyEnergyDepositionPulseHeightSpectrum, since some
  variance reduction techniques may bias the calculated distribution. This
  is not a limitation of penEasy, but a consequence of the nature of the
  tally.

- Not all tallies are compatible with voxelized geometries. See the section
  on voxelized geometries for an enumeration of those that should not be
  used.

- The use of PENELOPE with external static electromagnetic fields
  (penfield.f in the PENELOPE package) is not supported. See the PENELOPE
  manual for a description of the changes needed to do so.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
8) WHERE DO I GET THE LATEST UPDATE?

penEasy can be freely downloaded from
http://www.upc.es/inte/downloads/penEasy.htm

It is strongly recommended that you get the latest versions of PENELOPE and
penEasy before attempting a new simulation project.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
9) WHO DO I COMPLAIN TO?

You can send your comments or questions to josep.sempau@upc.es.

Some parts of penEasy (notably penVox, the package that handles voxelized
geometries), have been developed in collaboration with Andreu Badal.

Before sending your query, please make sure that you have the latest
version. Only the newest version is up-to-date; old penEasy releases,
intended for obsolete PENELOPE versions, which may be available for
completeness on the web site, are NOT maintained.


>>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
