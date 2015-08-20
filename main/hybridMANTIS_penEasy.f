!*******************************************************************
!*                         PENEASY                                 *
!*                                                                 *
!* Short description:                                              *
!*   General-purpose PENELOPE main program. Please refer to the    *
!*   README.txt file for detailed instructions.                    *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routines CLEANS,START,SECPAR                               *
!*   from PENVARED:                                                *
!*   -> routine JUMPF                                              *
!*   from other penEasy files:                                     *
!*   -> routines in penaux.f, penpatch.f, penvox.f,                *
!*      sourceXX.f, tallyXX.f and timing.f                         *
!*   In particular, these routines supersede their PENELOPE        *
!*   equivalents:                                                  *
!*   -> KNOCKX (necessary for the fluence tally) replaces KNOCK;   *
!*      it is included in penpatch.f                               *
!*   -> STEPX (from penvox.f) tracks particle trajectories in      *
!*      quadric+voxel geometries, replacing STEP; it is included   *
!*      in penvox.f                                                *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2006                                                          *
!*                                                                 *
!*******************************************************************

!*******************************************************************
!*                                                                 *
!* This version has been derived from the penEasy main program.    *
!* For copyright, refer to penEasy copyright notice in Readme file *
!* under the penEasy_Imaging folder.                               *
!*                                                                 *
!*******************************************************************

!*******************************************************************
!*    Includes                                                     *
!*******************************************************************
      ! PENELOPE routines:     !! PENELOPE files included in Makefile to avoid recompilation
!       include 'penelope.f'
!       include 'pengeom.f'
!       include 'penvared.f'

      ! Auxiliary routines:
      include 'penaux.f'
      include 'penpatch.f'
      include 'penvox.f'
      include 'timing.f'

      ! Source models (see documentation for a detailed description):
      include 'sourceBoxIsotropicGaussSpectrum.f'
      include 'sourcePhaseSpaceFile.f'
      ! <you may add your own here and in the source routines below>

      ! Tallies (see documentation for a detailed description):
      include 'tallyVoxelDoseDistrib.f'
      include 'tallySpatialDoseDistrib.f'
      include 'tallyCylindricalDoseDistrib.f'
      include 'tallySphericalDoseDistrib.f'
      include 'tallyEnergyDepositionPulseHeightSpectrum.f'
      include 'tallyFluenceTrackLength.f'
      include 'tallyPhaseSpaceFile.f'
      include 'tallyParticleCurrentSpectrum.f'
      include 'tallyParticleTrackStructure.f'
      ! <you may add your own here and in the tally routines below>
      

!*******************************************************************
!*    MAIN                                                         *
!*******************************************************************
      program main
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical endsim,absorb,forcing,isforcing
      integer*4 ncross,icol,left
      real*8 n,ds,dsef,de,dsmax

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     real*8 ULDE, ULDV, ULEM, TMAX, emax_track  !!  EM field
              
!      real*8 Zmax, Zmin                     !!  EM field: Z extend of the field, passed from main to GETEMF
!      COMMON/UFIELD_Z/ Zmin, Zmax
      
C  ****  EM field.                          !!  EM field
!      real*8 EFX,EFY,EFZ,BFX,BFY,BFZ
!      COMMON/UFIELD/EFX,EFY,EFZ,BFX,BFY,BFZ

      integer*4 maxvalopt
      parameter (maxvalopt=2304000)
      integer*4 cpu_num_real
      real*8 xbufopt(maxvalopt),ybufopt(maxvalopt)
      real*8 zbufopt(maxvalopt)
      real*8 debufopt(maxvalopt)	!! EDE buffer to be passed to fastDETECT2
      integer*4 nbufopt(maxvalopt)	
      integer*4 myctropt						!! buffer index
      COMMON/optical/xbufopt,ybufopt,zbufopt,debufopt,
     &               nbufopt,myctropt,cpu_num_real


      real*8 glgen, gldetect, glabstop   !! global counters to get optical transport statistics
      real*8 glabsbulk
      real*8 gllost 
      real*8 gloutofcol 
      real*8 gltheta1
      real*8 glgputime
      COMMON/optstats/glgen,gldetect,glabstop,glabsbulk,
     &                gllost,gloutofcol,gltheta1,glgputime

      real*8 detx, dety, detheight
      real*8 detradius, detnC, detnIC
      real*8 dettop, detbulk, detbeta
      real*8 detdmin, detdmax
      real*8 detlboundx, detlboundy 
      real*8 detuboundx, detuboundy
      real*8 detyield
      real*8 detsensorRefl
      integer*4 detpixel
      integer*4 rungpu
      integer*4 machinenum
      integer*4 mynumhist
      integer*4 minphotons,maxphotons
      integer*4 mynumbins
      COMMON/inputargs/detx,dety,detheight,
     &              detradius,detnC,detnIC,
     &              dettop,detbulk,detbeta,detdmin,
     &              detdmax,detlboundx,
     &              detlboundy,detuboundx,
     &              detuboundy,detyield,
     &              detsensorRefl,detpixel,
     &              rungpu,machinenum,
     &              mynumhist,minphotons,
     &              maxphotons,mynumbins

      integer*4 detprimopt		
      parameter (detprimopt=1000)		!! number of BINS
      integer*4 gldetprimary(0:(detprimopt-1))	!! array storing histogram of # detected optical photons
      COMMON/outputdetprim/gldetprimary

      !! output image from fastDETECT2
      integer*8 newimageopt(0:500,0:500) 	
      integer*8 tempimageopt(0:500,0:500)
      COMMON/outputimage/newimageopt,tempimageopt

      integer gpucounter, lbctr
      integer cpuctr
      COMMON/gpuctr/gpucounter,lbctr,cpuctr

      integer*8 gpuimage, gpudetect, hosta, deva
      integer*8 devpitch
      integer*8 depener
      COMMON/gpumemaddr/gpuimage,gpudetect,hosta,deva,
     &                  devpitch,depener

      !! calculate time for fastDETECT2 calls
      real*8 realtime, penopticaltime
      real*8 optstarttime, optendtime
      real*8 gptime, cptime		!! timings for load balancing
      COMMON/tpenoptical/penopticaltime
      COMMON/lbtime/gptime,cptime

      character*10 cdetx
      character*10 cdety, cdetz, cdetrad
      character*10 cnC, cnIC, cTabs, cMUabs
      character*10 cdetrough, cdetdmin, cdetdmax
      character*10 cdetLx,cdetLy,cdetUx,cdetUy
      character*10 cdetyield
      character*10 cdetpp, cdetsensorRefl
      character*10 crungpu,cmachinenum
      character*10 cmynumbins,cmynumhist
      character*10 cminphotons,cmaxphotons
      COMMON/cmdargs/cdetz,cdetrad,cTabs,cMUabs,cdetrough,cdetdmax,
     &               cdetsensorRefl,cdetyield,cdetpp,cmachinenum

      !! load balancing variables
      real*8 factorLB
      integer factorGPU,factorCPU 
      COMMON/lbfactor/factorLB,factorGPU,factorCPU 

      integer*4 newI,newJ
      integer*4 pixelX,pixelY
      character*150 mystring(24)	!! 24 input arguments
      character*150 newstring
      integer myi, myk
  
      !! Read input simulation parameters into strings
      open (11, FILE='hybridMANTIS_input.in',STATUS='OLD')
      myk = 1
      do
        read(11,'(a150)',end=2222) newstring
           newstring = adjustl(newstring)			!! left adjust, removing any leading spaces
        if ((newstring(1:1).ne.'#').and.(newstring(1:1).ne.'')) then  !! Not a comment line
           newstring = newstring(1:scan(newstring,' ')) 	!! Clip at 1st blank, removing any trailing spaces
           mystring(myk) = newstring        
           myk = myk+1
        endif

      enddo

      close(11)
               
      !! convert string to float or integers   
 2222 read(mystring(1),*) mynumhist
      read(mystring(2),*) minphotons  
      read(mystring(3),*) maxphotons
      read(mystring(4),*) mynumbins
      read(mystring(5),*) detx
      read(mystring(6),*) dety
      read(mystring(7),*) detheight
      read(mystring(8),*) detradius
      read(mystring(9),*) detnC
      read(mystring(10),*) detnIC
      read(mystring(11),*) dettop
      read(mystring(12),*) detbulk
      read(mystring(13),*) detbeta
      read(mystring(14),*) detdmin
      read(mystring(15),*) detdmax
      read(mystring(16),*) detlboundx
      read(mystring(17),*) detlboundy
      read(mystring(18),*) detuboundx
      read(mystring(19),*) detuboundy
      read(mystring(20),*) detyield
      read(mystring(21),*) detpixel
      read(mystring(22),*) detsensorRefl
      read(mystring(23),*) rungpu
      read(mystring(24),*) machinenum

     
      !! Write to output file
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'Input arguments:'
      write(*,'(a)')     '' 
      write(*,'(a,a)') 
     & ' NUMBER HIST: ', mystring(1)
      write(*,'(a,a)') 
     & ' MIN OPT PHOT: ', mystring(2)
      write(*,'(a,a)') 
     & ' MAX OPT PHOT: ', mystring(3)
      write(*,'(a,a)') 
     & ' NUMBER BINS: ', mystring(4)
      write(*,'(a,a)') 
     & ' DET X (um): ', mystring(5)
      write(*,'(a,a)') 
     & ' DET Y (um): ', mystring(6)
      write(*,'(a,a)') 
     & ' DET Z (um): ', mystring(7)
      write(*,'(a,a)') 
     & ' RADIUS (um): ', mystring(8)
      write(*,'(a,a)') 
     & ' REF INDEX COL: ', mystring(9)
      write(*,'(a,a)') 
     & ' REF INDEX IC: ', mystring(10)
      write(*,'(a,a)') 
     & ' TOP ABS FRAC: ', mystring(11)
      write(*,'(a,a)') 
     & ' BULK ABS COEFF (um^-1): ',mystring(12)
      write(*,'(a,a)') 
     & ' ROUGHNESS COEFF: ', mystring(13)
      write(*,'(a,a)') 
     & ' MIN DIST TO NEXT COL (um): ',mystring(14)
      write(*,'(a,a)') 
     & ' MAX DIST TO NEXT COL (um): ',mystring(15)
      write(*,'(a,a)') 
     & ' PRF LOWERBOUND X (um): ',mystring(16)
      write(*,'(a,a)') 
     & ' PRF LOWERBOUND Y (um): ', mystring(17)
      write(*,'(a,a)') 
     & ' PRF UPPERBOUND X (um): ',mystring(18)
      write(*,'(a,a)') 
     & ' PRF UPPERBOUND Y (um): ', mystring(19)
      write(*,'(a,a)') 
     & ' YIELD (/eV): ', mystring(20)
      write(*,'(a,a)') 
     & ' PIXEL PITCH (um): ', mystring(21)
      write(*,'(a,a)') 
     & ' NON-IDEAL SENSOR REFL: ', mystring(22)
      write(*,'(a,a)') 
     & ' Run on GPU: ', mystring(23)
      write(*,'(a,a)') 
     & ' Machine num: ', mystring(24)
       write(*,'(a)')     '' 
       write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
     
     !! copy to common block 'cmdargs' to append the arguments in output file names
      cdetz = mystring(7)(1:10)
      cdetrad = mystring(8)(1:10)
      cTabs = mystring(11)(1:10)
      cMUabs = mystring(12)(1:10)
      cdetrough = mystring(13)(1:10)
      cdetdmax = mystring(15)(1:10)
      cdetsensorRefl = mystring(22)(1:10)
      cdetyield = mystring(20)(1:10)
      cdetpp = mystring(21)(1:10)
      cmachinenum = mystring(24)(1:10)

      !! initialize
      myctropt = 1
      cpu_num_real = 0
      glgen = 0.0
      gldetect = 0.0
      glabstop = 0.0
      glabsbulk = 0.0
      gllost = 0.0
      gloutofcol = 0.0
      gltheta1 = 0.0
      glgputime = 0.0     
      gpucounter = 99	!! counter used for gpu initialization(99), run optical transport(100) or run last optical and return final results(101)
      lbctr = 0		!! counter for obtaining the load balancing parameters(lbctr=0), once this is done lbctr=1 and simulation start acc to lb rule
      cpuctr = 0 	!! counter for running factorGPU+1 to factorCPU buffer on CPU
      gpuimage = 0
      gpudetect = 0
      hosta = 0
      deva = 0
      devpitch = 0
      optstarttime = 0.0
      optendtime = 0.0
      penopticaltime = 0.0 
      factorLB = 0.0  
      factorGPU = 0
      factorCPU = 0
     
      pixelX = int((detuboundx-detlboundx)/detpixel)-1
      pixelY = int((detuboundy-detlboundy)/detpixel)-1
      
      DO 12011, newI = 0, pixelX	
       DO 2011, newJ = 0, pixelY
         newimageopt(newI,newJ) = 0
         tempimageopt(newI,newJ) = 0
2011   CONTINUE
12011 CONTINUE

       DO newJ = 0, (mynumbins-1)	
         gldetprimary(newJ) = 0
       ENDDO

       

C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


      write(*,'(a)')' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)')
     & '>>>> This is penEasy v.2008-06-15 >>>>'
      write(*,'(a)')                                  
     & '>>>>   !!NO  EM FIELDS!!  >>>>'       
      write(*,'(a)')                           
     & '>>>> !Voxel geometry disregarded! >>>>' 
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)')' '

      call initime ! Write date on the screen
      
      call init    ! Initialize the penEasy+PENELOPE+PENGEOM system
      call treset  ! Reset simulation timer
      n = 0.0d0    ! Reset history counter

      history: do              ! Each iteration simulates a new history
        n = n+1.0d0            ! Update history counter
        call noverflow(n)      ! Check that N does not overflow
        call cleans            ! Empty the stack, just in case
        call tally(1,n)        ! The simulation of this history begins
        call source(n)         ! Put primary particles (from the same history) in stack

        particle: do                       ! Each iteration simulates a new particle
          call secpar(left)                ! Retrieve a particle from the stack
          if (left.eq.0) exit particle     ! Stack was empty
          call tally(-99,-e)               ! The simulation of this particle begins

!         if(kpar.eq.1 .and. e.gt.emax_track) emax_track=e   !! EM field: tally the maximum electron energy in the track

          if (absorb()) cycle particle     ! Check particle absorption
          call start                       ! Reset transport mechanics
          forcing = isforcing()            ! Set interaction forcing (variance reduct)

          interact: do                     ! Each iteration simulates an interaction
            if (absorb()) exit interact    ! Check particle absorption


!            call TPEMF0(ULDV,ULDE,ULEM,TMAX)     !! EM field. Output: TMAX = maximum allowed step length.
            if (forcing) then
!              call jumpf(min(TMAX,dsmax()),ds)   
                 call jumpf(dsmax(),ds)       
            else
!              call jump(min(TMAX,dsmax()),ds)    
                 call jump(dsmax(),ds)        
            endif
!            call TPEMF1(ds,dsef,ncross)          !! EM field. TPEMF1 will call step; voxels not considered
               call stepx(ds,dsef,ncross)

!            if(kpar.eq.1 .and. e.gt.emax_track) emax_track=e   !! EM field: the electric field may increase the electron's energy!
              
            if (absorb()) exit interact          !! EM field: e- moving against the field will reduce energy, eventually below Eabs.

            if (ncross.eq.0) then
              call tally(3,ds)             ! Moved a distance DS, no interface crossed
            else
              call tally(4,dsef)           ! Moved a distance DSEF, interface found
              if (mat.eq.0) exit interact  ! New material is vacuum => gone
              call start                   ! New material => reset transport mechanics
              forcing = isforcing()        ! Set interaction forcing (variance reduct)
              cycle interact
            endif
            if (forcing) then
              call knockfx(de,icol)        ! Interaction forcing (see PENELOPE manual)
            else
              call knockx(de,icol)         ! Simulate an interaction
            endif
            call tally(-int(icol),de)      ! Tally kinetic energy released
          enddo interact
        enddo particle

        call tally(6,n)                    ! End-of-history bookkeeping


        if (endsim(n)) then			!! if endsim is true, to 0 and call fastDETECT2

          DO newI = myctropt, maxvalopt		!! re-initialize the remaining x,y,z,de arrays
                xbufopt(newI)  = 0.d0
                ybufopt(newI)  = 0.d0
                zbufopt(newI)  = 0.d0
                debufopt(newI) = 0.d0
                nbufopt(newI)  = 0
          ENDDO
       
            if(rungpu.eq.1) then	!! Run on GPU
             gpucounter = 101
             optstarttime = realtime()                   
             call gpuoptical(gpucounter,myctropt)	!! call fastDETECT2 for final remaining buffer (gpucounter = 101)
             optendtime = realtime()

             penopticaltime = penopticaltime +
     &              (optendtime - optstarttime)
           else 			!! Run on CPU
             optstarttime = realtime()                   
             call cpuoptical(0)	
             optendtime = realtime()

             penopticaltime = penopticaltime +
     &              (optendtime - optstarttime)
           endif

          exit history        !! Simulation is finished
        endif

      enddo history

      call report(n)                       !! Write final report

      end


      subroutine init
!*******************************************************************
!*    Initializes                                                  *
!*******************************************************************
      implicit none
      integer*4 nmat
      real*8 emax,realtime,cputime

      call iniconfig              ! Simulation config
      call inisource(emax)        ! Source models
      call inigeo(nmat)           ! Geometry: PENGEOM & penVox


      !!!emax = 100.0d0*emax  !! EM field:  increase emax: e- will get energy from the E field!
      
      call inipen(emax,nmat)      ! PENELOPE
      call initally               ! Tallies
      call iniforce(emax)         ! Interaction forcing
      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'init: INITIALIZATION ENDED'
      write(*,'(a,f9.2,a)') 'Elapsed real time:',realtime(),' s'
      write(*,'(a,f9.2,a)') 'Elapsed CPU time :',cputime(),' s'
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      end


      subroutine report(n)
!*******************************************************************
!*    Reports final results                                        *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*******************************************************************
      implicit none
      integer unc
      real*8 n,cputime,realtime,nowcpu,nowcpu2
      integer*4 seed1,seed2
      common/rseed/seed1,seed2


      real*8 gptime, cptime		
      COMMON/lbtime/gptime,cptime

      nowcpu = cputime()
      call tallyreport(n,nowcpu,unc)

      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'report: SIMULATION ENDED'
      write(*,'(a)')
     & 'Results have been written to the corresponding DAT files.'
      select case(unc)
      case(0)
        write(*,'(a)')
     &   'The requested uncertainty has NOT been reached.'
      case(2)
        write(*,'(a)')
     &   'The requested uncertainty has been reached.'
      end select
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      ! Write generic report to the screen:
      write(*,*) ''
      write(*,'(a)') 'Last random seeds:'
      write(*,'(2(1x,i0))') seed1,seed2

      write(*,'(a)') 'Elapsed real time (s), excluding init:'
      write(*,'(1x,1pe12.5)') realtime()

      nowcpu2 = cputime()
      write(*,'(a)') 'Elapsed CPU time (s), excluding init:'
      write(*,'(1x,1pe12.5)') nowcpu2

      write(*,'(a)') 'Each report update took (in CPU s):'
      write(*,'(1x,1pe12.5)') nowcpu2-nowcpu

      write(*,'(a)') 'No. of histories simulated:'
      write(*,'(1x,f18.0)') n

      if (nowcpu.gt.0.0d0) then
        write(*,'(a)') 'CPU Speed (histories/s):'
        write(*,'(1x,1pe12.5)') n/nowcpu

        write(*,'(a)') 
     &  'CPU Speed without load balancing time (histories/s):'
        write(*,'(1x,1pe12.5)') n/(nowcpu-cptime)
      endif

      call endtime  ! Report timings
      end


!*******************************************************************
!*******************************************************************
!*    Source routines start here.                                  *
!*    Add your own models or delete the unwated ones.              *
!*    Source models usually require:                               *
!*     i) an initialization routine that must be called by         *
!*        INISOURCE                                                *
!*    ii) a particle generation routine that must be called        *
!*        by SOURCE                                                *
!*******************************************************************
!*******************************************************************

      subroutine inisource(emax)
!*******************************************************************
!*    Init routines for source models                              *
!*                                                                 *
!*    Output:                                                      *
!*      emax -> max source energy (eV)                             *
!*******************************************************************
      implicit none
      real*8 emax
      logical active
      integer nsrc

      nsrc = 0
      call BIGSinisrc(active,emax)
      if (active) nsrc = nsrc+1
      call PSFinisrc(active,emax)
      if (active) nsrc = nsrc+1

      if (nsrc.ne.1) then
        write(*,*) ''
        write(*,'(a)')
     &    'inisource:ERROR: There must be ONE active source'
        stop
      endif
      end


      subroutine source(n)
!*******************************************************************
!*    Source models                                                *
!*                                                                 *
!*    Input:                                                       *
!*      n -> top history counter                                   *
!*******************************************************************
      implicit none
      real*8 n

      call BIGSsource(n)
      call PSFsource(n)
      end


!*******************************************************************
!*******************************************************************
!*    Tally routines start here.                                   *
!*    Add your own tallies or delete the unwated ones              *
!*    Tallies usually require:                                     *
!*    i) an initialization routine that must be called by INITALLY *
!*    ii) a tally routine that must be called by TALLY             *
!*    iii) a reporting routine that must be called by TALLYREPORT  *
!*                                                                 *
!*    Notice that the ordering of the tally initialization routines*
!*    must coincide with the ordering of the corresponding sections*
!*    in the input file.                                           *
!*******************************************************************
!*******************************************************************

      subroutine initally
!*******************************************************************
!*    Init tallying routines.                                      *
!*                                                                 *
!*    Comments:                                                    *
!*      -> VDDinitally sets variables that are needed for          *
!*         proper particle transport in voxelized geometries.      *
!*         Therefore, it should not be deactivated even if the     *
!*         tally is not used.                                      *
!*******************************************************************
      implicit none
      call VDDinitally
      call SDDinitally
      call CDDinitally
      call SPDinitally
      call EPSinitally
      call FTLinitally
      call PSFinitally
      call PCSinitally
      call PTSinitally

      call EDEinitally  !! Calling the new EDE tally!

      end


      subroutine tally(mode,arg)
!*******************************************************************
!*    Tallying routines.                                           *
!*                                                                 *
!*    Comments:                                                    *
!*      -> VDDtally sets state variables that are needed for       *
!*         proper particle transport in voxelized geometries.      *
!*         Therefore, it should not be deactivated even if the     *
!*         tally is not used.                                      *
!*      -> Furthermore, these variables could be used by other     *
!*         tallies and, in consequence, VDDtally should be the     *
!*         first tally to be called.                               *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      call VDDtally(mode,arg)  ! Must be 1st tally to be called
      call SDDtally(mode,arg)
      call CDDtally(mode,arg)
      call SPDtally(mode,arg)
      call EPStally(mode,arg)
      call FTLtally(mode,arg)
      call PSFtally(mode,arg)
      call PCStally(mode)
      call PTStally(mode,arg)

      call EDEtally(mode,arg)  !! Calling the new EDE tally!
      
      end


      subroutine tallyreport(n,cputim,unc)
!*******************************************************************
!*    Calls report routines for all tallies                        *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      simtim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      unc -> larger than 1 if requested uncert has been reached  *
!*******************************************************************
      integer unc,uncdone
      real*8 n,cputim

      ! Write partial reports to corresponding data files:
      unc = 1
      call VDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call SDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call CDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call SPDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call EPSreport(n,cputim,uncdone)
      unc = unc*uncdone
      call FTLreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PSFreport(n,cputim)           ! No uncertainty for this tally
      call PCSreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PTSreport                     ! No arguments for this tally

      call EDEreport(n)  !! Calling the new EDE tally! 
      
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
