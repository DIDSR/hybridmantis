!*******************************************************************
!*                          TALLY                                  *
!*                 ENERGY DEPOSITION EVENTS                        *
!*                                                                 *
!* Short description:                                              *
!*                                                                 *
!*                                                                 *
!*                                                                 *
!*   Writes to a file the position and value of energy loss events *
!*   so that particle tracks can be visualized.                    *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   from other penEasy files:                                     *
!*   -> routine GETLINE,FINDUF                                     *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*******************************************************************

      subroutine EDEtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      implicit none
      integer mode, i
      real*8 arg, de
      real*8 n

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)

C  ****  Secondary stack (PENELOPE).
      integer*4 KS,IBODYS,MS,ILBS,NSEC,NMS
      real*8 ES,XS,YS,ZS,US,VS,WS,WGHTS
      PARAMETER (NMS=2000)
      COMMON/SECST/ES(NMS),XS(NMS),YS(NMS),ZS(NMS),US(NMS),
     1   VS(NMS),WS(NMS),WGHTS(NMS),KS(NMS),IBODYS(NMS),MS(NMS),
     2   ILBS(5,NMS),NSEC


      logical active
      integer edeunit
      integer*4 num_ede, num_real, mat_detector, secondaries_bef_knock
      real*8 energy_before_knock, total_edep_history, distance_step
      common /comEDE/ energy_before_knock, total_edep_history,
     &                distance_step, secondaries_bef_knock,
     &                num_ede, num_real, mat_detector, edeunit, active

      integer*4 maxvalopt
      parameter (maxvalopt=2304000)
      integer*4 lbbuf
      parameter (lbbuf=100000)
      integer*4 cpubufsize
      parameter (cpubufsize=100000)	!! lbbuf should be >= cpubufsize
      integer*4 cpu_num_real
      real*8 xbufopt(maxvalopt),ybufopt(maxvalopt)
      real*8 zbufopt(maxvalopt)
      real*8 debufopt(maxvalopt)	!! EDE buffer to be passed to fastDETECT2 routines
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
      integer*4 gldetprimary(0:(detprimopt-1))	!! array for storing histogram of # detected optical photons
      COMMON/outputdetprim/gldetprimary

      !! output PRF image
      integer*8 newimageopt(0:500,0:500) 	
      integer*8 tempimageopt(0:500,0:500)
      COMMON/outputimage/newimageopt,tempimageopt

      integer gpucounter, lbctr
      integer cpuctr
      COMMON/gpuctr/gpucounter,lbctr,cpuctr

      integer*8 gpuimage, gpudetect, hosta, deva
      integer*8 devpitch
      COMMON/gpumemaddr/gpuimage,gpudetect,hosta,deva,
     &                  devpitch

      !! calculate time for fastDETECT2 calls
      real*8 realtime, penopticaltime
      real*8 optstarttime, optendtime
      real*8 gptime, cptime		!! timings for load balancing
      real*8 bufspeed,cpuspeed,gpuspeed	
      COMMON/tpenoptical/penopticaltime
      COMMON/lbtime/gptime,cptime

      !! load balancing variables
      real*8 factorLB
      real*8 factor11
      integer*4 factorGPU,factorCPU 
      COMMON/lbfactor/factorLB,factorGPU,factorCPU     
 

       if (.not.active) return

      if (mode.eq.-99) then
        ! ** Particle track begins (primary or secondary):        
        if (ilb(3).eq.5) then
          ! -- This particle was created in an inner shell ionization event: do not simulate.   
          e = 0.0
        else          
!          write(edeunit,'(a)')" "     ! Separate particles
        endif  

      else if (mode.eq.3) then
        ! ** An interaction will take place:
        secondaries_bef_knock = NSEC   ! Store the number of particles in secondary stack
        energy_before_knock = e
        distance_step = arg
        
      else if (mode.lt.0     .and.
     &          arg.gt.0.0d0 .and.           ! (Neglect elastic and inner shell events)
     &          mat.eq.mat_detector) then
        ! ** An interaction happened with energy deposition inside the detector:
        num_ede = num_ede + 1
        
        ! Subtracting the energy of the new secondaries from ARG:
        de = arg
        do i = (secondaries_bef_knock+1), NSEC
          de = de - ES(i)   ! Subtract the energy of the new secondaries            
        enddo

        total_edep_history = total_edep_history + de
        

		!- if the interaction created secondary particles but did not deposit energy 
		if (de.gt.0.0d0) then 

!	        write(edeunit,'(3(1x,1pe16.9),1x,1pe12.4,i8)')
!     &    x,y,z, de, num_real+1
  
                  xbufopt(myctropt)  = x		!! copying to EDE buffer to be passed to fastDETECT2 routines
                  ybufopt(myctropt)  = y
                  zbufopt(myctropt)  = z
                  debufopt(myctropt) = de
                  nbufopt(myctropt)  = num_real+1

                   if(myctropt.eq.cpubufsize) then
                      cpu_num_real = num_real + 1
                   endif

                   if((rungpu.eq.1).and.(lbbuf.eq.0)
     &                 .and.(cpubufsize.eq.0)) then	!! NO LOAD BALANCING - everything goes to GPU 
			 factorGPU = maxvalopt
			 factorCPU = 0
			 lbctr = 1		!! Run on the GPU
			 cpuctr = 0
		   endif


                   if((rungpu.eq.1).and.(lbctr.eq.0)
     &                 .and.(cpuctr.eq.0)) then	!! Load Balancing - only when rungpu is ON         

                       if (myctropt.eq.lbbuf) then	!! buffer equals gpu buffer size, then call gpuoptical !!
                            optstarttime = realtime()

                            call gpuoptlb(gptime,lbbuf)			!! run optical transport(lbbuf bufsize) on GPU
                            call cpuoptlb(cptime,cpubufsize)		!! run optical transport(with cpubufsize) on CPU (LBBUF >= CPUBUFSIZE)

                            optendtime = realtime()
                            penopticaltime = penopticaltime + 
     &                            (optendtime - optstarttime)

                            bufspeed = 794		!! Hard coded - running only PENELOPE on 1 cpu core
                            gpuspeed = (num_real+1)/gptime
			    cpuspeed = (cpu_num_real)/cptime

                           !! LOAD BALANCING - USING 1 CPU CORE + 1 GPU
                           factorLB = ((cpuspeed/bufspeed) + 1)/
     &                            ((cpuspeed/gpuspeed) + 1)

			   if(factorLB.gt.1.0) factorLB = 1.0 	!! factorLB cannot be > 1.0, because gpu cannot be greater than buffer size
                        		                    
                        		                    
                           factor11 = factorLB * maxvalopt
                           factorGPU = ceiling(factor11)
                           factorCPU = maxvalopt-factorGPU

                           write(*,'(a)')' '
                           write(*,'(a)')
     & 			   '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                           write(*,'(a)')
     & 			   '    LOAD BALANCING STATISTICS'
                           write(*,'(a)')' '
		           write(*,'(a)') 'GPU time (s)'
    	                   write(*,'(1x,F12.2)') gptime
		           write(*,'(a)') 'CPU time (s)'
    	                   write(*,'(1x,F12.2)') cptime
		           write(*,'(a)') 'PENELOPE speed (hist/sec)'
    	                   write(*,'(1x,F12.2)') bufspeed
		           write(*,'(a)') 'fastDETECT2 GPU speed (hist/sec)'
    	                   write(*,'(1x,F12.2)') gpuspeed
		           write(*,'(a)') 'fastDETECT2 CPU speed (hist/sec)'
    	                   write(*,'(1x,F12.2)') cpuspeed
		           write(*,'(a)') 'Load Balancing factor'
    	                   write(*,'(1x,F12.2)') factorLB
		           write(*,'(a)') 'Buffer sent to GPU (%)'
    	                   write(*,'(1x,I8)') (factorGPU/maxvalopt)*100
		           write(*,'(a)') 'Buffer sent to CPU (%)'
    	                   write(*,'(1x,I8)') (factorCPU/maxvalopt)*100

                            myctropt = 1		!! reset !!
                            lbctr = 1			!! Run rest of the simulation with LB rule
                       else
                            myctropt = myctropt + 1
                       endif  

                   else if((rungpu.eq.1).and.(lbctr.eq.1)
     &                  .and.(cpuctr.eq.0)) then	!! Run on GPU according to LB rule

                       if (gpucounter.eq.99) then	!! calling optical first time - for initializing gpu !!
                            optstarttime = realtime()
                            call gpuoptical(gpucounter,factorGPU)	!! initialize gpu
                            optendtime = realtime()
                            penopticaltime = penopticaltime + 
     &                            (optendtime - optstarttime)

                            myctropt = myctropt + 1
                            gpucounter = 100	!! now that the gpu is initialized, next time it should run optical transport
                       else if (myctropt.eq.maxvalopt) then	!! buffer full, call gpuoptical !!
                            optstarttime = realtime()		
                            call gpuoptical(gpucounter,factorGPU)	!! run optical transport [1,factorGPU]
                            optendtime = realtime()
                            penopticaltime = penopticaltime + 
     &                            (optendtime - optstarttime) 

			    !! myctropt will not be reset because the rest of buffer will be run on CPU
                            cpuctr = 1
                       else
                            myctropt = myctropt + 1
                       endif

		   endif !! if rungpu=1,lbctr=0/1,cpuctr=0 loop ends

                   if((rungpu.eq.1).and.(lbctr.eq.1)
     &                   .and.(cpuctr.eq.1)) then			!! Run on CPU according to LB rule

                       if ((myctropt.eq.maxvalopt).and.
     &                    (factorCPU.gt.0)) then	!! calling C optical transport function !!
                            optstarttime = realtime()
                            call cpuoptical(factorGPU) !! C has to run hist from [factorGPU+1,maxvalopt]
                            optendtime = realtime()

                            penopticaltime = penopticaltime + 
     &                            (optendtime - optstarttime)

                            myctropt = 1	!! reset !!
                            cpuctr = 0		!! reset this so for next buffer it goes first to GPU anf then come here
                       else			!! enter here when factorCPU=0, the entire buffer was processed at gpuoptical
                            myctropt = 1	
			    cpuctr = 0
                       endif

		   endif !! if rungpu=1,lbctr=1,cpuctr=1 loop ends			

                   if (rungpu.eq.0) then			!! Run optical on the CPU - no GPUs involved

                       if (myctropt.eq.maxvalopt) then	!! calling C optical transport function !!
                            optstarttime = realtime()
                            call cpuoptical(0)	!! C has to run hist from [1,maxvalopt]
                            optendtime = realtime()

                            penopticaltime = penopticaltime + 
     &                            (optendtime - optstarttime)

                            myctropt = 1	!! reset !!
                       else
                            myctropt = myctropt + 1
                       endif

                    endif !! if rungpu=0 loop ends

                endif !! (de.gt.0.0d0) ends
	
	     
      else if (mode.eq.6) then
   		num_real = num_real + 1
     
        write(edeunit,'(a,i8,a,1pe16.4)')"#  Past hist = ",
     &        num_real, " ; Total edep = ", total_edep_history  
        write(edeunit,'(a)')" "
        write(edeunit,'(a)')" "
        
2010    total_edep_history = 0.0d0    ! Reset the history edep counter
     

      endif
      
      end


      subroutine EDEreport(n)
!*******************************************************************
!*    Report number of deposition events per primary               *
!*******************************************************************
      implicit none

      real*8 n, inv_n

      logical active
      integer edeunit
      integer*4 num_ede, num_real, mat_detector, secondaries_bef_knock
      real*8 energy_before_knock, total_edep_history, distance_step
      common /comEDE/ energy_before_knock, total_edep_history,
     &                distance_step, secondaries_bef_knock,
     &                num_ede, num_real, mat_detector, edeunit, active

      real*8 glgen, gldetect, glabstop   !! global counters to get optical transport statistics !!
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
      parameter (detprimopt=1000)		!! number of BINS MAX
      integer*4 gldetprimary(0:(detprimopt-1))	!! storing histogram of # detected optical photons
      COMMON/outputdetprim/gldetprimary

      !! PRF image
      integer*8 newimageopt(0:500,0:500) 	
      integer*8 tempimageopt(0:500,0:500)
      COMMON/outputimage/newimageopt,tempimageopt

      !! timings for fastDETECT2 calls
      real*8 penopticaltime
      COMMON/tpenoptical/penopticaltime

      character*10 cdetz, cdetrad
      character*10 cTabs, cMUabs
      character*10 cdetrough, cdetdmax
      character*10 cdetsensorRefl
      character*10 cdetyield, cdetpp
      character*10 cmachinenum
      COMMON/cmdargs/cdetz,cdetrad,cTabs,cMUabs,cdetrough,cdetdmax,
     &               cdetsensorRefl,cdetyield,cdetpp,cmachinenum

      !! load balancing variables
      real*8 factorLB
      integer*4 factorGPU,factorCPU 
      COMMON/lbfactor/factorLB,factorGPU,factorCPU 

      integer*4 newI,newJ
      integer*4 pixelX,pixelY
      character*100 prfname, detname     
    
      pixelX = int((detuboundx-detlboundx)/detpixel)-1
      pixelY = int((detuboundy-detlboundy)/detpixel)-1

      prfname = 'myimage_'//trim(cdetz)//'_'//trim(cdetrad)//'_'//
     &          trim(cTabs)//'_'//trim(cMUabs)//'_'//trim(cdetrough)//
     &          '_'//trim(cdetdmax)//'_'//trim(cdetsensorRefl)//
     &          '_'//trim(cdetyield)//'_'//trim(cdetpp)//
     &          '_'//trim(cmachinenum)//'.dat'

      detname = 'detected_'//trim(cdetz)//'_'//trim(cdetrad)//'_'//
     &          trim(cTabs)//'_'//trim(cMUabs)//'_'//trim(cdetrough)//
     &          '_'//trim(cdetdmax)//'_'//trim(cdetsensorRefl)//
     &          '_'//trim(cdetyield)//'_'//trim(cdetpp)//
     &          '_'//trim(cmachinenum)//'.dat'


      open (2, FILE=prfname,STATUS='NEW')	      !! open file for writing output image
      open (3, FILE=detname,STATUS='NEW')	      
    

      write(*,'(a)') " "
      write(*,'(a)')   ">>>>>>>>> EDEreport >>>>>>>>>"
      write(*,'(a,i8)')"    Number of energy deposition"//
     &                 " events written to the output file: ", num_ede
      write(*,'(a,1pe12.5)')"    Depositions per primary = ",
     &                           dble(num_ede)/n

      inv_n = 1.0d0/n

!      write out the final image (prf) to the file
        do 1010 newJ = 0, pixelX		!! write output image to file
          do 2020 newI = 0, pixelY
            write(2,'(1x,F12.4)') newimageopt(newJ,newI) * inv_n    !! not normalizing it here, because it does not get properly normalized (max~1.08)
2020       continue
            write(2,'(a/)')
1010    continue

!	write the # detected per primary
          do newI = 0, (mynumbins-1)
            write(3,'(1x,I8)') gldetprimary(newI)   
          enddo

!      report the Optical Transport statistics
        write(*,'(a)')' '
        write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(*,'(a)')
     &	'    OPTICAL TRANSPORT STATISTICS'
        write(*,'(a)')' '

      if(rungpu.eq.1) then
        write(*,'(a)')
     &  '   !! Simulation in the GPU !!'
      else
        write(*,'(a)')
     &  '   !! Simulation in the CPU !!'
      endif

        write(*,'(a)')' '
        write(*,'(a)') 'Total # optical photons:'
        write(*,'(a)') 'Generated'
	write(*,'(1x,F12.2)') glgen 
	write(*,'(a)') 'Detected'
	write(*,'(1x,F12.2)') gldetect
        write(*,'(a)') 'Absorbed at top'
        write(*,'(1x,F12.2)') glabstop
        write(*,'(a)') 'Absorbed in bulk'
        write(*,'(1x,F12.2)') glabsbulk
        write(*,'(a)') 'Lost'
        write(*,'(1x,F12.2)') gllost
        write(*,'(a)') 'Out of Column'
        write(*,'(1x,F12.2)') gloutofcol
	write(*,'(a)') 'Total GPU or CPU time (s)'
	write(*,'(1x,F12.2)') glgputime*0.001
	write(*,'(a)') 'Total optical real time from PENELOPE (s)'
	write(*,'(1x,F12.2)') penopticaltime
        write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(*,'(a)')' '

        close (2)	!! close files
        close (3)
        close (4)
     
      end


      subroutine EDEinitally
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      implicit none
      
      logical active
      integer edeunit
      integer*4 num_ede, num_real, mat_detector, secondaries_bef_knock
      real*8 energy_before_knock, total_edep_history, distance_step
      common /comEDE/ energy_before_knock, total_edep_history,
     &                distance_step, secondaries_bef_knock,
     &                num_ede, num_real, mat_detector, edeunit, active

      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY ENERGY DEPOSITION EVENTS v.2010-06-03]')
      parameter (eos='[END OF EDE SECTION]')
      character*80 buffer
      integer finduf,error

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'EDEinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(1x,a3)') buffer
      if (adjustl(buffer(1:3)).eq.'ON') then
        active = .true.
      else if (buffer(1:3).eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Particle Track Structure is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'EDEinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'EDEinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Init counters:
      total_edep_history  = 0.0d0
      energy_before_knock = 0.0d0
      secondaries_bef_knock = 0
      num_ede = 0
      num_real = 0

      write(*,'(a)') ' Detector sensitive material number:'
      read(*,*) mat_detector
      write(*,'(2x,i3)') mat_detector

      write(*,'(a)') ' EDE output file name:'
      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      write(*,'(1x,a)') buffer

      edeunit = finduf()
      open(edeunit,file=buffer,iostat=error)
      if (error.ne.0) then
        write(*,'(a)')
     &    'EDEinitally:ERROR: cannot open output data file:'
        write(*,'(a80)') buffer
        stop
      endif
      
      write(edeunit,'(a)') "# Output penEasy:"
!     &     " [SECTION TALLY ENERGY DEPOSITION EVENTS v.2010-06-03]"
!      write(edeunit,'(a)') "# "
!      write(edeunit,'(a)') 
!     &   "#  x : y : z : ds : e0_primary : e_second "

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'EDEinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> EDE tally initialization finished >>>>'
      
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

