!*******************************************************************
!*                          TALLY                                  *
!*         ENERGY DEPOSITION AND PULSE HEIGHT SPECTRUM             *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Calculates the spectrum of the energy deposited in a given    *
!*   material and the total energy deposition. The so-called       *
!*   collision estimator is employed.                              *
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


      subroutine EPStally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      eloss -> energy deposition                                 *
!*******************************************************************
      implicit none
      integer mode
      real*8 eloss

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 edptmp,count,emin,ebin,iebin,edep,edep2,unclimit
      common /scoeps/ count(maxbin),edptmp,edep,edep2,emin,ebin,
     &                iebin,unclimit,nbin,detmat,active
      integer bin

      if (.not.active) return

      if (mode.le.0) then
        if (mat.ne.detmat) return
        edptmp = edptmp+eloss*wght

      else if (mode.eq.6) then
        edep   = edep +edptmp
        edep2  = edep2+edptmp**2
        bin = (edptmp-emin)*iebin+1.0d0
        edptmp = 0.0d0
        if (bin.lt.1.or.bin.gt.nbin) return
        count(bin) = count(bin)+1.0d0

      endif
      end


      subroutine EPSreport(n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      uncdone -> 2 if uncert reached, 1 if not defined, 0 else   *
!*******************************************************************
      implicit none
      integer uncdone
      real*8 n,cputim

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 edptmp,count,emin,ebin,iebin,edep,edep2,unclimit
      common /scoeps/ count(maxbin),edptmp,edep,edep2,emin,ebin,
     &                iebin,unclimit,nbin,detmat,active
      integer i,nchan,out,finduf,error
      real*8 q,q2,sigma,eff,avesig,maxq,emiddle,elow,uncert,invn

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      !*************************************
      !*** Total energy deposited report ***
      !*************************************

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyEnergyDeposition.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'EPSreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY DEPOSITION]'
      write(out,'(a)') '# Total Energy Deposition(eV) : +-2sigma'
      q = edep*invn
      q2 = edep2*invn
      sigma = (q2-q**2)*invn
      sigma = sqrt(max(sigma,0.0d0))
      write(out,'(1pe12.5,1x,1pe8.1)') q,2.0d0*sigma
      ! Evaluate rel. uncertainty:
      uncert = 200.0d0
      if (q.gt.0.0d0) uncert = 200.0d0*sigma/q

      ! Generic report:
      write(out,'(a)') ' '
      write(out,'(a)') '# Performance report'
      write(out,'(a)') '#   Random seeds:'
      write(out,'(a,i10)') '#   ',seed1
      write(out,'(a,i10)') '#   ',seed2
      write(out,'(a)') '#   No. of histories simulated [N]:'
      write(out,'(a,f18.0)') '#   ',n
      write(out,'(a)') '#   CPU time [t] (s):'
      write(out,'(a,1pe12.5)') '#   ',cputim
      if (cputim.gt.0.0d0) then
        write(out,'(a)') '#   Speed (histories/s):'
        write(out,'(a,1pe12.5)') '#   ',n/cputim
      endif
      write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
      write(out,'(a,1pe12.5)') '#   ',uncert
      eff = n*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      eff = cputim*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      close(out)

      !***********************
      !*** Spectrum report ***
      !***********************

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyPulseHeightSpectrum.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'EPSreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Evaluate 1/2 of the max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do i=1,nbin
        maxq = max(maxq,count(i))
      enddo
      maxq = 0.5d0*maxq

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PULSE HEIGHT SPECTRUM]'
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# For plotting purposes, two energies per bin are given,'
      write(out,'(a)')
     &  '#   namely, the low end and the middle point of each bin'
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# Elow(eV) : Emiddle(eV) : counts(1/eV) : +-2sigma'
      do i=1,nbin
        q = count(i)*invn
        sigma = q*(1.0d0-q)*invn
        sigma = sqrt(max(sigma,0.0d0))
        q = q*iebin
        sigma = sigma*iebin
        elow = emin+ebin*(i-1)
        emiddle = elow+ebin*0.5d0
        write(out,'(3(1x,1pe12.5),1x,1pe8.1)')
     &    elow,emiddle,q,2.0d0*sigma
        ! Evaluate average uncertainty for scores above 1/2 max score:
        if (count(i).gt.maxq) then
          avesig = avesig+(sigma/q)**2
          nchan = nchan+1
        endif
      enddo
      uncdone = 0
      if (nchan.gt.0) then
        uncert = 200.0d0*sqrt(avesig/nchan)
        if (uncert.lt.unclimit) uncdone = 2  ! Uncertainty reached
      else
        uncert = 0.0d0  ! Uncertainty assumed not reached when score is nil
      endif

      ! Generic report:
      write(out,'(a)') ' '
      write(out,'(a)') '# Performance report'
      write(out,'(a)') '#   Random seeds:'
      write(out,'(a,i10)') '#   ',seed1
      write(out,'(a,i10)') '#   ',seed2
      write(out,'(a)') '#   No. of histories simulated [N]:'
      write(out,'(a,f18.0)') '#   ',n
      write(out,'(a)') '#   CPU time [t] (s):'
      write(out,'(a,1pe12.5)') '#   ',cputim
      if (cputim.gt.0.0d0) then
        write(out,'(a)') '#   Speed (histories/s):'
        write(out,'(a,1pe12.5)') '#   ',n/cputim
      endif
      write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
      write(out,'(a,1pe12.5)') '#   ',uncert
      eff = n*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      eff = cputim*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      close(out)
      end


      subroutine EPSinitally
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*******************************************************************
      implicit none
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 edptmp,count,emin,ebin,iebin,edep,edep2,unclimit
      common /scoeps/ count(maxbin),edptmp,edep,edep2,emin,ebin,
     &                iebin,unclimit,nbin,detmat,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY ENERGY DEPOSITION PULSE SPECTRUM v.2006-08-01]')
      parameter (eos='[END OF EPS SECTION]')
      character*80 buffer
      integer i,error
      real*8 emax

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'EPSinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Energy Deposition Pulse Height Spectrum is OFF'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'EPSinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'EPSinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) detmat
      write(*,'(i3)') detmat
      write(*,'(a)')
     &  'Spectrum limits [Emin,Emax] (eV) and no. bins:'
      read(*,*) emin,emax,nbin
      write(*,'(2(1x,1pe12.5),1x,i5)') emin,emax,nbin
      ebin = (emax-emin)/nbin
      iebin = 1.0d0/ebin
      ! Add one bin to store, separately, counts with E=Emax:
      nbin = nbin+1
      if (nbin.gt.maxbin) then
        write(*,*)
     &    'EPSinitally:ERROR: Too many bins; increase MAXBIN'
        stop
      endif
      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,1pe12.5)') unclimit

      ! Clear counters:
      edptmp = 0.0d0
      edep   = 0.0d0
      edep2  = 0.0d0
      do i=1,nbin
        count(i) = 0.0d0
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'EPSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> EPS tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
