!*******************************************************************
!*                          TALLY                                  *
!*                    PARTICLE CURRENT SPECTRUM                    *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Determines the spectrum and total number of particles (grouped*
!*   according to their electrical charge) that hit a given        *
!*   material.                                                     *
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


      subroutine PCStally(mode)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*******************************************************************
      implicit none
      integer mode

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 countmp,count,count2,tcount2
      real*8 edeptmp,edep,edep2,tedep2,emin,ebin,iebin,unclimit
      common /scopcs/
     & countmp(3,maxbin),count(3,maxbin),count2(3,maxbin),tcount2(3),
     & edeptmp(3),edep(3),edep2(3),tedep2,
     & emin,ebin,iebin,unclimit,nbin,detmat,active
      integer bin,k,i
      real*8 totedep,totcount

      if (.not.active) return

      if (mode.eq.4.or.mode.eq.-99) then
        if (mat.ne.detmat) return
        edeptmp(kpar) = edeptmp(kpar)+e*wght
        bin = (e-emin)*iebin+1.0d0
        if (bin.lt.1.or.bin.gt.nbin) return
        countmp(kpar,bin) = countmp(kpar,bin)+wght

      else if (mode.eq.6) then
        totedep = 0.0d0
        do k=1,3
          totcount = 0.0d0
          do i=1,nbin
            totcount    = totcount   +countmp(k,i)
            count(k,i)  = count(k,i) +countmp(k,i)
            count2(k,i) = count2(k,i)+countmp(k,i)**2
            countmp(k,i) = 0.0d0
          enddo
          tcount2(k) = tcount2(k)+totcount**2
          totedep  = totedep +edeptmp(k)
          edep(k)  = edep(k) +edeptmp(k)
          edep2(k) = edep2(k)+edeptmp(k)**2
          edeptmp(k) = 0.0d0
        enddo
        tedep2 = tedep2+totedep**2

      endif
      end


      subroutine PCSreport(n,cputim,uncdone)
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
      real*8 countmp,count,count2,tcount2
      real*8 edeptmp,edep,edep2,tedep2,emin,ebin,iebin,unclimit
      common /scopcs/
     & countmp(3,maxbin),count(3,maxbin),count2(3,maxbin),tcount2(3),
     & edeptmp(3),edep(3),edep2(3),tedep2,
     & emin,ebin,iebin,unclimit,nbin,detmat,active
      character*80 buffer
      integer kpar,nchan,out,finduf,error,i
      real*8 q,q2,sigma,eff,avesig,maxq,emiddle,elow,uncert,qtot,invn

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      !*********************************
      !*** Particle Current Spectrum ***
      !*********************************

      do kpar=3,1,-1    ! One report for each particle type
        ! Prepare output files:
        out = finduf()
        if (kpar.eq.3)
     &  buffer = 'tallyParticleCurrentSpectrum-positron.dat'
        if (kpar.eq.2)
     &  buffer = 'tallyParticleCurrentSpectrum-photon.dat'
        if (kpar.eq.1)
     &  buffer = 'tallyParticleCurrentSpectrum-electron.dat'
        open(out,file=buffer,iostat=error)
        if (error.ne.0) then
          write(*,*)
          write(*,'(a)')
     &      '*********************************************'
          write(*,'(a)')
     &      'PCSreport:ERROR: cannot open output data file'
          write(*,'(a)')
     &      '*********************************************'
          close(out)  ! Just in case
          return
        endif

        ! Evaluate 1/2 of the max fluence score:
        avesig = 0.0d0
        nchan = 0
        maxq = 0.0d0
        do i=1,nbin
          maxq = max(maxq,count(kpar,i))
        enddo
        maxq = 0.5d0*maxq

        write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)')
     &    '# [SECTION REPORT PARTICLE CURRENT SPECTRUM]'
        if (kpar.eq.3) buffer = '# POSITRONS'
        if (kpar.eq.2) buffer = '# PHOTONS'
        if (kpar.eq.1) buffer = '# ELECTRONS'
        write(out,'(a)') buffer
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# Energy spectrum, per primary history, of particles'
        write(out,'(a)')
     &    '#   entering the detector'
        write(out,'(a,i3)') '# Detection material: ',detmat
        write(out,'(a)') '#'
        write(out,'(a)')
     &  '# For plotting purposes, two energies per bin are given,'
        write(out,'(a)')
     &  '#   namely, the low end and the middle point of each bin'
        write(out,'(a)') '#'
        write(out,'(a)')
     &  '# Elow(eV) : Emiddle(eV) : counts(1/eV) : +-2sigma'
        qtot = 0.0d0
        do i=1,nbin
          q = count(kpar,i)*invn
          qtot = qtot+q
          q2 = count2(kpar,i)*invn
          sigma = (q2-q**2)*invn
          sigma = sqrt(max(sigma,0.0d0))
          q = q*iebin
          sigma = sigma*iebin
          elow = emin+ebin*(i-1)
          emiddle = elow+ebin*0.5d0
          write(out,'(3(1x,1pe12.5),1x,1pe8.1)')
     &      elow,emiddle,q,2.0d0*sigma
          ! Evaluate average uncertainty for scores above 1/2 max score:
          if (count(kpar,i).gt.maxq) then
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

        ! Total (integrated) counter:
        q2 = tcount2(kpar)*invn
        sigma = (q2-qtot**2)*invn
        sigma = sqrt(max(sigma,0.0d0))
        write(out,'(a)') '#'
        write(out,'(a)')
     &  '# Total particle current (per history) : +-2sigma:'
        write(out,'(a,1pe12.5,1x,1pe8.1)') '# ',qtot,2.0d0*sigma

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
      enddo

      !**********************
      !*** Energy current ***
      !**********************

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyParticleCurrentSpectrum-totalE.dat',
     &     iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'PCSreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY CURRENT]'
      write(out,'(a)') '# Description:'
      write(out,'(a)')
     &  '#   Sum of energies, per primary history, of particles'
      write(out,'(a)')
     &  '#   entering the detector'
      write(out,'(a,i3)') '# Detection material: ',detmat
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# Particle_charge=-1,0,+1 and 3 indicate electrons, photons,'
      write(out,'(a)')
     &  '#   positrons and the sum, respectively.'
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# Particle_charge : Total_energy(eV) : +-2sigma'
      qtot  = 0.0d0
      do kpar=1,3
        q = edep(kpar)*invn
        qtot = qtot+q
        q2 = edep2(kpar)*invn
        sigma = (q2-q**2)*invn
        sigma = sqrt(max(sigma,0.0d0))
        write(out,'(i2,1x,1pe12.5,1x,1pe8.1)') kpar-2,q,2.0d0*sigma
      enddo

      ! Total (integrated) counter:
      q2 = tedep2*invn
      sigma = (q2-qtot**2)*invn
      sigma = sqrt(max(sigma,0.0d0))
      write(out,'(i2,1x,1pe12.5,1x,1pe8.1)') 3,qtot,2.0d0*sigma
      ! Evaluate rel. uncertainty:
      uncert = 200.0d0
      if (qtot.gt.0.0d0) uncert = 200.0d0*sigma/qtot

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
     & '#   Average uncertainty in % [uncert]:'
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


      subroutine PCSinitally
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*******************************************************************
      implicit none
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 countmp,count,count2,tcount2
      real*8 edeptmp,edep,edep2,tedep2,emin,ebin,iebin,unclimit
      common /scopcs/
     & countmp(3,maxbin),count(3,maxbin),count2(3,maxbin),tcount2(3),
     & edeptmp(3),edep(3),edep2(3),tedep2,
     & emin,ebin,iebin,unclimit,nbin,detmat,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PARTICLE CURRENT SPECTRUM v.2006-08-01]')
      parameter (eos='[END OF PCS SECTION]')
      character*80 buffer
      integer i,kpar,error
      real*8 emax

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PCSinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Particle Current Spectrum is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PCSinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PCSinitally:ERROR: expecting to find ON or OFF'
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
     &    'PCSinitally:ERROR: Too many bins; increase MAXBIN'
        stop
      endif
      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,1pe12.5)') unclimit

      ! Clear counters:
      tedep2   = 0.0d0
      do kpar=1,3
        edeptmp(kpar) = 0.0d0
        edep(kpar)    = 0.0d0
        edep2(kpar)   = 0.0d0
        tcount2(kpar)  = 0.0d0
        do i=1,nbin
          countmp(kpar,i) = 0.0d0
          count(kpar,i)   = 0.0d0
          count2(kpar,i)  = 0.0d0
        enddo
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PCSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PCS tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

