!*******************************************************************
!*                          TALLY                                  *
!*                     FLUENCE SPECTRUM                            *
!*          AND ENERGY DEPOSITION (TRACK LENGTH ESTIMATOR)         *
!*                                                                 *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Tallies track lengths in order to estimate the fluence spectr *
!*   (of each type of particle) in a given material.               *
!*                                                                 *
!*   Simultaneously, the energy deposition in the same material is *
!*   obtained by means of the same e+- track length estimator, i.e.*
!*   multiplying e+- track lengths by the stopping power restricted*
!*   below WCC,WCR. Direct energy depositions associated to photon *
!*   interactions, e+- track-ends and other "collision" events are *
!*   reported separately.                                          *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   -> common /CSIMPA/                                            *
!*   -> common /COMPOS/                                            *
!*   -> common /CEBR01/                                            *
!*   -> common /CJUMP0/                                            *
!*   -> routines EINaT,EBRaT,PINaT,PBRaT,SPLINE                    *
!*   -> EINF is set to the minimum energy PENELOPE can simulate    *
!*   from PENGEOM:                                                 *
!*   -> routine STEP                                               *
!*   from other penEasy files:                                     *
!*   -> routine GETLINE,FINDUF                                     *
!*   comments:                                                     *
!*   -> use of KNOCKX instead of KNOCK in MAIN is required.        *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2006                                                          *
!*******************************************************************


      subroutine FTLtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> identifies the state of the calling procedure      *
!*      arg -> different meanings for different modes; see MAIN    *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer maxmat
      real*8 eabs,c1,c2,wcc,wcr
      parameter (maxmat=10)
      common/csimpa/eabs(3,maxmat),c1(maxmat),c2(maxmat),wcc(maxmat),
     1  wcr(maxmat)
      real*8 p,st,dst,ds1,w1,w2,t1,t2
      common/cjump0/p(8),st,dst,ds1,w1,w2,t1,t2
      logical active
      integer maxbin,nbin,detmat
      parameter (maxbin=1000)
      real*8 flutmp,flu,flu2,energy,denerg,lemin,lebin,ilebin,eratio
      real*8 unclimit
      common /scoftl/ flutmp(3,maxbin),flu(3,maxbin),flu2(3,maxbin),
     &                energy(maxbin),denerg(maxbin),lemin,lebin,
     &                ilebin,eratio,unclimit,nbin,detmat,active
      logical called
      real*8 derest
      common /scofl0/ derest,called
      logical inside
      real*8 eini,ds,dssoft,desoft
      common /scofl1/ eini,ds,dssoft,desoft,inside
      real*8 erstmp,eres,eres2,edptmp,edep,edep2
      common /scofl2/ erstmp,eres,eres2,edptmp,edep,edep2
      integer i,k
      integer*4 ncross
      real*8 dseff,rstpow,deadd

      if (.not.active) return

      if (mode.eq.3) then
        ! STEP with no interface crossing, store JUMP distance for later use:
        ds = arg

      else if (mode.lt.0.and.mode.gt.-10) then  ! Modes reserved for KNOCK
        ! Check KNOCKX usage:
        if (.not.called) then
          write(*,*)
     &      'FTLtally:ERROR: fluence calculations require KNOCKX;'
          write(*,*)
     &      '  please substitute CALL KNOCK by CALL KNOCKX'
          write(*,*)
     &      '  in the main program'
          stop
        endif
        ! A KNOCK has occurred, check if inside detection material:
        if (.not.inside) return
        if (mode.eq.-1.and.kpar.ne.2) then
          ! Soft event; store partial jump distance and E-loss:
          dssoft = ds
          desoft = derest
          if (e.lt.1.0d0) then
            ! e+- absorption; correct restricted energy loss assuming
            ! the e+- travelled just until E=Eabs:
            derest = eini-eabs(kpar,mat)
            ! Scale down the step length just until E=Eabs
            ! (variable DS1 contains the remaining substep length--from JUMP):
            ds1 = (dssoft+ds1)*derest/desoft-dssoft
            if (ds1.gt.0.0d0) then
              call step(ds1,dseff,ncross)
            else
              ! Scaled-down step shorter than 1st substep -> 2nd substep <0:
              dseff = ds1
            endif
            ! Tally step track length and restricted energy loss:
            ds = (dssoft+dseff)*wght
            call trackl(kpar,eini,derest,ds)
            edptmp = edptmp+ds*rstpow(kpar,eini-0.5d0*derest)
          endif
        else
          ! Tally 2nd part of a completed soft step; or a jump for photons:
          ds = (dssoft+ds)*wght
          call trackl(kpar,eini,desoft,ds)
          if (kpar.ne.2)
     &      edptmp = edptmp+ds*rstpow(kpar,eini-0.5d0*desoft)
          ! Reset accumulated soft vars after a hard collision:
          eini = e
          dssoft = 0.0d0
          desoft = 0.0d0
        endif
        ! Tally residual energy loss (includes e+- track-ends):
        deadd = (arg-derest)*wght
        erstmp = erstmp+deadd
        edptmp = edptmp+deadd

      else if (mode.eq.4) then
        ! Interface crossing; treated as a hard interaction:
        if (inside) then
          ! Tally step track length and restricted energy loss:
          ds = (dssoft+arg)*wght
          call trackl(kpar,eini,desoft,ds)
          if (kpar.ne.2)
     &      edptmp = edptmp+ds*rstpow(kpar,eini-0.5d0*desoft)
          ! if we were inside the tally material, now we must be outside, so:
          inside = .false.
          return
        endif
        ! If we were not in the tally material, find out if we are now:
        if (mat.eq.detmat) then
          inside = .true.
          eini = e
          dssoft = 0.0d0
          desoft = 0.0d0
        endif

      else if (mode.eq.-99.or.mode.eq.0) then
        ! Stack was modified: a new particle retrieved (-99) or put by SOURCE (0):
        if (mat.eq.detmat) then
          inside = .true.
          eini = e
          dssoft = 0.0d0
          desoft = 0.0d0
          deadd = arg*wght
          erstmp = erstmp+deadd
          edptmp = edptmp+deadd
        else
          inside = .false.
        endif

      else if (mode.eq.6) then
        ! End-of-History; dump and clear counters:
        ! Fluence:
        do i=1,nbin
          do k=1,3
            flu(k,i)    = flu(k,i) +flutmp(k,i)
            flu2(k,i)   = flu2(k,i)+flutmp(k,i)**2
            flutmp(k,i) = 0.0d0
          enddo
        enddo
        ! Residual energy deposition:
        eres  = eres +erstmp
        eres2 = eres2+erstmp**2
        erstmp = 0.0d0
        ! Total energy deposition:
        edep   = edep+edptmp
        edep2  = edep2+edptmp**2
        edptmp = 0.0d0

      else if (mode.eq.-98.and.e.gt.eabs(kpar,mat).and.mat.eq.detmat)
     & then
        write(*,*)
     &    'FTLtally:ERROR: Particle absorbed inside the detector'
        write(*,*)
     &'             but its energy was higher than EABS; this'
        write(*,*)
     &'             prevents the calculation of the fluence down to'
        write(*,*)
     &'             EABS. Reformulate your particle rejection'
        write(*,*)
     &'             strategy appropriately.'
        stop

      endif
      end


      subroutine FTLreport(n,cputim,uncdone)
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
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer maxbin,nbin,detmat
      parameter (maxbin=1000)
      real*8 flutmp,flu,flu2,energy,denerg,lemin,lebin,ilebin,eratio
      real*8 unclimit
      common /scoftl/ flutmp(3,maxbin),flu(3,maxbin),flu2(3,maxbin),
     &                energy(maxbin),denerg(maxbin),lemin,lebin,
     &                ilebin,eratio,unclimit,nbin,detmat,active
      real*8 erstmp,eres,eres2,edptmp,edep,edep2
      common /scofl2/ erstmp,eres,eres2,edptmp,edep,edep2
      character*80 buffer
      integer i,nchan,out,finduf,error
      real*8 q,q2,sigma,eff,avesig,maxq,fact,emiddle,uncert,invn

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      !*******************************
      !*** Particle fluence report ***
      !*******************************

      do kpar=3,1,-1    ! One report for each particle type
        ! Prepare output files:
        out = finduf()
        if (kpar.eq.3)
     &  buffer = 'tallyFluenceTrackLength-positron.dat'
        if (kpar.eq.2)
     &  buffer = 'tallyFluenceTrackLength-photon.dat'
        if (kpar.eq.1)
     &  buffer = 'tallyFluenceTrackLength-electron.dat'
        open(out,file=buffer,iostat=error)
        if (error.ne.0) then
          write(*,*)
          write(*,'(a)')
     &      '*********************************************'
          write(*,'(a)')
     &      'FTLreport:ERROR: cannot open output data file'
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
          q = flu(kpar,i)/denerg(i)
          maxq = max(maxq,q)
        enddo
        maxq = 0.5d0*maxq

        write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)') '# [SECTION REPORT FLUENCE SPECTRUM]'
        if (kpar.eq.3) buffer = '# POSITRONS'
        if (kpar.eq.2) buffer = '# PHOTONS'
        if (kpar.eq.1) buffer = '# ELECTRONS'
        write(out,'(a)') buffer
        write(out,'(a)') '#'
        write(out,'(a)')
     &  '# The differential fluence is integrated over the detector'
        write(out,'(a)')
     &  '#   volume, hence its units are cm^3/(cm^2*eV)=cm/eV'
        write(out,'(a)')
     &  '# For plotting purposes, two energies per bin are given,'
        write(out,'(a)')
     &  '#   namely, the low end and the middle point of each bin'
        write(out,'(a)') '#'
        write(out,'(a)')
     &  '# Elow(eV) : Emiddle(eV) : dFluence(cm/eV) : +-2sigma'
        do i=1,nbin
          q = flu(kpar,i)*invn
          q2 = flu2(kpar,i)*invn
          sigma = (q2-q**2)*invn
          sigma = sqrt(max(sigma,0.0d0))
          fact = 1.0d0/denerg(i)
          q = q*fact
          sigma = sigma*fact
          emiddle = energy(i)*(eratio+1.0d0)*0.5d0  ! Mean energy in the bin
          write(out,'(3(1x,1pe12.5),1x,1pe8.1)')
     &      energy(i),emiddle,q,2.0d0*sigma
          ! Evaluate average uncertainty for scores above 1/2 max score:
          if (flu(kpar,i)*fact.gt.maxq) then
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
      enddo

      !*******************************
      !*** Energy deposited report ***
      !*******************************

      ! Prepare output files:
      out = finduf()
      open(out,file=
     &    'tallyFluenceTrackLength-totalE.dat',
     &     iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'FTLreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY DEPOSITION]'
      write(out,'(a)') '# (using e+- track length estimators)'
      write(out,'(a)') '#'
      write(out,'(a)')
     & '# Energy_deposition(eV) : +-2sigma : Residual(eV) : +-2sigma'
      q = edep*invn
      q2 = edep2*invn
      sigma = (q2-q**2)*invn
      sigma = sqrt(max(sigma,0.0d0))
      write(out,'(1pe12.5,1x,1pe8.1,$)') q,2.0d0*sigma
      uncert = 200.0d0
      if (q.gt.0.0d0) uncert = 200.0d0*sigma/q    ! Returned overall uncertainty
      q = eres*invn
      q2 = eres2*invn
      sigma = (q2-q**2)*invn
      sigma = sqrt(max(sigma,0.0d0))
      write(out,'(1x,1pe12.5,1x,1pe8.1)') q,2.0d0*sigma
      write(out,'(a)') ' '
      write(out,'(a)')
     &'#   Note: the residual contribution to the total energy'
      write(out,'(a)')
     &'#   deposition corresponds to track-ends and other events'
      write(out,'(a)')
     &'#   that contribute by mechanisms other than the product'
      write(out,'(a)')
     &'#    trackLength x restrictedStoppingPower.'

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


      subroutine FTLinitally
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Comments:                                                    *
!*      -> Must be called *AFTER* PEINIT                           *
!*******************************************************************
      implicit none
      integer maxmat
      parameter (maxmat=10)
      real*8 eabs,c1,c2,wcc,wcr
      common/csimpa/eabs(3,maxmat),c1(maxmat),c2(maxmat),wcc(maxmat),
     1  wcr(maxmat)
      logical active
      integer maxbin,nbin,detmat
      parameter (maxbin=1000)
      real*8 flutmp,flu,flu2,energy,denerg,lemin,lebin,ilebin,eratio
      real*8 unclimit
      common /scoftl/ flutmp(3,maxbin),flu(3,maxbin),flu2(3,maxbin),
     &                energy(maxbin),denerg(maxbin),lemin,lebin,
     &                ilebin,eratio,unclimit,nbin,detmat,active
      logical called
      real*8 derest
      common /scofl0/ derest,called
      real*8 erstmp,eres,eres2,edptmp,edep,edep2
      common /scofl2/ erstmp,eres,eres2,edptmp,edep,edep2
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY FLUENCE TRACK LENGTH v.2006-08-01]')
      parameter (eos='[END OF FTL SECTION]')
      character*80 buffer
      integer i,k,error
      real*8 emin,emax,lener,einf
      parameter (einf=50.0d0)  ! Minimum energy that PENELOPE can simulate

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'FTLinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Fluence Track Length is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'FTLinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'FTLinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Read data:
      write(*,'(a)') 'Detection material set to:'
      read(*,*) detmat
      write(*,'(i3)') detmat
      write(*,'(a)')
     &  'Spectrum limits [Emin,Emax] (eV) and no. bins:'
      read(*,*) emin,emax,nbin
      write(*,'(2(1x,1pe12.5),1x,i5)') emin,emax,nbin
      lemin = log(emin)
      lebin = (log(emax)-lemin)/nbin
      ilebin = 1.0d0/lebin
      eratio = (emax/emin)**(1.0d0/nbin)
      ! Add one bin to store, separately, counts with E=Emax:
      nbin = nbin+1
      if (nbin.gt.maxbin) then
        write(*,*)
     &    'FTLinitally:ERROR: Too many bins; increase MAXBIN'
        stop
      endif
      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,1pe12.5)') unclimit

      ! Clear fluence counters and set energy intervals width:
      do i=1,nbin
        do k=1,3
          flu(k,i)    = 0.0d0
          flu2(k,i)   = 0.0d0
          flutmp(k,i) = 0.0d0
        enddo
        lener = lemin+(i-1)*lebin
        energy(i) = exp(lener)
        denerg(i) = energy(i)*(eratio-1.0d0)
      enddo

      ! Clear residual energy counters:
      erstmp = 0.0d0
      eres   = 0.0d0
      eres2  = 0.0d0

      ! Clear energy deposition counters:
      edptmp = 0.0d0
      edep   = 0.0d0
      edep2  = 0.0d0

      ! Init restricted stopping power:
      call inirsp(max(eabs(1,detmat),wcc(detmat),einf),
     &            max(eabs(2,detmat),wcr(detmat),einf),detmat)

      ! Set KNOCKX test flag:
      called = .false.
      derest = 0.0d0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'FTLinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> FTL tally initialization finished >>>>'
      end


      subroutine trackl(kpar,eini,de,ds)
!*******************************************************************
!*    Tallies track lengths (from where fluences are calculated).  *
!*    The energy grid is logarithmic.                              *
!*                                                                 *
!*    Input:                                                       *
!*      kpar -> particle type                                      *
!*      eini -> initial particle energy (eV)                       *
!*      de -> energy loss along the step (eV)                      *
!*      ds -> track length (cm)                                    *
!*    Notes:                                                       *
!*      -> for an e+-, DE and DS refer to a whole hard-to-hard step*
!*******************************************************************
      implicit none
      integer*4 kpar
      real*8 eini,de,ds

      logical active
      integer maxbin,nbin,detmat
      parameter (maxbin=1000)
      real*8 flutmp,flu,flu2,energy,denerg,lemin,lebin,ilebin,eratio
      real*8 unclimit
      common /scoftl/ flutmp(3,maxbin),flu(3,maxbin),flu2(3,maxbin),
     &                energy(maxbin),denerg(maxbin),lemin,lebin,
     &                ilebin,eratio,unclimit,nbin,detmat,active
      integer i,binhi,binlo
      real*8 efin,istpav

      binhi = (log(eini)-lemin)*ilebin+1.0d0
      if (de.gt.0.0d0) then
        ! Average inverse stopping power:
        istpav = ds/de
      else
        ! No need to spread over various bins:
        if (binhi.lt.1.or.binhi.gt.nbin) return
        flutmp(kpar,binhi) = flutmp(kpar,binhi)+ds
        return
      endif

      efin = eini-de
      binlo = (log(efin)-lemin)*ilebin+1.0d0
      ! Spread over all bins from Efin to Eini, skipping invalid bins:
      do i=max(binlo,1),min(binhi,nbin)
        flutmp(kpar,i) = flutmp(kpar,i)+denerg(i)*istpav
      enddo

      ! Subtract excess in limiting bins, if needed:
      if (binhi.gt.0.and.binhi.le.nbin)
     &  flutmp(kpar,binhi) = flutmp(kpar,binhi)-
     &                       (energy(binhi)*eratio-eini)*istpav
      if (binlo.gt.0.and.binlo.le.nbin)
     &  flutmp(kpar,binlo) = flutmp(kpar,binlo)-
     &                       (efin-energy(binlo))*istpav
      end


      real*8 function rstpow(kpar,e)
!*******************************************************************
!*    Residual stopping power (spline interpolation). It is defined*
!*    similarly to the restricted stopping power, but the upper    *
!*    integration limit equals max(Wc,Eabs) instead of Wc.         *
!*                                                                 *
!*    Input:                                                       *
!*      kpar -> particle type (e- or e+)                           *
!*      e -> Kinetic energy in eV                                  *
!*    Output:                                                      *
!*      -> eV/cm                                                   *
!*    Comments:                                                    *
!*      -> It is assumed that kpar is never a photon.              *
!*******************************************************************
      implicit none
      integer*4 kpar
      real*8 e

      integer ndata
      parameter (ndata=1000)
      real*8 x,ye,ae,be,ce,de,yp,ap,bp,cp,dp,xmin,xmax,dx,idx
      common /stp001/
     &  x(ndata),xmin,xmax,dx,idx,
     &  ye(ndata),ae(ndata),be(ndata),ce(ndata),de(ndata),
     &  yp(ndata),ap(ndata),bp(ndata),cp(ndata),dp(ndata)
      integer i
      real*8 xdat

      xdat = log(e)
      i = (xdat-xmin)*idx+1.0d0
      if (kpar.eq.1) then
        rstpow = exp(ae(i)+xdat*(be(i)+xdat*(ce(i)+xdat*de(i))))
      else
        rstpow = exp(ap(i)+xdat*(bp(i)+xdat*(cp(i)+xdat*dp(i))))
      endif
      end


      subroutine inirsp(deltae,deltap,mat)
!*******************************************************************
!*    Initializes residual stopping power                          *
!*                                                                 *
!*    Input:                                                       *
!*      deltae -> cutoff for restricted delta-ray production       *
!*      deltap -> cutoff for restricted bremsstrahlung production  *
!*      mat -> material                                            *
!*    Comments:                                                    *
!*      -> Must be called *AFTER* PEINIT                           *
!*      -> Resorts to some 'deep' PENELOPE routines to obtain soft *
!*         stopping powers                                         *
!*      -> lines with 'c-dbg' can be activated to obtain restricted*
!*         stopping power tables of e+ and e- for debugging        *
!*         purposes                                                *
!*******************************************************************
      implicit none
      integer mat
      real*8 deltae,deltap

      ! PENELOPE 'deep' commons:
      integer*4 maxmat
      parameter (maxmat=10)
      integer*4 iz,nelem
      real*8 stf,zt,at,rho,vmol
      common/compos/stf(maxmat,30),zt(maxmat),at(maxmat),rho(maxmat),
     1  vmol(maxmat),iz(maxmat,30),nelem(maxmat)
      integer nbe,nbw,negp
      parameter (nbe=57,nbw=32,negp=200)
      real*8 ebt,xs,txs,xprima,yprima
      common/cebr01/ebt(nbe),xs(nbe,nbw),txs(nbe),
     &  xprima(nbe),yprima(nbe)
      integer ndata
      parameter (ndata=1000)
      real*8 x,ye,ae,be,ce,de,yp,ap,bp,cp,dp,xmin,xmax,dx,idx
      common /stp001/
     &  x(ndata),xmin,xmax,dx,idx,
     &  ye(ndata),ae(ndata),be(ndata),ce(ndata),de(ndata),
     &  yp(ndata),ap(ndata),bp(ndata),cp(ndata),dp(ndata)
      integer i
      integer*4 mlast
      real*8 xh0,xh1,xh2,xs0,xs1,xs2,xt1,xt2,delta
      real*8 spsin,spsbr,energy

      write(*,*) ' '
      write(*,'(a)')
     &  'inirsp: Restricted stopping powers initialization'
      write(*,'(a)')
     &  '  Cutoff energies (eV) for delta rays and bremss:'
      write(*,'(2(1x,1pe12.5))') deltae,deltap

      mlast = mat  ! Transform MAT to int*4:
      ! init 3-spline interpol:
      xmin = log(1.0d2)
      xmax = log(1.0d9)
      dx = (xmax-xmin)/(ndata-1)
      idx = 1.0d0/dx
      !-dbg write(*,*)
      !-dbg write(*,'(a)') 'init: Restricted Stopping Power (eV/cm)'
      !-dbg write(*,'(a)') 'E(eV) : RSP(-) : RBrem(-) : tot(-)'//
      !-dbg&           ': RSP(+) : RBrem(+) : tot(+)'
      do i=1,ndata
        x(i) = xmin+dx*(i-1)
        energy = exp(x(i))
        ! Electrons:
        call einat(energy,deltae,xh0,xh1,xh2,xs0,xs1,xs2,xt1,xt2,
     &             delta,mlast)
        !-unused sphin = xh1*vmol(mlast)
        spsin = xs1*vmol(mlast)
        call ebrat(energy,deltap,xh0,xh1,xh2,xs1,xs2,mlast)
        !-unused sphbr = xh1*vmol(mlast)
        spsbr = xs1*vmol(mlast)
        !-dbg   write(*,'(4(1x,1pe12.5),$)') energy,spsin,spsbr,spsin+spsbr
        ye(i) = log(spsin+spsbr)
        ! Positrons:
        call pinat(energy,deltae,xh0,xh1,xh2,xs0,xs1,xs2,xt1,xt2,
     &             delta,mlast)
        !-unused sphin = xh1*vmol(mlast)
        spsin = xs1*vmol(mlast)
        call pbrat(energy,deltap,xh0,xh1,xh2,xs1,xs2,mlast)
        !-unused sphbr = xh1*vmol(mlast)
        spsbr = xs1*vmol(mlast)
        !-dbg   write(*,'(3(1x,1pe12.5))') spsin,spsbr,spsin+spsbr
        yp(i) = log(spsin+spsbr)
      enddo
      call spline(x,ye,ae,be,ce,de,0.0d0,0.0d0,ndata)
      call spline(x,yp,ap,bp,cp,dp,0.0d0,0.0d0,ndata)
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

