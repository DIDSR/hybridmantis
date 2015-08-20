!*******************************************************************
!*                          TALLY                                  *
!*                 SPHERICAL DOSE DISTRIBUTION                     *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Dose distribution in a set of spherical shells. This set is   *
!*   superimposed on the solid-body PENGEOM geometry.              *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   -> common /COMPOS/                                            *
!*   from PENGEOM:                                                 *
!*   -> routine LOCATE                                             *
!*   from other penEasy files:                                     *
!*   -> routine GETLINE,FINDUF                                     *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*******************************************************************


      subroutine SPDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer prtxyz
      integer*4 binmax,nr
      parameter(binmax=1000)
      real*8 edptmp,edep,edep2,idens,rmin,dr,idr,nlast,nhist,unclimit
      common /scospd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),rmin,dr,idr,
     &                unclimit,nhist,nr,prtxyz,active
      integer*4 i

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0d0) return    ! Nothing to deposit
        ! Check if particle is inside tally region:
        i = (sqrt(x**2+y**2+z**2)-rmin)*idr+1.0d0
        if (i.lt.1.or.i.gt.nr) return
        ! Transfer partial tally to totals only when a new history visits:
        if (nhist.gt.nlast(i)) then
          edep(i)  = edep(i) +edptmp(i)
          edep2(i) = edep2(i)+edptmp(i)**2
          edptmp(i)= arg*wght
          nlast(i) = nhist+0.5d0
        else
          edptmp(i) = edptmp(i)+arg*wght
        endif

      else if (mode.eq.1.or.mode.eq.2) then  ! New history or hist. modified
        nhist = arg

      endif
      end


      subroutine SPDreport(n,cputim,uncdone)
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
      integer prtxyz
      integer*4 binmax,nr
      parameter(binmax=1000)
      real*8 edptmp,edep,edep2,idens,rmin,dr,idr,nlast,nhist,unclimit
      common /scospd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),rmin,dr,idr,
     &                unclimit,nhist,nr,prtxyz,active
      integer nchan,out,finduf,error
      integer*4 i
      real*8 q,sigma,eff,avesig,maxq,fact,uncert,r,rave,rnext
      real*8 threefour,invdr3,pi,inv43pi,invn
      parameter (threefour=3.0d0/4.0d0)
      parameter (pi=3.1415926535897932d0,inv43pi=3.0d0/(4.0d0*pi))

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      ! Prepare output files:
      out = finduf()
      open(out,file='tallySphericalDoseDistrib.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'SPDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Dump counters and obtain max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do i=1,nr
        if (nlast(i).gt.0.5d0) then
          edep(i)  = edep(i) +edptmp(i)
          edep2(i) = edep2(i)+edptmp(i)**2
          edptmp(i)= 0.0d0
          nlast(i) = 0.0d0
        endif
        maxq = max(maxq,edep(i))    ! 1/2 of the max score
      enddo
      maxq = 0.5d0*maxq

      ! Write header:
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT SPHERICAL DOSE DISTRIB]'
      write(out,'(a)') '# Dose units are: eV/g per history'
      write(out,'(a)') '# No. of radial (r) bins:'
      write(out,'(a,1x,i0)') '# ',nr
      write(out,'(a)') '# Min r and bin width (cm):'
      write(out,'(a,2(1x,1pe12.5))') '# ',rmin,dr
      if (prtxyz.eq.1) then
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# For plotting purposes, two radii r for each bin'//
     &    ' are given, namely, the low end'
        write(out,'(a)')
     &    '#   and an average radius; the latter is'//
     &    ' a weighted average with a weight'
        write(out,'(a)') '#   proportional to r squared'
      endif
      write(out,'(a)') '#'
      write(out,'(a,$)') '# '
      if (prtxyz.eq.1) then
        write(out,'(a,$)') 'rBinIndex : rLow(cm) : rAve(cm) : '
      endif
      write(out,'(a)') 'dose : +-2sigma'

      ! Write data:
      do i=1,nr
        r = rmin+dr*(i-1)
        rnext = r+dr
        invdr3 = 1.0d0/(rnext**3-r**3)
        fact = idens(i)*inv43pi*invdr3  ! This is 1/Delta_mass
        if (prtxyz.eq.1) then
          if (r.gt.0.0d0) then
            rave = threefour*(rnext**4-r**4)*invdr3  ! Average with weight(r)~r^2
          else
            rave = 0.0d0  ! If first bin is a sphere, the average r is zero
          endif
          write(out,'(1x,i5,2(1x,1pe12.5),$)') i,r,rave
        endif
        q = edep(i)*invn
        sigma = sqrt(max((edep2(i)*invn-q**2)*invn,0.0d0))*fact
        q = q*fact
        write(out,'(1x,1pe12.5,1x,1pe8.1)') q,2.0d0*sigma
        ! Evaluate average uncertainty for scores above 1/2 max score:
        if (edep(i).gt.maxq) then
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


      subroutine SPDinitally
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer maxmat
      parameter (maxmat=10)
      integer*4 iz,nelem
      real*8 stf,zt,at,rho,vmol
      common/compos/stf(maxmat,30),zt(maxmat),at(maxmat),rho(maxmat),
     1  vmol(maxmat),iz(maxmat,30),nelem(maxmat)
      logical active
      integer prtxyz
      integer*4 binmax,nr
      parameter(binmax=1000)
      real*8 edptmp,edep,edep2,idens,rmin,dr,idr,nlast,nhist,unclimit
      common /scospd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),rmin,dr,idr,
     &                unclimit,nhist,nr,prtxyz,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY SPHERICAL DOSE DISTRIB v.2006-08-01]')
      parameter (eos='[END OF SPD SECTION]')
      character*80 buffer
      integer error
      integer*4 i
      real*8 rmax,dens

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'SPDinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Spherical Dose Distrib is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'SPDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'SPDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Read geometry parameters:
      write(*,'(a)') 'rmin,rmax,nr:'
      read(*,*) rmin,rmax,nr
      write(*,'(2(1x,1pe12.5),1x,i5)') rmin,rmax,nr
      if (rmin.lt.0.0d0.or.rmin.ge.rmax.or.nr.lt.1) then
        write(*,'(a,a)') 'SPDinitally:ERROR: the condition ',
     &    '0 <= rmin < rmax  and  nr > 0 is not satisfied'
        stop
      endif
      dr = (rmax-rmin)/nr
      idr = 1.0d0/dr
      if (nr.gt.binmax) then
        write(*,*)
     &    'SPDinitally:ERROR: Too many bins; increase binmax.'
        stop
      endif

      write(*,'(a)') 'Print coordinates:'
      read(*,*) prtxyz
      if (prtxyz.eq.1) then
        write(*,'(a)') ' yes'
      else
        write(*,'(a)') ' no'
      endif

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,1pe12.5)') unclimit

      ! Init arrays:
      do i=1,nr
        x = 0.0d0
        y = 0.0d0
        z = rmin+dr*(i-0.5d0)  ! This is to locate a point and find its material
        edptmp(i) = 0.0d0
        edep(i)   = 0.0d0
        edep2(i)  = 0.0d0
        nlast(i)  = 0.0d0
        u = 0.0d0
        v = 0.0d0
        w = 1.0d0
        call locate
        dens = 0.0d0
        idens(i) = 0.0d0
        if (mat.ne.0) dens = rho(mat)
        if (dens.gt.0.0d0) idens(i) = 1.0d0/dens
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'SPDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> SPD tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
