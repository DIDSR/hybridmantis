!*******************************************************************
!*                          TALLY                                  *
!*            CYLINDRICAL (rho-z) DOSE DISTRIBUTION                *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Dose distribution in a grid of (rho,z) volumetric bins,       *
!*   where rho^2 is (x^2+y^2). The grid of bins is superimposed on *
!*   the solid-body PENGEOM geometry.                              *
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


      subroutine CDDtally(mode,arg)
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
      integer*4 binmax
      parameter(binmax=100000)
      integer*4 nr,nz
      real*8 edptmp,edep,edep2,idens,rmin,zmin
      real*8 dr,dz,idr,idz,nlast,nhist,unclimit
      common /scocdd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),rmin,zmin,dr,dz,
     &                idr,idz,unclimit,nhist,nr,nz,prtxyz,active
      integer*4 bin,i,k
      real*8 oneplus
      parameter(oneplus=1.0000000000001d0)

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0d0) return  ! Nothing to deposit
        ! Check if particle is inside tally region:
        i = (sqrt(x**2+y**2)-rmin)*idr+1.0d0
        if (i.lt.1.or.i.gt.nr) return
        k = (z-zmin)*idz+oneplus  ! ONEPLUS ensures proper truncation when IDZ=0
        if (k.lt.1.or.k.gt.nz) return
        bin = k+(i-1)*nz
        ! Transfer partial score to totals only when a new history visits:
        if (nhist.gt.nlast(bin)) then
          edep(bin)  = edep(bin) +edptmp(bin)
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= arg*wght
          nlast(bin) = nhist+0.5d0
        else
          edptmp(bin) = edptmp(bin)+arg*wght
        endif

      else if (mode.eq.1.or.mode.eq.2) then  ! New history or hist. modified
        nhist = arg

      endif
      end


      subroutine CDDreport(n,cputim,uncdone)
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
      integer*4 binmax
      parameter(binmax=100000)
      integer*4 nr,nz
      real*8 edptmp,edep,edep2,idens,rmin,zmin
      real*8 dr,dz,idr,idz,nlast,nhist,unclimit
      common /scocdd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),rmin,zmin,dr,dz,
     &                idr,idz,unclimit,nhist,nr,nz,prtxyz,active
      integer nchan,out,finduf,error
      integer*4 bin,i,k,nzef
      real*8 q,sigma,eff,avesig,maxq,fact,r,z,uncert,zmiddle,rave
      real*8 twothird,invdr2,pi,invpi,invn,factidz
      parameter (twothird=2.0d0/3.0d0)
      parameter (pi=3.1415926535897932d0,invpi=1.0d0/pi)

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyCylindricalDoseDistrib.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'CDDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Dump counters and obtain max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do bin=1,nr*nz
        if (nlast(bin).gt.0.5d0) then
          edep(bin)  = edep(bin) +edptmp(bin)
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= 0.0d0
          nlast(bin) = 0.0d0
        endif
        maxq = max(maxq,edep(bin))    ! 1/2 of the max score
      enddo
      maxq = 0.5d0*maxq

      ! Prepare z factors:
      if (idz.eq.0.0d0) then
        nzef = 0
        factidz = 1.0d0
      else
        nzef = nz
        factidz = idz
      endif

      ! Write header:
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT CYLINDRICAL DOSE DISTRIB]'
      if (nzef.ne.0) then
        write(out,'(a)') '# Dose units are eV/g per history'
      else
        write(out,'(a)') '# Dose units are eV.cm/g per history'
      endif
      write(out,'(a)') '# No. of radial (r) and depth (z) bins:'
      write(out,'(a,2(1x,i0))') '# ',nr,nzef
      write(out,'(a)') '# Min value and bin width for r and z(cm):'
      write(out,'(a,2(2x,1pe12.5,1x,1pe12.5))')
     &  '# ',rmin,dr,zmin,dz
      if (prtxyz.eq.1) then
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# For plotting purposes, two values per bin '//
     &    'coordinate are given, namely,'
        write(out,'(a)')
     &    '#   the low end and an average value.'
        write(out,'(a)')
     &    '# For the z coordinate the average is the middle'//
     &    ' point of the bin.'
        write(out,'(a)')
     &    '# For the r coordinate it is a weighted average'//
     &    ' with a weight proportional'
        write(out,'(a)') '#   to the radius r.'
      endif
      write(out,'(a)') '#'
      write(out,'(a,$)') '# '
      if (prtxyz.eq.1) then
        write(out,'(a,$)') 'rBinIndex : rLow(cm) : rAve(cm) : '
        if (idz.gt.0.0d0)
     &    write(out,'(a,$)') 'zBinIndex : zLow(cm) : zMiddle(cm) : '
      endif
      write(out,'(a)') 'dose : +-2sigma'

      ! Write data:
      do k=1,nz
        z = zmin+dz*(k-1)
        zmiddle = z+dz*0.5d0
        if (idz.gt.0.0d0.and.prtxyz.ne.1)  ! Since z is not written, give at least a summary
     &    write(out,'(a,i0,a,1pe12.5)')
     &    '# zBinIndex=',k,' zMiddle(cm)=',zmiddle

        do i=1,nr
          r = rmin+dr*(i-1)
          invdr2 = 1.0d0/(dr*(dr+2.0d0*r))  ! Note: dr2 = (r+dr)^2-r^2
          if (prtxyz.eq.1) then
            if (r.gt.0.0d0) then
              rave = twothird*((r+dr)**3-r**3)*invdr2  ! Average with weight(r)~r
            else
              rave = 0.0d0  ! If first bin is a circle, the average r is zero
            endif
            write(out,'(1x,i5,2(1x,1pe12.5),$)') i,r,rave
          endif
          if (idz.gt.0.0d0.and.prtxyz.eq.1)
     &        write(out,'(1x,i5,2(1x,1pe12.5),$)') k,z,zmiddle
          bin = k+(i-1)*nz
          fact = factidz*idens(bin)*invpi*invdr2
          q = edep(bin)*invn
          sigma = sqrt(max((edep2(bin)*invn-q**2)*invn,0.0d0))*fact
          q = q*fact
          write(out,'(1x,1pe12.5,1x,1pe8.1)') q,2.0d0*sigma

          ! Evaluate average uncertainty for scores above 1/2 max score:
          if (edep(bin).gt.maxq) then
            avesig = avesig+(sigma/q)**2
            nchan = nchan+1
          endif
        enddo

        if (nr.gt.1) write(out,*) ' '  ! Separate data blocks
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


      subroutine CDDinitally
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
      integer*4 binmax
      parameter(binmax=100000)
      integer*4 nr,nz
      real*8 edptmp,edep,edep2,idens,rmin,zmin
      real*8 dr,dz,idr,idz,nlast,nhist,unclimit
      common /scocdd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),rmin,zmin,dr,dz,
     &                idr,idz,unclimit,nhist,nr,nz,prtxyz,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY CYLINDRICAL DOSE DISTRIB v.2006-08-01]')
      parameter (eos='[END OF CDD SECTION]')
      character*80 buffer
      integer error
      integer*4 i,k,bin
      real*8 rmax,zmax,dens,nobins

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'CDDinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Cylindrical Dose Distrib is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'CDDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'CDDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Read geometry parameters:
      write(*,'(a)') 'rmin,rmax,nr:'
      read(*,*) rmin,rmax,nr
      write(*,'(2(1x,1pe12.5),1x,i5)') rmin,rmax,nr
      if (rmin.lt.0.0d0.or.rmin.ge.rmax.or.nr.lt.1) then
        write(*,'(a,a)') 'CDDinitally:ERROR: the condition ',
     &    '0 <= rmin < rmax  and  nr > 0 is not satisfied'
        stop
      endif
      dr = (rmax-rmin)/nr
      idr = 1.0d0/dr
      write(*,'(a)') 'zmin,zmax,nz:'
      read(*,*) zmin,zmax,nz
      write(*,'(2(1x,1pe12.5),1x,i5)') zmin,zmax,nz
      if (zmin.gt.zmax.or.nz.lt.0.or.(nz.gt.0.and.zmin.eq.zmax)) then
        write(*,*) 'CDDinitally:ERROR: zmin >= zmax  or  nz < 0'
        stop
      endif
      if (nz.eq.0) then
        dz = zmax-zmin
        idz = 0.0d0
        nz = 1
      else
        dz = (zmax-zmin)/nz
        idz = 1.0d0/dz
      endif
      nobins = dble(nr)*dble(nz)
      if (nobins.gt.binmax) then
        write(*,*)
     &    'CDDinitally:ERROR: Too many bins; increase binmax.'
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
        x = rmin+dr*(i-0.5d0)  ! This is to locate a point and find its material
        y = 0.0d0
        do k=1,nz
          bin = k+(i-1)*nz
          edptmp(bin) = 0.0d0
          edep(bin)   = 0.0d0
          edep2(bin)  = 0.0d0
          nlast(bin)  = 0.0d0
          z = zmin+dz*(k-0.5d0)
          u = 0.0d0
          v = 0.0d0
          w = 1.0d0
          call locate
          dens = 0.0d0
          idens(bin) = 0.0d0
          if (mat.ne.0) dens = rho(mat)
          if (dens.gt.0.0d0) idens(bin) = 1.0d0/dens
        enddo
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'CDDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> CDD tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
