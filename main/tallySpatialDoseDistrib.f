!*******************************************************************
!*                          TALLY                                  *
!*             SPATIAL (3D) DOSE DISTRIBUTION                      *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Dose distribution in a grid of volumetric bins                *
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


      subroutine SDDtally(mode,arg)
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
      parameter(binmax=1000000)
      integer*4 nx,ny,nz,nynz
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz,nlast,nhist,unclimit
      common /scosdd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,unclimit,nhist,
     &                nx,ny,nz,nynz,prtxyz,active
      integer*4 bin,i,j,k
      real*8 oneplus
      parameter(oneplus=1.0000000000001d0)

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0d0) return       ! Nothing to deposit
        ! Check if particle is inside tally region:
        i = (x-xmin)*idx+oneplus       ! ONEPLUS ensures proper truncation when IDX=0
        if (i.lt.1.or.i.gt.nx) return
        j = (y-ymin)*idy+oneplus
        if (j.lt.1.or.j.gt.ny) return
        k = (z-zmin)*idz+oneplus
        if (k.lt.1.or.k.gt.nz) return
        bin = k+(j-1)*nz+(i-1)*nynz    ! Map i,j,k into a single index
        ! Transfer partial tally to totals only when a new history visits:
        if (nhist.gt.nlast(bin)) then  ! Visit of a new history
          edep(bin)  = edep(bin) +edptmp(bin)
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= arg*wght
          nlast(bin) = nhist+0.5d0     ! Avoid round-off problems
        else
          edptmp(bin) = edptmp(bin)+arg*wght
        endif

      else if (mode.eq.1.or.mode.eq.2) then  ! New history or hist. modified
        nhist = arg

      endif
      end


      subroutine SDDreport(n,cputim,uncdone)
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
      parameter(binmax=1000000)
      integer*4 nx,ny,nz,nynz
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz,nlast,nhist,unclimit
      common /scosdd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,unclimit,nhist,
     &                nx,ny,nz,nynz,prtxyz,active
      logical printx,printy,printz
      integer nchan,out,finduf,error,ninteg
      integer*4 bin,i,j,k,nxef,nyef,nzef
      real*8 q,sigma,eff,avesig,maxq,fact,factint,x,y,z,uncert
      real*8 xmiddle,ymiddle,zmiddle,invn

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      ! Prepare output files:
      out = finduf()
      open(out,file='tallySpatialDoseDistrib.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'SDDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Dump counters and obtain max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do bin=1,nx*ny*nz
        if (nlast(bin).gt.0.5d0) then  ! At least visited once
          edep(bin)  = edep(bin) +edptmp(bin)
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= 0.0d0
          nlast(bin) = 0.0d0
        endif
        maxq = max(maxq,edep(bin))     ! 1/2 of the max score
      enddo
      maxq = 0.5d0*maxq

      ! Prepare x,y,z factors:
      ninteg = 0                       ! No. of dimensions to integrate
      factint = 1.0d0
      printx = .false.
      printy = .false.
      printz = .false.
      nxef = nx
      if (idx.eq.0.0d0) then
        nxef = 0
        ninteg = ninteg+1
      else
        factint = factint*idx
        if (prtxyz.eq.1) printx = .true.
      endif
      nyef = ny
      if (idy.eq.0.0d0) then
        nyef = 0
        ninteg = ninteg+1
      else
        factint = factint*idy
        if (prtxyz.eq.1) printy = .true.
      endif
      nzef = nz
      if (idz.eq.0.0d0) then
        nzef = 0
        ninteg = ninteg+1
      else
        factint = factint*idz
        if (prtxyz.eq.1) printz = .true.
      endif

      ! Write header:
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT SPATIAL DOSE DISTRIB]'
      if (ninteg.eq.0) write(out,'(a)')
     &  '# Dose units are: eV/g per history'
      if (ninteg.eq.1) write(out,'(a)')
     &  '# Dose units are: eV.cm/g per history'
      if (ninteg.eq.2) write(out,'(a)')
     &  '# Dose units are: eV.cm^2/g per history'
      if (ninteg.eq.3) write(out,'(a)')
     &  '# Dose units are: eV.cm^3/g per history'
      write(out,'(a)') '# No. of bins in x,y,z directions:'
      write(out,'(a,3(1x,i0))') '# ',nxef,nyef,nzef
      write(out,'(a)') '# Min values and bin widths for x,y,z(cm):'
      write(out,'(a,3(2x,1pe12.5,1x,1pe12.5))')
     &  '# ',xmin,dx,ymin,dy,zmin,dz
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# For plotting purposes, two values per bin '//
     &  'coordinate are given, namely,'
      write(out,'(a)')
     &  '#   the low end and the middle point of each bin'
      write(out,'(a)') '#'
      write(out,'(a,$)') '# '
      if (printx)
     &    write(out,'(a,$)') 'xBinIndex : xLow(cm) : xMiddle(cm) : '
      if (printy)
     &    write(out,'(a,$)') 'yBinIndex : yLow(cm) : yMiddle(cm) : '
      if (printz)
     &    write(out,'(a,$)') 'zBinIndex : zLow(cm) : zMiddle(cm) : '
      write(out,'(a)') 'dose : +-2sigma'

      ! Write data:
      do k=1,nz
        z = zmin+dz*(k-1)
        zmiddle = z+dz*0.5d0
        if (idz.gt.0.0d0.and..not.printz)  ! Since z is not written, give at least a summary
     &    write(out,'(a,i0,a,1pe12.5)')
     &    '# zBinIndex=',k,' zMiddle(cm)=',zmiddle

        do j=1,ny
          y = ymin+dy*(j-1)
          ymiddle = y+dy*0.5d0
          if (idy.gt.0.0d0.and..not.printy)  ! Since y is not written, give at least a summary
     &      write(out,'(a,i0,a,1pe12.5)')
     &      '# yBinIndex=',j,' yMiddle(cm)=',ymiddle

          do i=1,nx
            x = xmin+dx*(i-1)
            xmiddle = x+dx*0.5d0
            if (printx)
     &        write(out,'(1x,i5,2(1x,1pe12.5),$)') i,x,xmiddle
            if (printy)
     &        write(out,'(1x,i5,2(1x,1pe12.5),$)') j,y,ymiddle
            if (printz)
     &        write(out,'(1x,i5,2(1x,1pe12.5),$)') k,z,zmiddle
            bin = k+(j-1)*nz+(i-1)*nynz  ! Map i,j,k into a single index
            fact = factint*idens(bin)
            q = edep(bin)*invn
            sigma = sqrt(max((edep2(bin)*invn-q**2)*invn,0.0d0))*fact
            q = q*fact
            write(out,'(1x,1pe12.5,1x,1pe7.1)') q,2.0d0*sigma
            ! Evaluate average uncertainty for scores above 1/2 max score:
            if (edep(bin).gt.maxq) then
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
          enddo

          if (nx.gt.1) write(out,*) ' '  ! Separate 2D data blocks
        enddo

        if (ny.gt.1) write(out,*) ' '  ! Separate 3D data blocks
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


      subroutine SDDinitally
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
      parameter(binmax=1000000)
      integer*4 nx,ny,nz,nynz
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz,nlast,nhist,unclimit
      common /scosdd/ edptmp(binmax),edep(binmax),edep2(binmax),
     &                idens(binmax),nlast(binmax),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,unclimit,nhist,
     &                nx,ny,nz,nynz,prtxyz,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY SPATIAL DOSE DISTRIB v.2006-08-01]')
      parameter (eos='[END OF SDD SECTION]')
      character*80 buffer
      integer error
      integer*4 i,j,k,bin
      real*8 xmax,ymax,zmax,dens,nobins

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'SDDinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Spatial Dose Distrib is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'SDDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'SDDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Read geometry parameters:
      write(*,'(a)') 'xmin,xmax,nx:'
      read(*,*) xmin,xmax,nx
      write(*,'(2(1x,1pe12.5),1x,i5)') xmin,xmax,nx
      if (xmin.gt.xmax.or.nx.lt.0.or.(nx.gt.0.and.xmin.eq.xmax)) then
        write(*,*) 'SDDinitally:ERROR: xmin >= xmax  or  nx < 0'
        stop
      endif
      if (nx.eq.0) then
        ! Redefine parameters when dose is integrated along x-axis:
        dx = xmax-xmin
        idx = 0.0d0
        nx = 1
      else
        dx = (xmax-xmin)/nx
        idx = 1.0d0/dx
      endif
      write(*,'(a)') 'ymin,ymax,ny:'
      read(*,*) ymin,ymax,ny
      write(*,'(2(1x,1pe12.5),1x,i5)') ymin,ymax,ny
      if (ymin.gt.ymax.or.ny.lt.0.or.(ny.gt.0.and.ymin.eq.ymax)) then
        write(*,*) 'SDDinitally:ERROR: ymin >= ymax  or  ny < 0'
        stop
      endif
      if (ny.eq.0) then
        dy = ymax-ymin
        idy = 0.0d0
        ny = 1
      else
        dy = (ymax-ymin)/ny
        idy = 1.0d0/dy
      endif
      write(*,'(a)') 'zmin,zmax,nz:'
      read(*,*) zmin,zmax,nz
      write(*,'(2(1x,1pe12.5),1x,i5)') zmin,zmax,nz
      if (zmin.gt.zmax.or.nz.lt.0.or.(nz.gt.0.and.zmin.eq.zmax)) then
        write(*,*) 'SDDinitally:ERROR: zmin >= zmax  or  nz < 0'
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
      nobins = dble(nx)*dble(ny)*dble(nz)
      if (nobins.gt.binmax) then
        write(*,*)
     &    'SDDinitally:ERROR: Too many bins; increase binmax.'
        stop
      endif
      nynz = ny*nz

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
      do i=1,nx
        x = xmin+dx*(i-0.5d0)
        do j=1,ny
          y = ymin+dy*(j-0.5d0)
          do k=1,nz
            bin = k+(j-1)*nz+(i-1)*nynz  ! Map i,j,k into a single index
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
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'SDDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> SDD tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
