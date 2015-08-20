!*******************************************************************
!*               TALLY VOXEL DOSE DISTRIBUTION                     *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Dose distribution in the volume elements (voxels) defined in  *
!*   a PENVOX geometry file.                                       *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/, /COMPOS/, /RSEED/                          *
!*   from other penEasy files:                                     *
!*   -> common /GEOQUAD/                                           *
!*   -> routine FINDUF                                             *
!*   from PENVOX:                                                  *
!*   -> commons /GEOVOX/,/PARTVOX/                                 *
!*   -> routines LOCATEVOX                                         *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*                                                                 *
!* Authors:                                                        *
!*     Andreu Badal & Josep Sempau                                 *
!*     Universitat Politecnica de Catalunya, Barcelona, Spain      *
!*                                                                 *
!* SEE COPYRIGHT NOTICE IN README.txt under folder penEasy_Imaging *
!*******************************************************************


      subroutine VDDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*    Comments: ** IMPORTANT**                                     *
!*      -> Since this tally modifies the particle state variables  *
!*         through LOCATEVOX, it is required to be the 1st tally   *
!*         to be called by the main program.                       *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      character*80 voxfilen
      logical isFullvol
      integer matvox,granul,maxGranul
      parameter (maxGranul=1000)
      integer*4 nx,ny,nz,nxy,maxvox,bodymask,matmask
      parameter (maxvox=10000000)  ! Max no. of voxels
      real densvox,idensvox
      real*8 massvox,dx,dy,dz,idx,idy,idz,vbb
      common /geovox/ massvox(maxvox),densvox(maxvox),idensvox(maxvox),
     &  matvox(maxvox),dx,dy,dz,idx,idy,idz,vbb(3),nx,ny,nz,nxy,
     &  bodymask,matmask,granul,isFullvol,voxfilen
      integer*4 xvox,yvox,zvox,absvox
      real*8 uold,vold,ivx,ivy,ivz,sdx,sdy,sdz
      common /partvox/ uold,vold,ivx,ivy,ivz,sdx,sdy,sdz,
     &                 xvox,yvox,zvox,absvox
      logical active
      integer prtxyz,prtdens
      integer*4 xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax
      real*8 nlast,nhist,unclimit,edptmp,edep,edep2
      common /scovdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                nlast(maxvox),unclimit,nhist,
     &                xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax,
     &                prtxyz,prtdens,active
      integer*4 vox

      if (mode.eq.-99.and.ilb(1).ne.1) then      ! Locate 2nd particle
        if (ibody.eq.bodymask) then
          call locatevox                         ! Check voxels
        else                                     ! Inside quadrics
          absvox = 0                             ! Mark as quadrics
        endif
      endif

      if (.not.active) return                    ! Voxel doses not tallied

      if (mode.le.0) then                        ! Deposit energy
        if (arg.eq.0.0d0) return                 ! Nothing to deposit

        if (absvox.eq.0) then                    ! Inside quadrics
          if (isFullvol) then                    ! Consider partial volume
            xvox = x*idx+1.0d0                   ! Compute indices
            if(xvox.lt.xvoxmin.or.xvox.gt.xvoxmax) return ! Not in ROI
            yvox = y*idy+1.0d0
            if(yvox.lt.yvoxmin.or.yvox.gt.yvoxmax) return ! Not in ROI
            zvox = z*idz+1.0d0
            if(zvox.lt.zvoxmin.or.zvox.gt.zvoxmax) return ! Not in ROI
            vox = xvox+(yvox-1)*nx+(zvox-1)*nxy  ! Absolute voxel index
          else                                   ! Ignore partial volume
            return                               ! Nothing to do
          endif
        else                                     ! Inside voxels
          vox = absvox                           ! Transfer abs voxel
        endif

        if (nhist.gt.nlast(vox)) then            ! Visit of a new history
          edep(vox)  = edep(vox) +edptmp(vox)    ! Transfer partial scoring to totals
          edep2(vox) = edep2(vox)+edptmp(vox)**2 ! Score squared values for variance
          edptmp(vox)= arg*wght                  ! And update temporary counter
          nlast(vox) = nhist+0.5d0               ! Update NLAST (+0.5 avoids prec. problems)
        else
          edptmp(vox) = edptmp(vox)+arg*wght     ! Same history: update temporary counter only
        endif

      else if (mode.eq.1.or.mode.eq.2) then      ! New history or hist. modified
        nhist = arg                              ! Update NHIST

      endif
      end


      subroutine VDDreport(n,cputim,uncdone)
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
      character*80 voxfilen
      logical isFullvol
      integer matvox,granul,maxGranul
      parameter (maxGranul=1000)
      integer*4 nx,ny,nz,nxy,maxvox,bodymask,matmask
      parameter (maxvox=10000000)  ! Max no. of voxels
      real densvox,idensvox
      real*8 massvox,dx,dy,dz,idx,idy,idz,vbb
      common /geovox/ massvox(maxvox),densvox(maxvox),idensvox(maxvox),
     &  matvox(maxvox),dx,dy,dz,idx,idy,idz,vbb(3),nx,ny,nz,nxy,
     &  bodymask,matmask,granul,isFullvol,voxfilen
      logical active
      integer prtxyz,prtdens
      integer*4 xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax
      real*8 nlast,nhist,unclimit,edptmp,edep,edep2
      common /scovdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                nlast(maxvox),unclimit,nhist,
     &                xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax,
     &                prtxyz,prtdens,active
      logical isQuad
      common /geoquad/ isQuad
      character typevox(3)
      integer nchan,out,finduf,error
      integer*4 vox,i,j,k
      real*8 q,sigma,eff,avesig,maxq,fact,x,y,z,uncert
      real*8 xmiddle,ymiddle,zmiddle,invn

      uncdone = 1
      if (.not.active) return
      invn = 1.0d0/n

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyVoxelDoseDistrib.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'VDDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Dump counters and obtain max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do k=zvoxmin,zvoxmax
        do j=yvoxmin,yvoxmax
          do i=xvoxmin,xvoxmax
            vox = i+(j-1)*nx+(k-1)*nxy     ! Absolute voxel index
            if (nlast(vox).gt.0.5d0) then  ! Do not update if nothing scored
              edep(vox)  = edep(vox) +edptmp(vox)
              edep2(vox) = edep2(vox)+edptmp(vox)**2
              edptmp(vox)= 0.0d0
              nlast(vox) = 0.0d0
            endif
            maxq = max(maxq,edep(vox))
          end do
        end do
      end do
      maxq = 0.5d0*maxq     ! 1/2 of the max score

      ! Write header:
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT VOXEL DOSE DISTRIB]'
      write(out,'(a)') '# Deposited energy over voxel mass'
      if (isFullvol) then
        write(out,'(a)') '#   (INCLUDING parts overlapped by quadrics)'
      else
        write(out,'(a)') '#   (EXCLUDING parts overlapped by quadrics)'
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Dose units are eV/g per primary history'
      write(out,'(a)') '# Number of voxels in x,y,z:'
      write(out,'(a,3(1x,i0))') '# ',nx,ny,nz
      write(out,'(a)') '# Voxels size dx,dy,dz (cm):'
      write(out,'(a,3(1x,es12.5))')'# ',dx,dy,dz
      if (isQuad) then
        write(out,'(a)')
     &    '# The granularity used to compute the voxels mass was:'
        write(out,'(a,i0)') '# ',granul
      endif
      if (prtxyz.eq.1) then
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# For plotting purposes, two values per voxel '//
     &    'coordinate are given, namely,'
        write(out,'(a)')
     &    '#   the low end and the middle point of each voxel'
        write(out,'(a)') '#'
      endif
      write(out,'(a,$)') '# '
      if (prtxyz.eq.1) then
        write(out,'(a,$)') 'xLow(cm) : xMiddle(cm) : '
        write(out,'(a,$)') 'yLow(cm) : yMiddle(cm) : '
        write(out,'(a,$)') 'zLow(cm) : zMiddle(cm) : '
      endif
      write(out,'(a,$)') 'dose (eV/g) : +-2sigma'
      if (prtdens.eq.1) write(out,'(a,$)')
     &  ' : voxel mass (g) : Pure(+)/Overlapped(-)'
      typevox(1) = '-'
      typevox(3) = '+'
      write(out,'(a)') ''  ! End of line

      ! Write data:
      do k=zvoxmin,zvoxmax
        z = dz*(k-1)
        zmiddle = z+dz*0.5d0
        write(out,'(a,i0,$)') '# zVoxIndex=',k
        if (prtxyz.ne.1)
     &    write(out,'(a,es12.5,$)') ' zMiddle(cm)=',zmiddle
        write(out,*) ''  ! EOL

        do j=yvoxmin,yvoxmax
          y = dy*(j-1)
          ymiddle = y+dy*0.5d0
          write(out,'(a,i0,$)') '# yVoxIndex=',j
          if (prtxyz.ne.1)
     &      write(out,'(a,es12.5,$)') ' yMiddle(cm)=',ymiddle
          write(out,*) ''  ! EOL

          do i=xvoxmin,xvoxmax
            vox = i+(j-1)*nx+(k-1)*nxy    ! Absolute voxel index
            if (prtxyz.eq.1) then
              x = dx*(i-1)
              xmiddle = x+dx*0.5d0
              write(out,'(6(1x,es12.5),$)')
     &          x,xmiddle,y,ymiddle,z,zmiddle
            endif
            fact = 0.0d0 ! Voxel mass (may be null if partial vol.):
            if (massvox(vox).gt.0.0d0) fact = 1.0d0/massvox(vox)
            q = edep(vox)*invn
            sigma = sqrt(max((edep2(vox)*invn-q**2)*invn,0.0d0))*fact
            q = q*fact
            write(out,'(1x,es12.5,1x,es7.1,$)') q,2.0d0*sigma
            if (prtdens.eq.1)             ! Print voxel partial mass
     &        write(out,'(1x,es12.5,1x,a1,$)') massvox(vox),
     &          typevox(2+sign(1,matvox(vox))) ! (1) or (3)
            write(out,*) ''               ! End of line
            ! Evaluate average uncertainty for scores above 1/2 max score:
            if (edep(vox).gt.maxq.and.fact.gt.0.0d0) then
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
          enddo

          if (xvoxmax.gt.xvoxmin) write(out,*) '' ! Separate 2D data blocks
        enddo

        if (yvoxmax.gt.yvoxmin) write(out,*) ''   ! Separate 3D data blocks
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
      write(out,'(a,es12.5)') '#   ',cputim
      if (cputim.gt.0.0d0) then
        write(out,'(a)') '#   Speed (histories/s):'
        write(out,'(a,es12.5)') '#   ',n/cputim
      endif
      write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
      write(out,'(a,es12.5)') '#   ',uncert
      eff = n*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
        write(out,'(a,es12.5)') '#   ',1.0d0/eff
      endif
      eff = cputim*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
        write(out,'(a,es12.5)') '#   ',1.0d0/eff
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      close(out)
      end


      subroutine VDDinitally
!*******************************************************************
!*    Initializes the Voxel Dose tally. To be called before TALLY. *
!*    The voxelized geometry must be active to use this tally.     *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 iz,nelem,maxmat
      parameter (maxmat=10)
      real*8 stf,zt,at,rho,vmol
      common/compos/stf(maxmat,30),zt(maxmat),at(maxmat),rho(maxmat),
     &               vmol(maxmat),iz(maxmat,30),nelem(maxmat)
      character*80 voxfilen
      logical isFullvol
      integer matvox,granul,maxGranul
      parameter (maxGranul=1000)
      integer*4 nx,ny,nz,nxy,maxvox,bodymask,matmask
      parameter (maxvox=10000000)  ! Max no. of voxels
      real densvox,idensvox
      real*8 massvox,dx,dy,dz,idx,idy,idz,vbb
      common /geovox/ massvox(maxvox),densvox(maxvox),idensvox(maxvox),
     &  matvox(maxvox),dx,dy,dz,idx,idy,idz,vbb(3),nx,ny,nz,nxy,
     &  bodymask,matmask,granul,isFullvol,voxfilen
      logical isQuad
      common /geoquad/ isQuad
      logical active
      integer prtxyz,prtdens
      integer*4 xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax
      real*8 nlast,nhist,unclimit,edptmp,edep,edep2
      common /scovdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                nlast(maxvox),unclimit,nhist,
     &                xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax,
     &                prtxyz,prtdens,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY VOXEL DOSE v.2008-06-01]')
      parameter (eos='[END OF VDD SECTION]')
      character*80 buffer
      integer error,answer
      integer*4 vox,dk,djk,i,j,k
      real*8 volvox

      write(*,*)
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'VDDinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(1x,a3)') buffer
      if (adjustl(buffer(1:3)).eq.'ON') then
        if (bodymask.ne.-1) then  ! Voxel geometry is ON
          active = .true.
        else
          active = .false.
          write(*,'(a)')
     &     'Voxel Dose tally switch was ON but no voxels are defined'
          write(*,'(a)')
     &     'so the switch has been turned OFF automatically.'
        endif
      else if (buffer(1:3).eq.'OFF') then
        active = .false.
      else
        write(*,'(a)')
     &    'VDDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      if (active) then
        write(*,'(a)') 'Region Of Interest set to:'
        read(*,*) xvoxmin,xvoxmax
        xvoxmin = max(1,xvoxmin)
        xvoxmax = min(nx,xvoxmax)
        if (xvoxmax.eq.0) xvoxmax = nx
        write(*,'(a,2(2x,i0))') '  x-voxel interval: ',xvoxmin,xvoxmax
        read(*,*) yvoxmin,yvoxmax
        yvoxmin = max(1,yvoxmin)
        yvoxmax = min(ny,yvoxmax)
        if (yvoxmax.eq.0) yvoxmax = ny
        write(*,'(a,2(2x,i0))') '  y-voxel interval: ',yvoxmin,yvoxmax
        read(*,*) zvoxmin, zvoxmax
        zvoxmin = max(1,zvoxmin)
        zvoxmax = min(nz,zvoxmax)
        if (zvoxmax.eq.0) zvoxmax = nz
        write(*,'(a,2(2x,i0))') '  z-voxel interval: ',zvoxmin,zvoxmax

        write(*,'(a)') 'Partial volume policy:'
        read(*,*) answer
        select case(answer)
        case(0)
          isFullvol = .false.
          write(*,'(a)')
     &      '  Voxels mass and doses EXCLUDE overlapping quadrics.'
        case(1)
          isFullvol = .true.
          write(*,'(a)')
     &      '  Voxels mass and doses INCLUDE overlapping quadrics.'
        case default
          write(*,'(a)') 'VDDinitally:ERROR: expecting 0 or 1.'
          stop
        end select

        write(*,'(a)') 'Print voxels mass:'
        read(*,*) prtdens
        if (prtdens.eq.1) then
          write(*,'(a)') ' yes'
        else
          write(*,'(a)') ' no'
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
        write(*,'(es12.5)') unclimit

        ! Check section integrity:
        read(*,'(a80)') buffer
        if (index(buffer,eos).eq.0) then
          write(*,'(a)') 'inivox ERROR: End-Of-Section mark not found'
          write(*,'(a,a)') '  expecting to find: ',eos
          write(*,'(a,a)') '  found instead:     ',buffer
          stop
        endif

      else                         ! Tally is inactive
        isFullvol = .false.        ! (in case voxels mass are computed)
      endif

      ! Init counters, compute voxels mass and sign:
      !  (if tally is inactive, only voxels sign matter)
      if (bodymask.ne.-1) then     ! Voxel geometry is ON
        write(*,'(a)') 'Computing voxels mass...'
        volvox = dx*dy*dz          ! Voxels volume
        do k=1,nz                  ! For all voxels
          dk = (k-1)*nxy
          do j=1,ny
            djk = (j-1)*nx+dk
            do i=1,nx
              vox = i+djk          ! Absolute voxel index
              edptmp(vox) = 0.0d0  ! Init deposited energy arrays
              edep(vox)   = 0.0d0
              edep2(vox)  = 0.0d0
              nlast(vox)  = 0.0d0
              if (isQuad) then     ! Quadrics present, init mass to nil
                massvox(vox) = 0.0d0
              else                 ! No quadrics, compute mass as vol*dens
                massvox(vox) = volvox*densvox(vox)
              endif
            enddo
          enddo
        enddo
        if (isQuad) then           ! Quadrics present
          call setMassvox          ! Compute mass and set (-) matvox's
          !! call writeMassvox     ! Write voxels mass and sign to a file
        endif                      !   (reserved for future versions)


        ! Redefine voxels density, normalizing by nominal mat density:
        do k=1,nz                  ! For all voxels
          dk = (k-1)*nxy
          do j=1,ny
            djk = (j-1)*nx+dk
            do i=1,nx
              vox = i+djk
              densvox(vox) = densvox(vox)/rho(abs(matvox(vox)))
              idensvox(vox) = 1.0/densvox(vox)
            enddo
          enddo
        enddo

        write(*,'(a)') 'Done.'
      endif

      if (active) then
        write(*,'(a)') '>>>> VDD tally initialization finished >>>>'
      else
        write(*, '(a)')
     &    '>>>> Tally Voxel Dose Distribution is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'VDDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      endif
      end


      subroutine setMassvox
!*******************************************************************
!*    Computes the mass of all voxels by integrating density over  *
!*    volume, taking into account whether the mode is fullvol or   *
!*    not. It also sets the sign of the voxel material, depending  *
!*    on whether the voxel is pure (+) or overlapped (-).          *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 iz,nelem,maxmat
      parameter (maxmat=10)
      real*8 stf,zt,at,rho,vmol
      common/compos/stf(maxmat,30),zt(maxmat),at(maxmat),rho(maxmat),
     &               vmol(maxmat),iz(maxmat,30),nelem(maxmat)
      character*80 voxfilen
      logical isFullvol
      integer matvox,granul,maxGranul
      parameter (maxGranul=1000)
      integer*4 nx,ny,nz,nxy,maxvox,bodymask,matmask
      parameter (maxvox=10000000)  ! Max no. of voxels
      real densvox,idensvox
      real*8 massvox,dx,dy,dz,idx,idy,idz,vbb
      common /geovox/ massvox(maxvox),densvox(maxvox),idensvox(maxvox),
     &  matvox(maxvox),dx,dy,dz,idx,idy,idz,vbb(3),nx,ny,nz,nxy,
     &  bodymask,matmask,granul,isFullvol,voxfilen
      integer*4 xvox,yvox,zvox,i,j,vox
      integer*4 voxIni,zvoxIni,zvoxFin,ncross,matold
      real*8 da(maxGranul,maxGranul),xg(maxGranul),yg(maxGranul)
      real*8 darea,zinf,xIni,yIni,zIni,dist,dsef,shift,delta,eps,oneeps
      parameter (delta=0.01d0)   ! A small fraction of unity
      parameter (eps=1.0d-10,oneeps=1.0d0+1.0d-12)  ! From PENGEOM

      darea = dx*dy/(granul-1)**2  ! Area of cells in which voxels area subdivided
      ! Set grid points and weights to sweep a voxel according to GRANUL:
      do j=1,granul
        do i=1,granul
          da(i,j) = darea
        enddo
      enddo
      do i=1,granul
        xg(i) = (i-1)*dx/(granul-1)
        da(i,1) = da(i,1)*0.5d0
        da(i,granul) = da(i,granul)*0.5d0  ! Halve weights
      enddo
      do j=1,granul
        yg(j) = (j-1)*dy/(granul-1)
        da(1,j) = da(1,j)*0.5d0
        da(granul,j) = da(granul,j)*0.5d0  ! Halve weights
      enddo
      ! Shift sides:
      shift = delta*dx/(granul-1)
      xg(1) = shift
      xg(granul) = dx-shift
      shift = delta*dy/(granul-1)
      yg(1) = shift
      yg(granul) = dy-shift
      ! ...Note that corners are shifted both in x and y and their areas halved twice

      zinf = 1.01d0*nz*dz ! Practical infinity (>VBB in z) avoids overflows
      u = 0.0d0
      v = 0.0d0
      w = 1.0d0                                ! Pointing upwards
      do yvox=1,ny                             ! For each voxel
        yIni = (yvox-1)*dy                     ! Coordinates of the voxel corner
        do xvox=1,nx
          xIni = (xvox-1)*dx
          voxIni = xvox+(yvox-1)*nx            ! Absolute voxel index for this column

          do j=1,granul                        ! For each ray (i.e. each cell) in a voxel
            do i=1,granul
              darea = da(i,j)                  ! dArea to compute voxel mass
              x = xIni+xg(i)                   ! Coordinates of the {x,y} ray
              y = yIni+yg(j)
              z = 0.0d0                        ! Initial z on the VBB
              call locate                      ! Set ibody and mat
              zvoxIni = 1                      ! Initial z index for this ray
              vox = voxIni                     ! Init voxel index for this ray
              zIni = z                         ! First raystep position

              raystep: do                      ! For each step in a ray
                matold = mat                   ! Save current mat to compute vox mass
                z = zIni                       ! Set start of next step for STEP to work
                call step(zinf,dsef,ncross)    ! Compute next step of current ray
                if (mat.eq.0) then             ! Gone, last step
                  if (matold.eq.0) then        ! From vac to vac, no matter
                    z = +zinf                  ! Set z to practical infinity
                  else
                    z = zIni+oneeps*dsef+eps   ! Reposition in close vacuum
                  endif
                else if (ncross.gt.1) then     ! Matter after vaccum gap
                  z = zIni+oneeps*dsef+eps     ! Reposition inside vac gap
                endif
                zvoxFin = int(z*idz+1.0d0)     ! Final voxel of the current step

                do zvox=max(1,zvoxIni),min(nz,zvoxFin)           ! For each visited voxel along z
                  ! Compute distance travelled inside this voxel:
                  dist = dz                                      ! For a whole voxel
                  if (zvox.eq.zvoxIni) dist = dist-zIni+(zvox-1)*dz  ! Correct start-of-step
                  if (zvox.eq.zvoxFin) dist = dist+z-zvox*dz     ! Correct end-of-step
                  ! Compute mass and mark overlapped voxels:
                  if (matold.eq.matmask) then                    ! Moved inside voxels
                    massvox(vox) = massvox(vox)+darea*dist*densvox(vox)
                  else                                           ! Moved in quadric geometry
                    matvox(vox) = -abs(matvox(vox))              ! Mark voxels as overlapped
                    if (matold.gt.0.and.isFullvol)               ! Protect vacuum and partial vol.
     &               massvox(vox) = massvox(vox)+darea*dist*rho(matold)
                  endif
                  vox = vox+nxy                                  ! Move one voxel up along z
                enddo
                vox = vox-nxy                                    ! Move 'pointer' back to last voxel

                zIni = z                         ! Update next step initial position
                zvoxIni = zvoxFin                ! Update next step initial index
                if (zvoxIni.gt.nz) exit raystep  ! All voxels along this ray have been completed
              enddo raystep                      ! Cycle next step
            enddo
          enddo
        enddo
      enddo
      end


      subroutine writeMassvox
!*******************************************************************
!*    Writes voxels mass and sign to an external file for reference*
!*                                                                 *
!*    Comments:                                                    *
!*    -> Currently unused, awaiting for future versions that would *
!*       allow to load the mass file instead of recalculating them.*
!*******************************************************************
      implicit none
      character*80 voxfilen
      logical isFullvol
      integer matvox,granul,maxGranul
      parameter (maxGranul=1000)
      integer*4 nx,ny,nz,nxy,maxvox,bodymask,matmask
      parameter (maxvox=10000000)  ! Max no. of voxels
      real densvox,idensvox
      real*8 massvox,dx,dy,dz,idx,idy,idz,vbb
      common /geovox/ massvox(maxvox),densvox(maxvox),idensvox(maxvox),
     &  matvox(maxvox),dx,dy,dz,idx,idy,idz,vbb(3),nx,ny,nz,nxy,
     &  bodymask,matmask,granul,isFullvol,voxfilen
      character*90 buffer
      integer ufile,finduf,error
      integer*4 vox,i,j,k,dk,djk

      buffer = trim(voxfilen)//'.mass'
      ufile = finduf()  ! Find a valid unit file
      open(ufile,file=buffer,iostat=error)
      if (error.ne.0) then
        write(*,'(a,a)') 'writeMassvox:ERROR: unable to open file: ',
     &    buffer
        stop
      endif

      write(ufile,'(a,a)')
     &  '# This file was automatically created by penEasy from '//
     &  'the voxels file named: ',trim(voxfilen)
      write(ufile,'(a)')
     &  '# It contains the list of voxels mass (in g) with sign,'
      write(ufile,'(a)')
     &  '#   (+) for pure voxels and (-) for voxels overlapped '//
     &  'by a quadric body.'
      if (isFullvol) then
        write(ufile,'(a)')
     &    '# The voxels mass INCLUDES the volume covered by quadrics.'
      else
        write(ufile,'(a)')
     &    '# The voxels mass EXCLUDES the volume covered by quadrics.'
      endif
      write(ufile,'(a)')
     &  '# This file is for your information; '//
     &  'it can be deleted at any time.'
      write(ufile,'(a,i0)')
     &  '# The granularity used to compute the mass was: ',granul
      write(ufile,'(a)') '#'

      ! Write data to file:
      do k=1,nz
        write(ufile,'(a,i0)') '# zVoxIndex=',k
        dk = (k-1)*nxy
        do j=1,ny
          write(ufile,'(a,i0)') '# yVoxIndex=',j
          djk = (j-1)*nx+dk
          do i=1,nx
            vox = i+djk
            write(ufile,'(es12.5)') sign(massvox(vox),dble(matvox(vox)))
          enddo
          if (nx.gt.1) write(ufile,*) '' ! Separate 2D data blocks
        enddo
        if (ny.gt.1) write(ufile,*) ''   ! Separate 3D data blocks
      enddo
      close(ufile)
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
