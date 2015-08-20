!*******************************************************************
!*                        PENVOX                                   *
!*                                                                 *
!* Short description:                                              *
!*   Performs the rectilinear transport of particles in a geometry *
!*   that combines quadric surfaces and voxels. Routines from      *
!*   PENGEOM are used to handle the quadric part of the geometry.  *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/, /QTREE/                                    *
!*   -> routines STEP,LOCATE                                       *
!*   from other penEasy files:                                     *
!*   -> common /GEOQUAD/                                           *
!*   -> routine FINDUF                                             *
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


      subroutine stepx(ds,dsef,ncross)
!*******************************************************************
!*    Substitute for PENGEOM's STEP that transports a particle     *
!*    across a geometry that combines quadric surfaces and voxels. *
!*    The particle will travel across multiple voxels in a single  *
!*    step if the material does not change. The flight stops if:   *
!*    (i) mat changes; (ii) the particle enters or exits the voxels*
!*    bounding box (VBB); (iii) similarly to STEP, when a vacuum   *
!*    gap is found, and the particle stops on the next non-vacuum  *
!*    mat; (iv) the particle escapes to infinity; and (v) the      *
!*    distance ds has been traversed.                              *
!*                                                                 *
!*    Input:                                                       *
!*      ds -> distance to travel (cm)                              *
!*    Output:                                                      *
!*      dsef -> effective distance travelled (cm)                  *
!*      ncross -> 0 if ds has been completed, >0 else              *
!*    Comments:                                                    *
!*      -> Values of ncross different from 0 may NOT have the same *
!*         meaning as in the original STEP.                        *
!*******************************************************************
      implicit none
      integer*4 ncross
      real*8 ds,dsef

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
      logical ismix
      integer*4 mat0,ncquad
      real*8 x0,y0,z0,dsvox,dsquad,dsvbb,fuel
      real*8 inf,eps,oneeps,distVBB
      parameter (inf=1.0d30)
      parameter (eps=1.0d-10,oneeps=1.0d0+1.0d-12)  ! From PENGEOM

      if (absvox.eq.0) then             ! Particle in quadric geometry
        if (ibody.ne.bodymask) then     !  and not in transparent body
          call step(ds,dsef,ncross)     !  so interrogate quadrics only
          if (ibody.eq.bodymask) call locatevox ! Check voxels
          return                        ! Stop flight
        endif

        ! Else, transparent mat but out of VBB:
        dsvbb = distVBB()*oneeps+eps    ! Distance (enlarged) to VBB
        x0 = x                          ! Save position before STEP
        y0 = y
        z0 = z
        call step(ds,dsef,ncross)       ! Distance to quadrics
        if (dsef.le.dsvbb) then         ! Quadrics are closer
          if (ncross.ne.0.and.ibody.eq.bodymask) call locatevox ! Vacuum gap
        else                            ! Entered voxels
          ibody = bodymask              ! So, restore transparent mat
          ncross = 1                    ! Return value
          dsef = dsvbb                  ! Ditto
          x = x0+dsvbb*u                ! Move particle up to VBB
          y = y0+dsvbb*v
          z = z0+dsvbb*w
          call locatevox                ! Set indices and mat
          if (absvox.eq.0) then         ! Consistency check, VBB missed
            call locate                 ! Reset ibody and mat
            !! write(*,'(a,3(1x,es12.5))')
            !!  'stepx:ERROR: inconsistency found at x,y,z:',x,y,z
            !! write(*,'(a)') '  Expected inside VBB, but it is not.'
            !! stop
          endif
        endif
        return                          ! Stop flight
      endif
      ! ...else particle already in voxels

      ! Get ready for the flight:
      fuel = ds                         ! 'Fuel' is distance in nominal mat density
      dsef = 0.0d0                      ! Reset travelled distance
      mat0 = mat                        ! Save original vox material
      ismix = matvox(absvox).lt.0       ! True if in an overlapped (-) vox

      ! Move particle across voxels of same mat:
      do                                ! Up to next (-) voxel or quadric
        call stepvox(fuel,dsvox,ncross) ! Find distance to (-) voxel walls
        if (ismix) then                 ! It's a (-) voxel
          x0 = x                        ! Save original position
          y0 = y                        !   (STEPVOX does not change it,
          z0 = z                        !    but STEP does)
          call step(+inf,dsquad,ncquad) ! Find distance to quads
          if (dsquad.lt.dsvox+eps) then ! Entering quad geom
            absvox = 0                  ! Mark as quadrics
            dsef = dsef+dsquad          ! Accumulate travelled distance
            ncross = ncquad             ! Set ncross as found by STEP
            if (ibody.eq.bodymask) call locatevox ! Re-entered after vac gap
            return                      ! xyz,ibody,mat set by STEP
          endif
          x = x0                        ! Restore original position ...
          y = y0
          z = z0
          ibody = bodymask              ! ... and body ...
          mat = mat0                    ! ... and material
        endif

        ! Voxels only, quadrics play no role:
        dsef = dsef+dsvox               ! Accumulate travelled distance
        if (ncross.eq.-9999) dsvox = dsvox*oneeps+eps ! Particle left the VBB, push out
        x = x+dsvox*u                   ! Move--no need to enlarge in voxs
        y = y+dsvox*v
        z = z+dsvox*w
        if (ncross.eq.0) return         ! Interaction inside current mat
        if (ncross.gt.0) then           ! Found voxel with new mat
          mat = ncross                  ! Update mat
          return
        endif ! ncross.lt.0             ! Overlap (-) voxel or out of VBB
        ncross = -ncross                ! Prevent returning neg values
        call locate                     ! Reset ibody and mat
        if (ibody.ne.bodymask) then     ! Entered quadrics
          absvox = 0                    ! Mark as quadrics
          if (mat.ne.0) return          ! Not in vacuum, so we're done
          call step(+inf,dsquad,ncquad) ! Skip vacuum
          if (ibody.eq.bodymask) call locatevox ! Check voxels
          return                        ! Truncate flight after a vacuum
        endif
        if (ncross.eq.9999) return      ! In bodymask, but out of VBB
        mat = ncross                    ! Entered a (-) voxel, restore mat
        if (mat.ne.mat0) return         ! Material changed
        ismix = .true.                  ! From 2nd cycle on, overlapped vox
      enddo
      end


      subroutine stepvox(fuel,dsef,ncross)
!*******************************************************************
!*    Starting from a point inside a voxel, this routine moves     *
!*    the particle across *pure* (+) voxels made of the same MAT.  *
!*    The flight stops when: (i) voxel material has changed; (ii)  *
!*    an interaction happens; (iii) entered an overlapped (-) voxel*
!*    or (iv) the particle has exited the voxels bounding box (VBB)*
!*                                                                 *
!*    Input:                                                       *
!*      fuel -> max flight distance if in nominal mat density      *
!*      Voxel indices {xyz}vox,absvox                              *
!*    Output:                                                      *
!*      ds -> remaining fuel (if new mat found, else undef)        *
!*      dsef -> actual distance travelled (cm)                     *
!*      ncross -> 0 if interaction, new material (+ or -), or -9999*
!*      Voxel indices {xyz}vox,absvox                              *
!*    Comments:                                                    *
!*      -> Only voxel indices are updated; neither xyz nor mat are.*
!*      -> If particle leaves the VBB, absvox is set to 0 and      *
!*         ncross is set to -9999.                                 *
!*******************************************************************
      implicit none
      integer*4 ncross
      real*8 fuel,dsef

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
      integer dxvox,dyvox,dzvox,mat0
      real*8 sx,sy,sz,dfuel,inf
      parameter (inf=1.0d30)

      ! Check if direction has changed since last call and store
      !   info to save time when going across multiple voxels:
      if (u.ne.uold.or.v.ne.vold) then
        if (u.ne.0.0d0) then
          ivx = 1.0d0/u                 ! Inverse of direction cosines
          sdx = abs(dx*ivx)             ! Distance to travel to span dx
        else
          ivx = inf
          sdx = inf
        endif
        if (v.ne.0.0d0) then
          ivy = 1.0d0/v
          sdy = abs(dy*ivy)
        else
          ivy = inf
          sdy = inf
        endif
        if (w.ne.0.0d0) then
          ivz = 1.0d0/w
          sdz = abs(dz*ivz)
        else
          ivz = inf
          sdz = inf
        endif
        uold = u                        ! Save new particle's direction
        vold = v
      endif

      ! Determine min distance to the three relevant voxel walls:
      if (ivx.gt.0.0d0) then
        sx = (xvox*dx-x)*ivx
        dxvox = +1                          ! Will move one voxel forward
      else
        sx = ((xvox-1)*dx-x)*ivx
        dxvox = -1                          ! Will move one voxel backwards
      endif
      if (ivy.gt.0.0d0) then
        sy = (yvox*dy-y)*ivy
        dyvox = +1
      else
        sy = ((yvox-1)*dy-y)*ivy
        dyvox = -1
      endif
      if (ivz.gt.0.0d0) then
        sz = (zvox*dz-z)*ivz
        dzvox = +1
      else
        sz = ((zvox-1)*dz-z)*ivz
        dzvox = -1
      endif

      dsef = 0.0d0                          ! Reset traveled distance
      mat0 = mat                            ! Store original material
      do                                    ! One voxel per cycle
        if (sx.lt.sy) then
          if (sx.lt.sz) then                ! Intersection with x wall

            sx = max(sx,0.0d0)              ! Never go backwards
            dfuel = sx*densvox(absvox)      ! Dist2Surf*dens(vox)/dens(nominal)
            if (dfuel.gt.fuel) then         ! Run out of fuel => interaction
              ncross = 0                    ! Signals an interaction
              dsef = dsef+fuel*idensvox(absvox) ! Distance up to interaction
              return                        ! Stop flight
            endif
            dsef = dsef+sx                  ! Update travelled distance
            xvox = xvox+dxvox               ! Update coordinate voxel index
            if (xvox.lt.1.or.xvox.gt.nx) then  ! Particle is gone
              ncross = -9999
              absvox = 0                    ! Mark as quadrics
              return
            endif
            fuel = fuel-dfuel               ! Remove spent fuel from tank
            absvox = absvox+dxvox           ! Update abs voxel index
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return      ! MAT changed or (-) voxel
            sy = sy-sx                      ! Update distances to walls
            sz = sz-sx
            sx = sdx                        ! Reset to full voxel width
            cycle                           ! Process next voxel

          else                              ! Intersection with z wall
                                            ! Repeat code above for z
            sz = max(sz,0.0d0)
            dfuel = sz*densvox(absvox)
            if (dfuel.gt.fuel) then
              ncross = 0
              dsef = dsef+fuel*idensvox(absvox)
              return
            endif
            dsef = dsef+sz
            zvox = zvox+dzvox
            if (zvox.lt.1.or.zvox.gt.nz) then
              ncross = -9999
              absvox = 0
              return
            endif
            fuel = fuel-dfuel
            absvox = absvox+dzvox*nxy
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return
            sx = sx-sz
            sy = sy-sz
            sz = sdz
            cycle

          endif
        else
          if (sy.lt.sz) then                ! Intersection with y wall
                                            ! Repeat code above for y
            sy = max(sy,0.0d0)
            dfuel = sy*densvox(absvox)
            if (dfuel.gt.fuel) then
              ncross = 0
              dsef = dsef+fuel*idensvox(absvox)
              return
            endif
            dsef = dsef+sy
            yvox = yvox+dyvox
            if (yvox.lt.1.or.yvox.gt.ny) then
              ncross = -9999
              absvox = 0
              return
            endif
            fuel = fuel-dfuel
            absvox = absvox+dyvox*nx
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return
            sz = sz-sy
            sx = sx-sy
            sy = sdy
            cycle

          else                              ! Intersection with z wall
                                            ! Repeat code above for z
            sz = max(sz,0.0d0)
            dfuel = sz*densvox(absvox)
            if (dfuel.gt.fuel) then
              ncross = 0
              dsef = dsef+fuel*idensvox(absvox)
              return
            endif
            dsef = dsef+sz
            zvox = zvox+dzvox
            if (zvox.lt.1.or.zvox.gt.nz) then
              ncross = -9999
              absvox = 0
              return
            endif
            fuel = fuel-dfuel
            absvox = absvox+dzvox*nxy
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return
            sx = sx-sz
            sy = sy-sz
            sz = sdz
            cycle

          endif
        endif
      enddo
      end


      real*8 function distVBB()
!*******************************************************************
!*    Returns distance to intersection between the particle's      *
!*    trajectory and the VBB. If there is no intersection, +inf is *
!*    returned.                                                    *
!*                                                                 *
!*    Input:                                                       *
!*      -> xyz, uvw, vbb needs to be defined.                      *
!*    Comments:                                                    *
!*      -> Max considered distance is inf=1e90 cm.                 *
!*******************************************************************
      implicit none
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
      integer i
      real*8 r(3),d(3),inf,slo,shi,invd
      parameter (inf=1.0d90)

      r = (/ x,y,z /)         ! Starting point
      d = (/ u,v,w /)         ! Direction
      slo = 0.0d0      ! Position 'low pointer' at start of segment
      shi = +inf       ! Position 'high pointer' at end of segment
      distVBB = +inf   ! Infinity signals 'no intersection'

      do i=1,3                         ! For x,y,z
        if (d(i).eq.0.0d0) then        ! Going parallel to surface
          if (r(i).lt.0.0d0.or.r(i).gt.vbb(i)) return ! No VBB inters.
          cycle                        ! Otherwise, no need to clip
        endif
        invd = -1.0d0/d(i)             ! Avoids repeating the division
        if (invd.lt.0.0d0) then        ! Going forward
          if (r(i).gt.vbb(i)) return ! No intersection
          shi = min(shi,(r(i)-vbb(i))*invd)         ! Clip high end
          if (r(i).lt.0.0d0) slo = max(slo,r(i)*invd) ! Clip low end
        else                           ! Going backwards
          if (r(i).lt.0.0d0) return    ! No intersection
          shi = min(shi,r(i)*invd)     ! Clip high end
          if (r(i).gt.vbb(i)) slo = max(slo,(r(i)-vbb(i))*invd) ! Low end
        endif
        if (shi.lt.slo) return         ! No intersection
      enddo
      distVBB = slo                    ! Distance to VBB
      end


      subroutine locatex
!*******************************************************************
!*    Substitute for PENGEOM's LOCATE that locates a particle      *
!*    in the geometry composed of quadric objects and voxels.      *
!*******************************************************************
      implicit none
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

      call locate                         ! Locate in quadrics
      if (ibody.eq.bodymask) then
        call locatevox                    ! Check voxels
      else
        absvox = 0                        ! Mark as quadrics
      endif
      end


      subroutine locatevox
!*******************************************************************
!*    Checks if inside voxels bounding box (VBB). Determines the   *
!*    voxel that contains the point with coordinates {x,y,z}       *
!*    and sets the voxel indices and the material. If the position *
!*    is outside the VBB, the absolute index is set to 0.          *
!*******************************************************************
      implicit none
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

      xvox = x*idx+1.0d0
      if (xvox.lt.1.or.xvox.gt.nx) then   ! Outside VBB
        absvox = 0                        ! Set mark to quadrics
        return                            ! Done
      endif
      yvox = y*idy+1.0d0
      if (yvox.lt.1.or.yvox.gt.ny) then
        absvox = 0
        return
      endif
      zvox = z*idz+1.0d0
      if (zvox.lt.1.or.zvox.gt.nz) then
        absvox = 0
        return
      endif
      absvox = xvox+(yvox-1)*nx+(zvox-1)*nxy  ! Abs voxel index
      mat = abs(matvox(absvox))
      end


      subroutine inivox(nmatvox)
!*******************************************************************
!*    Initializess the voxel geometry package, setting material and*
!*    density arrays.                                              *
!*                                                                 *
!*    Output:                                                      *
!*      nmatvox -> no. of materials in VOX file.                   *
!*    Comments:                                                    *
!*      -> Computation of voxels mass and sign NOT performed since *
!*         PENELOPE has not been initialized, and so the mass      *
!*         density for each material has yet to be defined. Mass   *
!*         and sign will be determined by the voxels dose routine. *
!*******************************************************************
      implicit none
      integer*4 nmatvox

      integer*4 gnb,gnx,gns ! PENGEOM changed to avoid conflicting names
      parameter (gns=10000,gnb=5000,gnx=250)
      integer*4 nbody,mater,kmoth,kdght,ksurf,kflag,kalias,kslast
      common/qtree/nbody,mater(gnb),kmoth(gnb),kdght(gnb,gnx),
     &    ksurf(gnb,gnx),kflag(gnb,gnx),kalias(gns),kslast
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
      logical isQuad
      common /geoquad/ isQuad
      character*80 buffer
      integer nmask
      integer*4 i

      ! Read voxels file, if there is one:
      read(*,'(1x,a30)') voxfilen ! Voxels filename
      write(*,*) ''
      if (len_trim(voxfilen).eq.0) then ! No voxels
        write(*,'(a)') 'No voxels file defined, using quadrics only.'
        read(*,'(a80)') buffer   ! Dummy (unused) lines from input
        read(*,'(a80)') buffer
        nmatvox = 0              ! No voxel materials
        bodymask = -1            ! Impossible body mask signals no voxels
        return                   ! Nothing else to do
      endif                      ! Else, a vox file exists
      write(*,'(a)')
     &  '>>>> Initializing voxelized geometry >>>>'
      write(*,'(a)')
     &  'NOTICE: when voxels are active, some tallies may not yield'
      write(*,'(a)')
     &  '  correct values if evaluated inside the voxelized region'
      write(*,'(a)')
     &  '  (e.g. fluence, cylindrical and spherical dose distrib).'
      write(*,'(a)')
     &  '  Read supplied documentation for more info.'
      call readvox(nmatvox)

      if (isQuad) then           ! Quadrics geom present
        write(*,'(a)') 'Quadric material that unveils voxels:'
        read(*,*) matmask
        write(*,'(1x,i0)') matmask
        if (matmask.lt.1) then
          write(*,'(a)') 'inivox:ERROR: mask material must be >0.'
          stop  ! Vacuum forbidden as transparent mat (but allowed elsewhere)
        endif   !   because it would be missed by STEP
        write(*,'(a)')
     &    'Searching for bodies containing the transparent material...'
        nmask = 0                ! No. of bodies with matmask
        do i=1,gnb
          if (mater(i).eq.matmask) then  ! This body contains the transparent mat
            bodymask = i
            nmask = nmask+1
            write(*,'(a,i0)')
     &        '  found in body no. (PENGEOM internal coding): ',i
          endif
        enddo
        if (nmask.ne.1) then
          write(*,'(a)')
     &      'inivox:ERROR: there must be one and only one such body.'
          stop
        endif

        write(*,'(a)') 'Voxel scan granularity:'
        read(*,*) granul  ! No. rays along voxel side (granul**2 per voxel)
        write(*,'(1x,i0)') granul
        if (granul.lt.2.or.granul.gt.maxGranul) then
          write(*,'(a,i0)')
     &      'inivox:ERROR: granularity must be between 2 and ',
     &      maxGranul
          stop
        endif

      else                       ! No quadric geometry
        write(*,'(a)') 'No quadrics present, voxels assumed in vacuum.'
        read(*,*) matmask        ! Read a dummy number
        matmask = 0              ! Unused, but just in case
        if (mater(1).eq.0) then  ! Only one body is defined in def GEO
          bodymask = 1
        else
          write(*,'(a,i0)')
     &      'inivox:ERROR: internal inconsistency; expecting '//
     &      'vacuum, found instead: ',mater(1)
          stop
        endif

        read(*,*) granul         ! Read a dummy number
        granul = 2               ! Unused, but just in case
      endif

      ! Initialize particle state variables with dummy values:
      uold = 1.0d0
      vold = 1.0d0
      ivx = 1.0d0
      ivy = 1.0d0
      ivz = 1.0d0

      ! Note that voxels mass and sign are not computed now;
      !   see header for comments on this.

      write(*,*) ''
      write(*,'(a)')
     &  '>>>> Voxelized geometry initialization finished >>>>'
      end


      subroutine readvox(nmat)
!*******************************************************************
!*    Reads the voxel geometry file.                               *
!*                                                                 *
!*    Output:                                                      *
!*      nmat -> no. of materials in VOX file.                      *
!*******************************************************************
      implicit none
      integer*4 nmat

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
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION VOXELS HEADER v.2008-04-13]')
      parameter (eos='[END OF VXH SECTION]')
      character*500 buffer
      integer error,finduf,ufile,mat,isblank
      integer colMat,colDens,maxCol,ncol,n
      parameter (maxCol=20)
      integer*4 i,j,k,line,vox
      real zero,dens,col(maxCol)
      real*8 rnvox,minvoxSide,eps
      parameter (eps=1.0d-10,minvoxSide=eps*1.0d4,zero=1.0e-30)

      nmat = 0
      write(*,*) ''
      write(*,'(a)') '>>>> Reading voxels file >>>>'
      write(*,'(a)') 'Opening file:'
      write(*,'(1x,a30)') voxfilen
      ufile = finduf()  ! Find a valid unit file
      open(ufile,file=voxfilen,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)') 'readvox:ERROR: unable to open voxels file'
        stop
      endif

      line = 0                           ! Clear file line number
      call getvoxline(ufile,line,buffer) ! Get a non-comment line
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'readvox:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)') 'No. of voxels in x,y,z and total:'
      line = line+1
      read(ufile,*) nx,ny,nz  ! Note that no comments are allowed inside section
      if (min(nx,ny,nz).lt.1) then
        write(*,'(a)') 'readvox:ERROR: no. voxels must be >0.'
        stop
      endif
      rnvox = dble(nx)*dble(ny)*dble(nz)
      write(*,'(4(1x,i0))') nx,ny,nz,int(rnvox)
      if (rnvox.gt.dble(maxvox)) then
        write(*,'(a)')
     &     'readvox:ERROR: Too many voxels, increase MAXVOX '//
     &     'and recompile.'
        stop
      endif
      nxy = nx*ny

      write(*,'(a)') 'Voxel dimensions in x,y,z (cm):'
      line = line+1
      read(ufile,*) dx,dy,dz
      write(*,'(3(1x,es12.5))') dx,dy,dz
      if (min(dx,dy,dz).lt.minvoxSide) then
        write(*,'(a)')
     &    'readvox:ERROR: voxel side too small, tracking algorithm'
        write(*,'(a,es12.5)')
     &    '  requires voxel sides to be larger than (cm):',minvoxSide
        stop
      endif
      write(*,'(a)') 'Voxels volume (cm^3):'
      write(*,'(1x,es12.5)') dx*dy*dz
      idx = 1.0d0/dx
      idy = 1.0d0/dy
      idz = 1.0d0/dz
      vbb = (/dx*nx, dy*ny, dz*nz /)
      write(*,'(a)') 'Size of voxels bounding box, VBB (cm):'
      write(*,'(3(1x,es12.5))') (vbb(i), i=1,3)
      if (maxval(vbb).gt.1.0d5) then
        write(*,'(a)') 'readvox:ERROR: VBB too large.'
        stop
      endif

      line = line+1
      read(ufile,*) colMat
      write(*,'(a)') 'Material ID expected in column:'
      write(*,'(1x,i0)') colMat
      line = line+1
      read(ufile,*) colDens
      write(*,'(a)') 'Mass density expected in column:'
      write(*,'(1x,i0)') colDens
      ncol = max(colMat,colDens)
      if (min(colMat,colDens).lt.1.or.ncol.gt.maxCol)
     &  then
        write(*,'(a,i0)')
     &    'readvox:ERROR: column numbers must be between 1 and ',
     &     maxCol
        stop
      endif

      line = line+1
      read(ufile,*) isblank
      if (isblank.eq.1) then
        write(*,'(a)')
     &    'Expecting blank lines to separate y,z-cycles'
      else
        write(*,'(a)')
     &    'NOT expecting blank lines to separate y,z-cycles'
      endif

      ! Check section integrity:
      line = line+1
      read(ufile,'(a500)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)') 'readvox ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)') 'Now loading mat and dens for each voxel...'
      do k=1,nz
        do j=1,ny
          do i=1,nx

            call getvoxline(ufile,line,buffer)
            read(buffer,*,iostat=error) (col(n), n=1,ncol)
            mat  = int(col(colMat)+0.1)  ! Round to integer
            dens = col(colDens)
            if (error.ne.0.or.mat.lt.1.or.dens.lt.0.0) then
              write(*,'(a,i0)')
     &          'readvox:ERROR: invalid material or density at line ',
     &          line
              write(*,'(a)') '  Line contents:'
              write(*,'(a)') buffer
              stop
            endif
            vox = i+(j-1)*nx+(k-1)*nxy        ! Absolute voxel index
            matvox(vox)  = mat
            densvox(vox) = max(zero,dens)     ! Min dens is +zero
            idensvox(vox)= 1.0/densvox(vox)   ! Inverse vox dens (neglects overlap)
            ! Note: densvox will be redefined latter by VDDinitally() to normalize
            !   to the nominal mass density, not known at this point.
            nmat         = max(nmat,mat)      ! Largest mat index in file
          enddo

          if (isblank.eq.1.and.(j.lt.ny.or.k.lt.nz)) then
            call getvoxline(ufile,line,buffer) ! Read one blank line except in the last line
            if (buffer.ne.' ') then
              write(*,'(a,i0)')
     &          'readvox:ERROR: Line should be blank, line ',line
              stop
            endif
          endif
        enddo

        if (isblank.eq.1.and.k.lt.nz) then     ! Read one blank line except in the last line
          call getvoxline(ufile,line,buffer)
          if (buffer.ne.' ') then
            write(*,'(a,i0)')
     &        'readvox:ERROR: Line should be blank, line ',line
            stop
          endif
        endif
      enddo

      close(ufile)
      write(*,'(a,i0,a)') 'Done. ',line,' lines processed.'
      end


      subroutine getvoxline(ufile,line,buffer)
!*******************************************************************
!*    Reads a new line from an external unit. The line is returned *
!*    only if it is not a comment line.                            *
!*                                                                 *
!*    Input:                                                       *
!*      ufile -> file unit to be read                              *
!*      line -> last line number read from that file               *
!*    Output:                                                      *
!*      line -> updated line number                                *
!*      buffer -> line read from file                              *
!*    Comments:                                                    *
!*      -> comment lines start with a '#'                          *
!*******************************************************************
      implicit none
      character*500 buffer
      integer ufile
      integer*4 line

      integer error

      do
        line = line+1
        read(ufile,'(a500)',iostat=error) buffer
        if (error.ne.0) then
          write(*,'(a)')
     &      'getvoxline:ERROR: unable to read vox file at line:'
          write(*,'(i0)') line
          stop
        endif
        if (buffer(1:1).ne.'#') exit  ! Else, a comment line, read another
      enddo
      end


      subroutine autogeo(filename)
!*******************************************************************
!*    Creates a PENGEOM file with an empty (filled with vacuum)    *
!*    geometry. Used when no GEO file has been provided.           *
!*                                                                 *
!*    Input:                                                       *
!*      filename -> name of the file to be created                 *
!*******************************************************************
      implicit none
      character*(*) filename

      character*80 linezero,lineend,linesurf,lineindex,linescale
      character*80 lineshift,linebody,linemat,linepoint
      parameter (linezero='0000000000000000000000000000000000000'//
     &                    '000000000000000000000000000')
      parameter (lineend ='END      0000000000000000000000000000'//
     &                    '000000000000000000000000000')
      parameter (linesurf  ='SURFACE (')
      parameter (lineindex ='INDICES=(')
      parameter (linescale ='-SCALE=(+')
      parameter (lineshift ='-SHIFT=(+')
      parameter (linebody  ='BODY    (')
      parameter (linemat   ='MATERIAL(')
      parameter (linepoint =', SIDE POINTER=(')
      integer ufile,finduf

      ufile = finduf()
      open(ufile,file=filename)
      write(ufile,'(a)')
     & '# This PENGEOM file was automatically created by penEasy.'
      write(ufile,'(a)')
     & '# It can be removed only after the execution is completed.'
      write(ufile,'(a)')
     & '# Do not edit this file while the code is running.'
      write(ufile,*) ''
      write(ufile,'(a80)') linezero         ! Start with a line of zeroes
      write(ufile,'(a9,i4,")")') linebody,1 ! Define a body
      write(ufile,'(a9,i4,")")') linemat,0  ! Filled with vacuum (no surfs)
      write(ufile,'(a80)') linezero         ! Another line of zeroes
      write(ufile,'(a80)') lineend          ! And the end line
      close(ufile)

      ! Other possible elements:
      !!! Pair of X planes:
      !!write(ufile,'(a9,i4,")")') linesurf,1
      !!write(ufile,'(a9,4(i2,","),i2,")")') lineindex,1,0,0,0,-1
      !!write(ufile,'(a1,a9,1pe21.15,",",i4,")")') 'X',
     &!!  linescale,xscale,0
      !!write(ufile,'(a1,a9,1pe21.15,",",i4,")")') 'X',
     &!!  lineshift,xscale,0
      !!write(ufile,'(a80)') linezero
      !!! Limiting quadrics:
      !!write(ufile,'(a9,i4,")",a16,i2,")")') linesurf,1,linepoint,-1
      !!write(ufile,'(a9,i4,")",a16,i2,")")') linesurf,2,linepoint,-1
      !!write(ufile,'(a9,i4,")",a16,i2,")")') linesurf,3,linepoint,-1
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
