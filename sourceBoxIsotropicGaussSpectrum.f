!*******************************************************************
!*                         SOURCE                                  *
!*               BOX-ISOTROPIC-GAUSS-SPECTRUM                      *
!*                                                                 *
!* Short description:                                              *
!*   Generation of primary particle states for radiation transport *
!*   calculations with PENELOPE.                                   *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routines STORES,RAND                                       *
!*   from PENVOX.F                                                 *
!*   -> common /PARTVOX/                                           *
!*   -> routine LOCATEX                                            *
!*   -> routine STEPX                                              *
!*   from other penEasy libraries:                                 *
!*   -> routines GETLINE,TALLY                                     *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*******************************************************************


      subroutine BIGSsource(n)
!*******************************************************************
!*    Input:                                                       *
!*      n -> History no.                                           *
!*    Output:                                                      *
!*      through /track/ and sec stack                              *
!*******************************************************************
      implicit none
      real*8 n

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 xvox,yvox,zvox,absvox
      real*8 uold,vold,ivx,ivy,ivz,sdx,sdy,sdz
      common /partvox/ uold,vold,ivx,ivy,ivz,sdx,sdy,sdz,
     &                 xvox,yvox,zvox,absvox
      logical warned,srcpoint,active
      integer parsrc,matsrc,nspc,dim
      parameter (dim=1000)
      real*8 shots,usrc,vsrc,wsrc,cossrc,espc,pspc,despc,rot
      real*8 xsrc,ysrc,zsrc,dxsrc,dysrc,dzsrc
      common /srcbig/ rot(3,3),espc(dim),pspc(dim),despc(dim),shots,
     &                cossrc,usrc,vsrc,wsrc,xsrc,ysrc,zsrc,dxsrc,
     &                dysrc,dzsrc,parsrc,matsrc,nspc,warned,srcpoint,
     &                active
      integer ntrial,i,seeki
      integer*4 ncross,bodytmp,mattmp,xvoxtmp,yvoxtmp,zvoxtmp,voxtmp
      real*8 costhe,phi,rand,infty,dsef,xrot,yrot,zrot
      real*8 pi,dospi,one3rd,randno
      parameter (pi=3.1415926535897932d0,dospi=2.0d0*pi,infty=1.0d30)
      parameter (one3rd=1.0d0/3.0d0)
      external rand

      if (.not.active) return

      ntrial = 0
      kpar = parsrc  ! Particle type

      direction: do
        u = usrc
        v = vsrc
        w = wsrc
        costhe = cossrc+(1.0d0-cossrc)*rand(1.0d0)
        phi = dospi*rand(2.0d0)
        call rotate(u,v,w,costhe,phi)

        position: do
          ntrial = ntrial+1
          if (n.lt.1.01d0.and.ntrial.gt.100) then
            write(*,*)
     &    'BIGSsource:ERROR: could not generate a valid particle '//
     &    'after 100 trials;'
            write(*,*)
     &       '  box enclosure does not contain the source '//
     &       'material or direction aims at vacuum'
            stop
          endif
          ! Sample coordinates relative to the box center:
          if (srcpoint) then  ! No need to sample for a point source
            x = xsrc
            y = ysrc
            z = zsrc
          else
            x = dxsrc*(rand(3.0d0)-0.5d0)
            y = dysrc*(rand(4.0d0)-0.5d0)
            z = dzsrc*(rand(5.0d0)-0.5d0)
            ! Rotate:
            xrot = rot(1,1)*x+rot(1,2)*y+rot(1,3)*z
            yrot = rot(2,1)*x+rot(2,2)*y+rot(2,3)*z
            zrot = rot(3,1)*x+rot(3,2)*y+rot(3,3)*z
            x = xrot+xsrc     ! Shift to the box center position
            y = yrot+ysrc
            z = zrot+zsrc
          endif
          call locatex        ! Finds body and material (needs U,V,W)

          ! Accept or reject?
          if (matsrc.eq.0) then                             ! Any material will do, even vacuum
            if (mat.eq.0) call stepx(infty,dsef,ncross)     ! Vacuum, advance up to object or infinity
            exit direction                                  ! Always accept
          else if (srcpoint) then                           ! Particle must aim at specified material
            if (mat.eq.matsrc) exit direction               ! Already inside specified material
            bodytmp = ibody                                 ! Store original position variables
            mattmp = mat
            xvoxtmp = xvox
            yvoxtmp = yvox
            zvoxtmp = zvox
            voxtmp  = absvox
            do
              call stepx(infty,dsef,ncross)                 ! Move forward to find specified material
              if (mat.eq.0.or.ncross.eq.0) cycle direction  ! Not aiming at specified material
              if (mat.eq.matsrc) then                       ! Accept this point
                x = xsrc                                    ! Restore original position
                y = ysrc
                z = zsrc
                ibody  = bodytmp
                mat    = mattmp
                xvox   = xvoxtmp
                yvox   = yvoxtmp
                zvox   = zvoxtmp
                absvox = voxtmp
                if (mat.eq.0) call stepx(infty,dsef,ncross) ! Advance up to object
                exit direction                              ! Accepted
              endif
            enddo
          else if (mat.eq.matsrc) then                      ! Materials match
            exit direction
          endif                                             ! Re-sample XYZ if materials don't match
        enddo position
      enddo direction

      ! Kinetic energy:
      if (nspc.eq.2.and.pspc(2).eq.0.0d0) then          ! Gaussian, sample using Box-Muller
        do
          e = espc(2)+espc(1)*                          ! Were set to Emean,sigma respectively
     &         sqrt(-2.0d0*log(rand(1.8d1)))*sin(dospi*rand(1.9d1))
          if (e.gt.espc(3).and.e.lt.espc(4)) exit       ! Gaussian between Emin and Emax
        enddo
      else                                              ! User-defined spectrum
        randno = rand(2.0d1)
        i = seeki(pspc,randno,nspc)
        e = espc(i)+(randno-pspc(i))*despc(i)
      endif

      ! Average no. of shots per call:
      shots = shots+ntrial
      if (.not.warned.and.shots.gt.n*10.0d0) then
        write(*,*) ' '
        write(*,'(a)')
     &    '***************'
        write(*,'(a)')
     &    'BIGSsource:WARNING: source effectiveness is too low !'
        write(*,'(a)')
     &    '  Redefine source parameters appropriately.'
        write(*,'(a)')
     &    '  Histories and shots per source call so far:'
        write(*,'(2x,f18.0,1x,1pe12.5)') n,shots/n
        write(*,'(a)')
     &    '***************'
        write(*,*) ' '
        warned = .true.
      endif

      wght = 1.0d0                              ! Init statistical weight
      ilb(1) = 1                                ! Tag as source particle (i.e. 1st generation)
      ilb(2) = 0                                ! Clear other labels
      ilb(3) = 0
      ilb(4) = 0
      ilb(5) = 0                                ! Optional label (transferred to descendants)

      call stores(e,x,y,z,u,v,w,wght,kpar,ilb)  ! Push particle to stack
      call tally(0,e)                           ! Deposit its kinetic E
      end


      subroutine BIGSinisrc(activated,emax)
!*******************************************************************
!*    Initializes the source.                                      *
!*                                                                 *
!*    Output:                                                      *
!*      activated -> TRUE if the source is active.                 *
!*      emax -> max source energy (eV)                             *
!*******************************************************************
      implicit none
      logical activated
      real*8 emax

      logical warned,srcpoint,active
      integer parsrc,matsrc,nspc,dim
      parameter (dim=1000)
      real*8 shots,usrc,vsrc,wsrc,cossrc,espc,pspc,despc,rot
      real*8 xsrc,ysrc,zsrc,dxsrc,dysrc,dzsrc
      common /srcbig/ rot(3,3),espc(dim),pspc(dim),despc(dim),shots,
     &                cossrc,usrc,vsrc,wsrc,xsrc,ysrc,zsrc,dxsrc,
     &                dysrc,dzsrc,parsrc,matsrc,nspc,warned,srcpoint,
     &                active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION SOURCE BOX ISOTROPIC GAUSS SPECTRUM v.2006-08-01]')
      parameter (eos='[END OF BIGS SECTION]')
      integer j,error
      real*8 omega,theta,phi,comega,somega,ctheta,stheta,cphi,sphi
      real*8 prob,pi,norm,fwhm2sig,emean,sigma,deg2rad,emin,kfact
      real*8 mc2,twomc2,zero
      parameter (pi=3.1415926535897932d0,deg2rad=pi/180.0d0)
      parameter (fwhm2sig=4.2466090014400952d-1,kfact=6.0d0)
      parameter (mc2=5.10998918d5,twomc2=2.0d0*mc2,zero=1.0d-30)

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'BIGSinisrc:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(1x,a3)') buffer
      if (adjustl(buffer(1:3)).eq.'ON') then
        active = .true.
        activated = active
      else if (buffer(1:3).eq.'OFF') then
        active = .false.
        activated = active
        write(*, '(a)')
     &    '>>>> Source Box Isotropic Gauss Spectrum is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'BIGSinisrc:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'BIGSinisrc:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Type of particle:
      write(*,'(a)') 'Particle type:'
      read(*,*) parsrc
      write(*,'(1x,i1)') parsrc
      if (parsrc.ne.1.and.parsrc.ne.2.and.parsrc.ne.3) then
        write(*,*) 'BIGSinisrc:ERROR: invalid particle type'
        stop
      endif

      ! Energy:
      write(*,'(a)') 'Reading energy spectrum'
      read(*,'(a80)') buffer
      write(*,'(a)') '  Energy(eV)  Relat.Probability   Bin#'
      nspc = 0
      pspc(1) = 0.0d0
      do
        nspc = nspc+1
        read(*,*) espc(nspc),prob
        write(*,'(2(1x,1pe12.5),1x,i5)') espc(nspc),prob,nspc
        if (espc(nspc).lt.0.0d0) then
          write(*,*) 'BIGSinisrc:ERROR: negative energy'
          stop
        else if (espc(nspc).lt.espc(max(nspc-1,1))) then
          write(*,*) 'BIGSinisrc:ERROR: decreasing energy'
          stop
        endif
        if (prob.lt.0.0d0) exit    ! End of spectrum
        if (nspc.ge.dim) then
          write(*,*) 'BIGSinisrc:ERROR: too many bins in spectrum;'
          write(*,*) '              enlarge DIM'
          stop
        endif
        pspc(nspc+1) = pspc(nspc)+prob
      enddo
      write(*,'(a)') 'No. of bins read:'
      write(*,'(1x,i0)') nspc-1
      if (nspc.lt.2) then
        write(*,*)
     &    'BIGSinisrc:ERROR: at least 1 bin must be defined'
        stop
      endif
      if (pspc(nspc).gt.0.0d0) then
        write(*,'(a)')
     &    'Sum of relative probabilities before normalization:'
        write(*,'(1pe12.5)') pspc(nspc)
        do j=1,nspc    ! Normalize to unity
          pspc(j) = pspc(j)/pspc(nspc)
        enddo
        do j=1,nspc-1    ! Prepare auxiliar array for sampling
          despc(j) = 0.0d0
          if (pspc(j+1).gt.pspc(j))
     &      despc(j) = (espc(j+1)-espc(j))/(pspc(j+1)-pspc(j))
        enddo
        emax = espc(nspc)  ! Set max energy
      else
        write(*,'(a)')
     &    'Null probability; assuming Gaussian spectrum'
        if (nspc.gt.2) then
          write(*,*)
     &      'BIGSinisrc:ERROR: Gaussian requires only 1 bin'
          stop
        endif
        emean = 0.5d0*(espc(1)+espc(2))
        sigma = fwhm2sig*(espc(2)-espc(1))
        emin = max(0.0d0,emean-kfact*sigma)
        emax = emean+kfact*sigma
        write(*,'(a)')
     &    'Mean energy, FWHM, sigma, Emin, Emax (eV) of Gaussian:'
        write(*,'(5(1x,1pe12.5))')
     &    emean,espc(2)-espc(1),sigma,emin,emax
        if (emean.lt.0.0d0) then
          write(*,*) 'BIGSinisrc:ERROR: negative mean energy.'
          stop
        endif
        espc(4) = emax
        espc(3) = emin
        espc(2) = emean
        espc(1) = sigma
      endif
      if (parsrc.eq.3) emax = emax+twomc2  ! Allowance for e+ annihilation

      ! Position:
      write(*,'(a)') 'Center coordinates (cm):'
      read(*,*) xsrc,ysrc,zsrc
      write(*,'(3(1x,1pe12.5))') xsrc,ysrc,zsrc
      write(*,'(a)') 'Box sides (cm):'
      read(*,*) dxsrc,dysrc,dzsrc
      write(*,'(3(1x,1pe12.5))') dxsrc,dysrc,dzsrc
      write(*,'(a)') 'Euler angles (deg):'
      read(*,*) omega,theta,phi
      write(*,'(3(1x,1pe12.5))') omega,theta,phi
      ! Rotation matrix:
      somega = sin(omega*deg2rad)
      comega = cos(omega*deg2rad)
      stheta = sin(theta*deg2rad)
      ctheta = cos(theta*deg2rad)
      sphi = sin(phi*deg2rad)
      cphi = cos(phi*deg2rad)
      rot(1,1) = cphi*ctheta*comega-sphi*somega
      rot(1,2) = -cphi*ctheta*somega-sphi*comega
      rot(1,3) = cphi*stheta
      rot(2,1) = sphi*ctheta*comega+cphi*somega
      rot(2,2) = -sphi*ctheta*somega+cphi*comega
      rot(2,3) = sphi*stheta
      rot(3,1) = -stheta*comega
      rot(3,2) = stheta*somega
      rot(3,3) = ctheta
      ! Material:
      if (dxsrc.eq.0.0d0.and.dysrc.eq.0.0d0.and.dzsrc.eq.0.0d0)
     &  srcpoint = .true.  ! Needed to interpret source material
      write(*,'(a)') 'Source material:'
      read(*,*) matsrc
      write(*,'(1x,i0)') matsrc
      if (srcpoint.and.matsrc.ne.0) write(*,'(a)')
     &  '  (interpreted as the material particles must be aiming at)'
      if (matsrc.lt.0) then
        write(*,*) 'BIGSinisrc:ERROR: negative materials are invalid'
        stop
      endif

      ! Direction:
      write(*,'(a)') 'Direction vector (u,v,w):'
      read(*,*) usrc,vsrc,wsrc
      write(*,'(3(1x,1pe12.5))') usrc,vsrc,wsrc
      write(*,'(a)') 'Aperture (deg) and cos():'
      read(*,*) theta
      cossrc = cos(theta*pi/180.0d0)
      write(*,'(2(1x,1pe12.5))') theta,cossrc
      if (theta.lt.0.0d0.or.theta.gt.180.0d0) then
        write(*,*) 'BIGSinisrc:ERROR: aperture must be in [0,180].'
        stop
      endif
      norm = sqrt(usrc**2+vsrc**2+wsrc**2)
      if (norm.lt.zero) then
        if (theta.ne.180.0d0) then
          write(*,*)
     &  'BIGSinisrc:ERROR: null direction only valid if aperture=180'
          stop
        endif
        usrc = 0.0d0
        vsrc = 0.0d0
        wsrc = 1.0d0
      else
        usrc = usrc/norm
        vsrc = vsrc/norm
        wsrc = wsrc/norm
      endif

      ! Init performance vars:
      shots = 0.0d0
      warned = .false.

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'BIGSinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> BIGS source initialization finished >>>>'
      end


      integer function seeki(x,xc,n)
!*******************************************************************
!*    Finds the interval (x(i),x(i+1)] containing the value xc.    *
!*                                                                 *
!*    Input:                                                       *
!*      x(1..n) -> data array                                      *
!*      xc -> point to be located                                  *
!*      n -> no. of data points                                    *
!*    Output:                                                      *
!*      index i of the semiopen interval where xc lies             *
!*    Comments:                                                    *
!*      -> If xc=x(1) then i=1 is returned.                        *
!*      -> If xc is outside the closed interval [x(1),x(n)]  the   *
!*         execution is aborted.                                   *
!*******************************************************************
      implicit none
      integer n
      real*8 xc,x(n)

      integer itop,imid

      if(xc.gt.x(n)) then
        write(*,*) 'seeki error: value outside range, xc>x(n):'
        write(*,*) xc,x(n)
        stop
      endif
      if(xc.lt.x(1)) then
        write(*,*) 'seeki error: value outside range, xc<x(1):'
        write(*,*) xc,x(1)
        stop
      endif

      seeki = 1
      itop = n
      do
        imid = (seeki+itop)/2
        if(xc.gt.x(imid)) then
          seeki = imid
        else
          itop = imid
        endif
        if(itop-seeki.le.1) exit
      enddo
      end


      subroutine rotate(u,v,w,costh,phi)
!*******************************************************************
!*    Rotates a vector; the rotation is specified by giving        *
!*    the polar and azimuthal angles in the "self-frame", as       *
!*    determined by the vector to be rotated.                      *
!*                                                                 *
!*    Input:                                                       *
!*      (u,v,w) -> input vector (=d) in the lab. frame             *
!*      costh -> cos(theta), angle between d before and after turn *
!*      phi -> azimuthal angle (rad) turned by d in its self-frame *
!*    Output:                                                      *
!*      (u,v,w) -> rotated vector components in the lab. frame     *
!*    Comments:                                                    *
!*      -> (u,v,w) should have norm=1 on input; if not, it is      *
!*         renormalized on output, provided norm>0.                *
!*      -> The algorithm is based on considering the turned vector *
!*         d' expressed in the self-frame S',                      *
!*           d' = (sin(th)cos(ph), sin(th)sin(ph), cos(th))        *
!*         and then apply a change of frame from S' to the lab     *
!*         frame. S' is defined as having its z' axis coincident   *
!*         with d, its y' axis perpendicular to z and z' and its   *
!*         x' axis equal to y'*z'. The matrix of the change is then*
!*                   / uv/rho    -v/rho    u \                     *
!*          S ->lab: | vw/rho     u/rho    v |  , rho=(u^2+v^2)^0.5*
!*                   \ -rho       0        w /                     *
!*      -> When rho=0 (w=1 or -1) z and z' are parallel and the y' *
!*         axis cannot be defined in this way. Instead y' is set to*
!*         y and therefore either x'=x (if w=1) or x'=-x (w=-1)    *
!*******************************************************************
      implicit none
      real*8 u,v,w,costh,phi

      real*8 rho2,sinphi,cosphi,sthrho,urho,vrho,sinth,norm
      real*8 szero,zero
      parameter (szero=1.0d-14,zero=1.0d-30)

      rho2 = u*u+v*v
      norm = rho2+w*w
      ! Check normalization:
      if (dabs(norm-1.0d0).gt.szero) then
        if (norm.lt.zero) then
          write(*,*)
     &      'rotate:ERROR: null vector cannot be renormalized'
          stop
        endif
        ! Renormalize:
        norm = 1.d0/dsqrt(norm)
        u = u*norm
        v = v*norm
        w = w*norm
      endif

      sinphi = dsin(phi)
      cosphi = dcos(phi)
      ! Case z' not= z:
      if (rho2.gt.zero) then
        sthrho = dsqrt((1.0d0-costh*costh)/rho2)
        urho =  u*sthrho
        vrho =  v*sthrho
        u = u*costh - vrho*sinphi +      w*urho*cosphi
        v = v*costh + urho*sinphi +      w*vrho*cosphi
        w = w*costh               - rho2*sthrho*cosphi
      else
        ! 2 especial cases when z'=z or z'=-z:
        sinth = dsqrt(1.0d0-costh*costh)
        v = sinth*sinphi
        if (w.gt.0.d0) then
          u = sinth*cosphi
          w = costh
        else
          u = -sinth*cosphi
          w = -costh
        endif
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

