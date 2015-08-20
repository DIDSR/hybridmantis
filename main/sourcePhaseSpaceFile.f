!*******************************************************************
!*                         SOURCE                                  *
!*                   PHASE SPACE FILE                              *
!*                                                                 *
!* Short description:                                              *
!*   Generation of primary particle states for radiation transport *
!*   calculations with PENELOPE.                                   *
!*                                                                 *
!*   Initial particle states are read from an external phase space *
!*   file (PSF). Notice that, due to the modification of history   *
!*   numbers according to the information in the PSF, this source  *
!*   is *incompatible* with other source models.                   *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routine STORES                                             *
!*   from PENAUX.F                                                 *
!*   -> routine DOSTOP                                             *
!*   -> routine FINDUF                                             *
!*   from PENVOX.F                                                 *
!*   -> routine LOCATEX                                            *
!*   -> routine STEPX                                              *
!*   from other penEasy libraries:                                 *
!*   -> routine GETLINE,TALLY                                      *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*******************************************************************


      subroutine PSFsource(n)
!*******************************************************************
!*    Input:                                                       *
!*      n -> top history counter                                   *
!*    Output:                                                      *
!*      through /track/ and sec stack                              *
!*      n -> top history counter                                   *
!*******************************************************************
      implicit none
      real*8 n

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns
      logical active
      integer in,split,formatiaea
      real*8 rot,xshift,yshift,zshift,nlin
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                formatiaea,active
      logical getpar
      integer i
      integer*4 ncross
      real*8 infty,dsef,norm
      parameter (infty=1.0d30)

      if (.not.active) return

      n = n+dble(dns-1) ! Corrects history counter from main program
      call tally(2,n)   ! Inform tallies about change of history number
      do
        ! Load previously stored particle into active common:
        kpar = kpars
        wght = wghts/split
        e = es
        ! Rotate and translate position:
        x = xshift+rot(1,1)*xs+rot(1,2)*ys+rot(1,3)*zs
        y = yshift+rot(2,1)*xs+rot(2,2)*ys+rot(2,3)*zs
        z = zshift+rot(3,1)*xs+rot(3,2)*ys+rot(3,3)*zs
        ! Rotate direction and renormalize to double precision:
        u = rot(1,1)*us+rot(1,2)*vs+rot(1,3)*ws
        v = rot(2,1)*us+rot(2,2)*vs+rot(2,3)*ws
        w = rot(3,1)*us+rot(3,2)*vs+rot(3,3)*ws
        norm = 1.0d0/sqrt(u**2+v**2+w**2)
        u = u*norm
        v = v*norm
        w = w*norm
        ilb(1) = ilbs(1)
        ilb(2) = ilbs(2)
        ilb(3) = ilbs(3)
        ilb(4) = ilbs(4)
        ilb(5) = ilbs(5)
        call locatex                                 ! Sets ibody and mat from /track/ values
        if (mat.eq.0) call stepx(infty,dsef,ncross)  ! Where is it aiming at
        do i=1,split                                 ! Particle splitting
          call stores(e,x,y,z,u,v,w,wght,kpar,ilb)   ! Store particle in stack
          call tally(0,e)                            ! Tally new particle
        enddo

        ! Read a new particle and store for later calls:
        if (.not.getpar()) then
          call dostop                                ! Forces simulation to stop
          write(*,*) ''
          write(*,'(a)')
     &      'PSFsource: PSF exhausted; simulation stopped forcefully'
          exit
        endif
        if (dns.ne.0) exit                           ! While top primary remains the same
      enddo
      end


      subroutine PSFinisrc(activated,emax)
!*******************************************************************
!*    Initializes. To be called before SOURCE.                     *
!*                                                                 *
!*    Output:                                                      *
!*      activated -> TRUE if the source is active.                 *
!*      emax -> max source energy (eV)                             *
!*******************************************************************
      implicit none
      logical activated
      real*8 emax

      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns
      logical active
      integer in,split,formatiaea
      real*8 rot,xshift,yshift,zshift,nlin
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                formatiaea,active
      logical getpar,checkedFormat
      character*80 psfnam,buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION SOURCE PHASE SPACE FILE v.2008-06-01]')
      parameter (eos='[END OF SPSF SECTION]')
      integer finduf,error,validate
      real*8 ntop,nele,npho,npos
      real*8 omega,theta,phi,comega,ctheta,cphi,somega,stheta,sphi
      real*8 pi,deg2rad,mc2,twomc2
      parameter (pi=3.1415926535897932d0,deg2rad=pi/180.0d0)
      parameter (mc2=5.10998918d5,twomc2=2.0d0*mc2)

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PSFinisrc:ERROR: incorrect section header;'
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
     &    '>>>> Source Phase Space File is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PSFinisrc:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PSFinisrc:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      read(*,*) formatiaea
      if (formatiaea.eq.0) then
        write(*,'(a)') 'PSF format: standard penEasy in ASCII.'
      else if (formatiaea.eq.1) then
        write(*,'(a)')
     &  'PSFinisrc:ERROR: IAEA PSF format requested but not available.'
        write(*,'(a)')
     &    '  Read the instructions provided in the README file to '//
     &    'activate this option.'
        stop
      else
        write(*,'(a)') 'PSFinisrc:ERROR: PSF format must be 0 or 1.'
        stop
      endif

      write(*,'(a)') 'PSF filename:'
      read(*,'(1x,a30)') psfnam
      write(*,'(1x,a)') psfnam

      write(*,'(a)') 'Splitting factor:'
      read(*,*) split
      write(*,'(1x,i0)') split
      if (split.lt.1) then
        write(*,'(a)') 'PSFinisrc:ERROR: split < 1'
        stop
      endif

      write(*,'(a)')
     & 'Euler angles (deg) to rotate position and direction vectors:'
      read(*,*) omega,theta,phi
      write(*,'(3(1x,1pe12.5))') omega,theta,phi
      omega = omega*deg2rad
      theta = theta*deg2rad
      phi = phi*deg2rad
      ! Calculate rotation matrix:
      sphi   = sin(phi)
      cphi   = cos(phi)
      stheta = sin(theta)
      ctheta = cos(theta)
      somega = sin(omega)
      comega = cos(omega)
      rot(1,1) = cphi*ctheta*comega-sphi*somega
      rot(1,2) = -cphi*ctheta*somega-sphi*comega
      rot(1,3) = cphi*stheta
      rot(2,1) = sphi*ctheta*comega+cphi*somega
      rot(2,2) = -sphi*ctheta*somega+cphi*comega
      rot(2,3) = sphi*stheta
      rot(3,1) = -stheta*comega
      rot(3,2) = stheta*somega
      rot(3,3) = ctheta

      write(*,'(a)') 'Cartesian components of position shift (cm):'
      read(*,*) xshift,yshift,zshift
      write(*,'(3(1x,1pe12.5))') xshift,yshift,zshift

      read(*,*) validate
      if (validate.eq.1) then
        ! Pre-read PSF, validate and statistics:
        write(*,'(a)') 'Starting PSF validation'
        in = finduf()
        open(in,file=psfnam,status='old',iostat=error)
        if (error.ne.0) then
          write(*,'(a)') 'PSFinisrc:ERROR: cannot open the PSF'
          stop
        endif
        nlin = 0
        call checkFormat  ! Find if it is 2008-compliant
        checkedFormat = .true.

        emax = 0.0d0
        ntop = 0
        nele = 0
        npho = 0
        npos = 0
        do
          if (.not.getpar()) exit  ! EOF reached
          ! Count particles in PSF:
          ntop = ntop+dns  ! Note that nlin is increased by GETPAR
          if (kpars.eq.1) then
            nele = nele+1
          else if (kpars.eq.2) then
            npho = npho+1
          else if (kpars.eq.3) then
            npos = npos+1
          else
            write(*,'(a)') 'PSFinisrc:ERROR: invalid KPAR: in line:'
            write(*,'(1x,i10,1x,f18.0)') kpars,nlin
            stop
          endif
          if (es.lt.0.0d0.or.es.gt.1.0d9) then
            write(*,'(a)')
     &        'PSFinisrc:ERROR: invalid energy(eV): in line:'
            write(*,'(1pe12.5,1x,f18.0)') es,nlin
            stop
          endif
          emax = max(emax,es)
          if (us**2+vs**2+ws**2.lt.1.0d-30) then
            write(*,'(a)')
     &        'PSFinisrc:ERROR: null vector direction found in line:'
            write(*,'(f18.0)') nlin
            stop
          endif
        enddo
        close(in)

        write(*,'(a)') 'PSF statistics:'
        write(*,'(a)') '  No. electrons:'
        write(*,'(2x,f18.0)') nele
        write(*,'(a)') '  No. photons:'
        write(*,'(2x,f18.0)') npho
        write(*,'(a)') '  No. positrons:'
        write(*,'(2x,f18.0)') npos
        write(*,'(a)') '  No. particles, total:'
        write(*,'(2x,f18.0)') nele+npho+npos
        write(*,'(a)') '  No. top primary histories:'
        write(*,'(a)') '  (may be less than actual number if '//
     &    'last histories did not contribute to PSF)'
        write(*,'(2x,f18.0)') ntop
        write(*,'(a)') '  Max energy(eV):'
        write(*,'(2x,1pe12.5)') emax
        write(*,'(a)') ' (max energy declared in input file ignored)'
        if (npos.gt.0) emax = emax+twomc2  ! Allow for e+ annihilation
        read(*,'(a80)') buffer  ! Dummy line

      else  ! Do not validate the PSF
        write(*,'(a)') '** User opted not to pre-validate PSF **'
        write(*,'(a)') 'Max energy (eV) taken from input file:'
        read(*,*) emax
        write(*,'(1pe12.5)') emax
        checkedFormat = .false.
      endif

      ! Prepare for 1st call to SOURCE:
      in = finduf()
      open(in,file=psfnam,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)') 'PSFinisrc:ERROR: cannot open the PSF'
        stop
      endif
      nlin = 0
      if (.not.checkedFormat) call checkFormat
      if (.not.getpar()) then
        write(*,'(a)') 'PSFinisrc:ERROR: PSF is empty'
        stop
      endif

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'PSFinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PSF source initialization finished >>>>'
      end


      logical function getpar()
!*******************************************************************
!*    Reads a new particle from the PSF.                           *
!*                                                                 *
!*    Output:                                                      *
!*      -> returns .false. if EOF has been reached, else .true.    *
!*      -> particle state in /srcpsf/                              *
!*******************************************************************
      implicit none
      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns
      logical active
      integer in,split,formatiaea
      real*8 rot,xshift,yshift,zshift,nlin
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                formatiaea,active
      logical format2008
      common /formatVer/ format2008
      character*256 buffer
      integer error
      real*8 maxn
      parameter (maxn=1.0d15)

      do
        read(in,'(a256)',end=10,iostat=error) buffer
        if (error.ne.0) then
          write(*,'(a)')
     &      'getpar:ERROR: unable to read PSF line; last line read:'
          write(*,'(f18.0)') nlin
          stop
        endif
        if (nlin.gt.-0.5d0) nlin = nlin+1  ! Lines read, if not overflowed yet
        if (nlin.gt.maxn) then  ! Counter overflow
          write(*,*) ''
          write(*,'(a)') '*****************************************'
          write(*,'(a,f18.0)')
     &    'getpar:WARNING: No. of lines in PSF exceeds ',maxn
          write(*,'(a)')
     &    '  This is too large for penEasy real*8 counters to handle'
          write(*,'(a)')
     &    '  so, continuing without counting lines.'
          write(*,'(a)') '*****************************************'
          write(*,*) ''
          nlin = -1  ! Mark as overflowed
        endif
        if (buffer(1:1).ne.'#') exit  ! A non-comment line was found
      enddo

      if (format2008) then  ! 2008-compliant format
        read(buffer,*,iostat=error)
     &    kpars,es,xs,ys,zs,us,vs,ws,wghts,dns,
     &    ilbs(1),ilbs(2),ilbs(3),ilbs(4),ilbs(5)
      else               ! Pre-2008 format
        read(buffer,*,iostat=error)
     &    kpars,es,xs,ys,zs,us,vs,ws,wghts,ilbs(5),dns
        ilbs(1) = 1      ! Assume all particles are primaries
        ilbs(2) = 0
        ilbs(3) = 0
        ilbs(4) = 0
      endif
      if (error.ne.0) then
        write(*,'(a)')
     &    'getpar:ERROR: invalid or missing datum in PSF line:'
        write(*,'(f20.0)') nlin
        write(*,'(a)') '  line contents:'
        write(*,'(a)') buffer
        stop
      endif
      getpar = .true.
      return

 10   getpar = .false.  ! EOF
      end


      subroutine checkFormat
!*******************************************************************
!*    Checks if PSF data format is 2008-compliant.                 *
!*******************************************************************
      implicit none
      logical active
      integer in,split,formatiaea
      real*8 rot,xshift,yshift,zshift,nlin
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                formatiaea,active
      logical format2008
      common /formatVer/ format2008
      character*256 buffer,id2008,id2006
      parameter (id2008=
     & '# [PHASE SPACE FILE FORMAT penEasy v.2008-05-15]')
      parameter (id2006=
     & '# kpar : e : x : y : z : u : v : w : wght : ilb(5) : DeltaN')
      integer error

      read(in,'(a256)',iostat=error) buffer
      if (error.ne.0) then
        write(*,'(a)')
     &    'ckeckFormat:ERROR: unable to read first PSF line.'
        stop
      endif
      nlin = nlin+1

      if (buffer.eq.id2008) then
        format2008 = .true.
      else if (buffer.eq.id2006) then
        format2008 = .false.
        write(*,*) ''
        write(*,'(a)')
     &    '********************'
        write(*,'(a)')
     &    'checkFormat:WARNING: header for 2008-compliant format '//
     &    'not found;'
        write(*,'(a)')
     &    '  assuming pre-2008 data format, which consists of:'
        write(*,'(a)')
     &    '  KPAR:E:X:Y:Z:U:V:W:WGHT:ILB(5):DeltaN'
        write(*,'(a)')
     &    '********************'
        write(*,*) ''
      else
        format2008 = .false.
        write(*,*) ''
        write(*,'(a)')
     &  'checkFormat:ERROR: unable to identify PSF format;'
        write(*,'(a)')
     &  '  expecting one of these two headers in first line of PSF:'
        write(*,'(a)') id2008
        write(*,'(a)') id2006
        write(*,'(a)')
     &    '  (for 2008 and 2006 formats, respectively) '//
     &    'but neither was found.'
        stop
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

