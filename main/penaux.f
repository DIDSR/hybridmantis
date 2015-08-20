!*******************************************************************
!*                           PENAUX                                *
!*                                                                 *
!* Short description:                                              *
!*   Miscelaneous routines that serve to initialize and monitor a  *
!*   PENELOPE main program.                                        *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   -> common /CSIMPA/ (includes MAXMAT, max number of materials) *
!*   -> routines PEINIT and PANR                                   *
!*   from PENGEOM:                                                 *
!*   -> routine GEOMIN                                             *
!*   -> common /QTREE/                                             *
!*   from PENVARED                                                 *
!*   -> common /CFORCE/                                            *
!*   from TIMING:                                                  *
!*   -> routines CPUTIME and REALTIME                              *
!*   from other penEasy files:                                     *
!*   -> routines READVOX, INIVOX, AUTOGEO                          *
!*   -> routines TALLY, TALLYREPORT                                *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2006                                                          *
!*******************************************************************


      subroutine iniconfig
!*******************************************************************
!*    Reads simulation configuration.                              *
!*                                                                 *
!*    Output:                                                      *
!*      /ctrsim/ /rseed/                                           *
!*******************************************************************
      implicit none
      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      real*8 nhist,time0,atime,refresh,ncalls,lastime
      common /ctrsim/ nhist,time0,atime,refresh,ncalls,lastime
      character*80 buffer
      character*(*) secid,eos
      parameter (secid='[SECTION CONFIG v.2006-08-01]')
      parameter (eos='[END OF CONFIG SECTION]')
      integer in,finduf
      real*8 maxn
      parameter (maxn=1.0d15)

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'iniconfig:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'No. of histories requested: '
      read(*,*) nhist
      write(*,'(1pe12.5)') nhist
      if (nhist.gt.maxn) then
        write(*,*) ''
        write(*,'(a,1pe12.5)')
     &    'iniconfig:ERROR: No. of requested histories exceeds ',maxn
        stop
      endif
      write(*,'(a)') 'Allotted time (s):'
      read(*,*) atime
      write(*,'(1pe12.5)') atime
      write(*,'(a)') 'Update interval:'
      read(*,*) refresh
      write(*,'(1pe12.5)') refresh
      if (refresh.gt.5.0d4) then
        write(*,'(a)') 'iniconfig:ERROR: update time must be < 50000'
        stop
      endif
      write(*,'(a)') 'Random seeds:'
      read(*,*) seed1,seed2
      if (seed1.eq.0.and.seed2.eq.0) then
        ! Reads seeds from an external file:
        write(*,'(a)') '  (reading seeds for external file):'
        read(*,'(1x,a30)') buffer
        write(*,'(2x,a)') buffer
        write(*,*) ''
        in = finduf()
        open(in,file=buffer,status='old')
        read(in,*) seed1,seed2
        close(in)
      endif
      write(*,'(2(1x,i0))') seed1,seed2
      ncalls = 0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'iniconfig:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> Config initialization finished >>>>'
      end


      subroutine inigeo(nmat)
!*******************************************************************
!*    Reads geometry section and initializes PENGEOM and penVox.   *
!*                                                                 *
!*    Output:                                                      *
!*      nmat -> largest material index in GEO & VOX files          *
!*******************************************************************
      implicit none
      integer*4 nmat

      logical isQuad
      common /geoquad/ isQuad
      character*80 buffer,filen
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION PENGEOM+PENVOX v.2008-06-01]')
      parameter (eos='[END OF GEO SECTION]')
      integer ufile,finduf,error
      integer*4 maxmat,nmatvox,nmatquad,nbody
      parameter (maxmat=10)
      real*8 parinp

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'inigeo:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Read quadrics file, if there is one:
      write(*,*) ''
      write(*,'(a)')
     &  '>>>> Initializing quadric geometry >>>>'
      read(*,'(1x,a30)') filen   ! Quadrics filename
      if (len_trim(filen).ne.0) then
        isQuad = .true.
        write(*,'(a)') 'Opening quadric geometry file:'
        write(*,'(1x,a30)') filen
      else  ! No quadrics
        isQuad = .false.
        filen = 'default.tmp.geo'
        write(*,'(a,a,a)')
     &    '(No quadrics; an empty quadric geometry '//
     &    'will be automatically created in file ',
     &    trim(filen),')'
        call autogeo(filen)
      endif

      ! Init quadric geom (PENGEOM):
      ufile = finduf()
      open(ufile,file=filen,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)') 'inigeo:ERROR: unable to open quadrics file'
        stop
      endif
      write(*,'(a)') 'Now calling PENGEOM...'
      write(*,*) ''
      call geomin(parinp,0,nmatquad,nbody,ufile,6)
      close(ufile)
      write(*,*) ''
      write(*,'(a)')
     &  '>>>> Quadric geometry initialization finished >>>>'

      ! Init voxel geom (penVox), if there is one:
      call inivox(nmatvox)
      ! Note: some additional init procedures are performed by VDD tally

      ! Set no. of materials in geometry files:
      nmat = max(nmatvox,nmatquad)
      write(*,*) ''
      write(*,'(a)') 'No. of materials referenced in the geometry:'
      write(*,'(1x,i0)') nmat
      if (nmat.gt.maxmat) then
        write(*,'(a)')
     &    'inigeo:ERROR: too many materials; '//
     &    'enlarge MAXMAT in Fortran code and recompile.'
        stop
      else if (nmat.eq.0) then
        write(*,'(a)') 'inigeo:ERROR: no geometry defined.'
        stop
      endif

      ! Check section integrity:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)') 'inigeo:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,*) ''
      write(*,'(a)') '>>>> Geometry initialization finished >>>>'
      end


      subroutine inipen(emax,nmatgeo)
!*******************************************************************
!*    Reads simulation parameters and initializes PENELOPE.        *
!*                                                                 *
!*    Input:                                                       *
!*      emax -> max source energy (eV)                             *
!*      nmatgeo -> max number of materials in GEO & VOX files      *
!*******************************************************************
      implicit none
      integer*4 nmatgeo
      real*8 emax

      integer*4 maxmat
      parameter (maxmat=10)
      real*8 eabs,c1,c2,wcc,wcr
      common/csimpa/eabs(3,maxmat),c1(maxmat),c2(maxmat),
     &              wcc(maxmat),wcr(maxmat)
      real*8 dsmaxval
      common /ctrsi1/ dsmaxval(maxmat)
      character*80 buffer,mfilen
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION PENELOPE v.2008-02-20]')
      parameter (eos='[END OF PEN SECTION]')
      integer ufile,finduf,imat,i,error
      integer*4 nmat
      real*8 infty
      parameter (infty=1.0d30)

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'inipen:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ufile = finduf()
      read(*,'(1x,a30)') mfilen  ! Material file name
      ! Read simulation parameters:
      write(*,*) ''
      write(*,'(a)') 'Simulation parameters'
      write(*,'(a)') 'No. materials expected in MAT file (0=auto):'
      read(*,*) nmat
      write(*,'(1x,i0)') nmat
      if (nmat.gt.maxmat) then
        write(*,'(a)')
     &    'inipen:ERROR: Too many materials; enlarge MAXMAT'
        stop
      endif
      write(*,'(a86)') 'MAT  EABS(e-)   EABS(ph)   EABS(e+)   C1'//
     &           '         C2        WCC         WCR        DSMAX'
      if (nmat.eq.0) then
        write(*,'(a)') '(default values requested)'
        nmat = nmatgeo
        do i=1,nmat
          eabs(1,i) = min(max(emax*1.0d-2,50.0d0),1.0d6)
          eabs(2,i) = min(max(emax*1.0d-3,50.0d0),1.0d6)
          eabs(3,i) = eabs(1,i)
          c1(i)     = 0.1d0
          c2(i)     = 0.1d0
          wcc(i)    = min(emax*1.0d-2,eabs(1,i))
          wcr(i)    = min(emax*1.0d-3,eabs(2,i))
          dsmaxval(i)  = +infty
          write(*,'(i3,8(1x,1pe10.3))') i,eabs(1,i),eabs(2,i),
     &      eabs(3,i),c1(i),c2(i),wcc(i),wcr(i),dsmaxval(i)
        enddo
      else
        ! User has given simulation parameters explicitly:
        if (nmatgeo.gt.nmat) then
          write(*,'(a)')
     &      'inipen:ERROR: There are more materials in the geometry'
          write(*,'(a)')
     &      '   file than expected in the materials file'
          stop
        endif
        read(*,'(a80)') buffer
        do i=1,nmat
          read(*,*) imat,eabs(1,i),eabs(2,i),eabs(3,i),
     &              c1(i),c2(i),wcc(i),wcr(i),dsmaxval(i)
          write(*,'(i2,8(1x,1pe10.3))') imat,eabs(1,i),eabs(2,i),
     &              eabs(3,i),c1(i),c2(i),wcc(i),wcr(i),dsmaxval(i)
          if (imat.ne.i) then
            write(*,'(a)')
     &      'inipen:ERROR: Materials must be ordered sequentially:'
            write(*,'(i3)') imat
            stop
          endif
        enddo
      endif

      ! Init PENELOPE kernel:
      write(*,'(a)') 'Opening material data file:'
      write(*,'(1x,a30)') mfilen
      open(ufile,file=mfilen,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)') 'inipen:ERROR: unable to open materials file'
        stop
      endif
      write(*,'(a)') 'Now calling PEINIT...'
      write(*,*) ''
      call peinit(emax*1.01d0,nmat,ufile,6,1)  ! Add 1% to EMAX to allow perfect absorbents
      close(ufile)

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)') 'inipen:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,*) ''
      write(*,'(a)') '>>>> PENELOPE initialization finished >>>>'
      end


      logical function endsim(n)
!*******************************************************************
!*    Determines whether the simulation is done or not and writes  *
!*    progress reports.                                            *
!*                                                                 *
!*    Input:                                                       *
!*      n -> current history no.                                   *
!*    Output:                                                      *
!*      .true. if done                                             *
!*******************************************************************
      implicit none
      real*8 n

      real*8 nhist,time0,atime,refresh,ncalls,lastime
      common /ctrsim/ nhist,time0,atime,refresh,ncalls,lastime
      logical needfresh
      integer unc
      real*8 realtime,cputime,nowreal,nowcpu

      nowreal = realtime() ! Ensures one call per day at least
      nowcpu  = cputime()
      endsim  = .false.

      ! Check no. of histories:
      if (n-nhist.gt.-0.5d0) endsim = .true.

      ! Check time:
      if (atime.gt.0.0d0) then
        if (nowreal.gt.atime) endsim = .true.  ! Real time
      else
        if (nowcpu.gt.-atime) endsim = .true.  ! User (CPU) time
      endif

      if (endsim) return

      ! Check if an update is needed:
      needfresh = .false.
      if (refresh.gt.0.0d0) then
        if (nowreal-lastime.gt.refresh) then
          lastime = nowreal
          needfresh = .true.
        endif
      else
        ncalls = ncalls+1
        if (ncalls.ge.-refresh) then
          ncalls = 0
          needfresh = .true.
        endif
      endif

      ! Update:
      if (needfresh) then
        call comand(n)
        call tallyreport(n,nowcpu,unc)  ! Writes progress report
        endsim = unc.gt.1               ! Required accuracy attained?
      endif
      end


      subroutine noverflow(n)
!*******************************************************************
!*    Determines whether the number of histories overflows the     *
!*    real*8 counter; if it does, if forces the simulation to stop.*
!*                                                                 *
!*    Input:                                                       *
!*      n -> current history no.                                   *
!*******************************************************************
      implicit none
      real*8 n

      real*8 maxn
      parameter (maxn=1.0d15)

      if (n.gt.maxn) then
        write(*,*) ''
        write(*,'(a)') '***************************************** '
        write(*,'(a,f20.0)')
     &  'noverflow:WARNING: No. of histories too large for'
        write(*,'(a)')
     &  '  penEasy real*8 counters to handle.'
        write(*,'(a)')
     &  '  Forcing the simulation to stop after current history.'
        write(*,'(a)') '***************************************** '
        write(*,*) ''
        call dostop
      endif
      end


      subroutine dostop
!*******************************************************************
!*    Forces the simulation to stop after completion of current    *
!*    history by re-setting the no. of histories to simulate to    *
!*    zero.                                                        *
!*******************************************************************
      implicit none
      real*8 nhist,time0,atime,refresh,ncalls,lastime
      common /ctrsim/ nhist,time0,atime,refresh,ncalls,lastime

      nhist = 0
      end


      real*8 function dsmax()
!*******************************************************************
!*    Maximum step length as passed to JUMP; this function is an   *
!*    interface to common /CTRSI1/                                 *
!*                                                                 *
!*    Input:                                                       *
!*      /TRACK/                                                    *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 maxmat
      parameter (maxmat=10)
      real*8 dsmaxval
      common /ctrsi1/ dsmaxval(maxmat)

      dsmax = dsmaxval(mat)
      end


      logical function absorb()
!*******************************************************************
!*    Checks whether or not a particle is absorbed given its       *
!*    current energy.                                              *
!*                                                                 *
!*    Input:                                                       *
!*      /TRACK/                                                    *
!*    Output:                                                      *
!*      .true. if absorbed, .false. else                           *
!*    Comments:                                                    *
!*      -> Tallies deposited energy when necessary; to do so it    *
!*         calls TALLY.                                            *
!*      -> Absorbs particles above 1 GeV; these can be created by  *
!*         SOURCE or be the result of the annihilation of a e+.    *
!*      -> This function is well suited to implement range         *
!*         rejection by making the returned value dependent        *
!*         on the particle's position, material, energy, etc.      *
!*         Notice that in order to preserve the calculation of     *
!*         the fluence unbiased an electron or positron should     *
!*         *never* be rejected *inside* the detector.              *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 maxmat
      parameter (maxmat=10)
      real*8 eabs,c1,c2,wcc,wcr
      common/csimpa/eabs(3,maxmat),c1(maxmat),c2(maxmat),
     &              wcc(maxmat),wcr(maxmat)
      real*8 mc2,twomc2
      parameter (mc2=5.10998918d5,twomc2=2.0d0*mc2)  ! Electron rest energy

      absorb = .false.
      if (mat.eq.0) then  ! vacuum
        absorb = .true.
        return
      endif
      if (e.gt.eabs(kpar,mat).and.e.lt.1.0d9) return  ! Don't absorb

      absorb = .true.
      call tally(-98,e)         ! Tallies remaining kinetic energy
      if (kpar.eq.3.and.e.gt.0.0d0) then  ! Precaution in case of positron
        call panar              ! Puts two annihilation photons in the stack
        call tally(-97,twomc2)  ! Tally the extra 2mc^2
      endif
      if (e.gt.1.0d9) then
        write(*,*) ''
        write(*,'(a)') '***************'
        write(*,'(a,1pe12.5,a)')
     &  'absorb:WARNING: particle with energy out of range: ',e,' eV'
        write(*,'(a)') '  (it has been locally absorbed)'
        write(*,'(a)') '***************'
        write(*,*) ''
      endif
      end


      subroutine comand(n)
!*******************************************************************
!*    Reads and executes commands from an external file, allowing  *
!*    in-flight steering of the simulation.                        *
!*                                                                 *
!*    Input:                                                       *
!*      n -> current no. of histories already simulated.           *
!*      -> Command is read from external file command.in           *
!*    Output:                                                      *
!*      -> Command is executed and file command.in reset to default*
!*******************************************************************
      implicit none
      real*8 n

      real*8 nhist,time0,atime,refresh,ncalls,lastime
      common /ctrsim/ nhist,time0,atime,refresh,ncalls,lastime
      character*80 buffer
      logical ffound
      integer com,in,finduf,error
      real*8 param,maxn
      parameter (maxn=1.0d15)

      inquire(file='command.in',exist=ffound)
      if (.not.ffound) return  ! File not found
      in = finduf()
      open(in,file='command.in',iostat=error)
      if (error.ne.0) return

      ! Read command file:
      iter: do
        read(in,*,iostat=error) com
        if (error.ne.0.or.com.eq.0) exit iter
        write(*,*) ''
        write(*,'(a,f18.0)')
     &   'comand:WARNING: command received when history number was: ',
     &    n
        select case (com)
        case (1)
          read(in,*,iostat=error) param
          if (error.ne.0) exit iter
          nhist = param
          if (nhist.gt.maxn) then
            write(*,*) ''
            write(*,'(a,f18.0)')
     &      'comand:WARNING: requested max No. of histories exceeds ',
     &        maxn
            nhist = maxn
          endif
          write(*,'(a,f18.0)') '  max No. of histories reset to: ',
     &      nhist
        case (2)
          read(in,*,iostat=error) param
          if (error.ne.0) exit iter
          atime = param
          write(*,'(a,es12.5)') '  max time reset to: ',atime
        case (3)
          read(in,*,iostat=error) param
          if (error.ne.0) exit iter
          refresh = param
          write(*,'(a,es12.5)') '  update interval reset to: ',refresh
        case (9)
          read(in,'(a80)',iostat=error) buffer
          if (error.ne.0) exit iter
          write(*,'(a)') '  message from command.in follows:'
          write(*,'(a80)') buffer
        case default
          write(*,'(a,i0,a)')
     &     '  unknown command: ',com,'; remaining commands ignored.'
          exit iter
        end select
      enddo iter
      close(in,iostat=error)

      ! Reset command file:
      open(in,file='command.in',iostat=error)
      if (error.ne.0) return
      write(in,'(a)',iostat=error)
     &  ' 0   <- write here the command code'
      write(in,'(a)',iostat=error)
     &  '     <- write here the command parameter, if any'
      write(in,'(a)',iostat=error)
     &  '     <- repeat previous lines as many times as needed'
      write(in,'(a)',iostat=error)
     &  ''
      write(in,'(a)',iostat=error)
     &  '>>>> END OF INPUT >>>>'
      write(in,'(a)',iostat=error)
     &  ''
      write(in,'(a)',iostat=error)
     &  'Code  Parameter  What it does'
      write(in,'(a)',iostat=error)
     &  '-----------------------------'
      write(in,'(a)',iostat=error)
     &  '0     None       Nothing, keep going'
      write(in,'(a)',iostat=error)
     &  '1     N          Reset No. of histories to N'
      write(in,'(a)',iostat=error)
     &  '2     t          Reset simulation time to t'
      write(in,'(a)',iostat=error)
     &  '3     U          Reset update interval to U'
      write(in,'(a)',iostat=error)
     &  '9     <string>   Write <string> (<80 chars) to output'
      write(in,'(a)',iostat=error)
     &  ''
      write(in,'(a)',iostat=error)
     &  '>>>> END OF FILE >>>>'
      write(in,'(a)',iostat=error)
     &  ''
      close(in,iostat=error)
      end


      integer function finduf()
!*******************************************************************
!*    Finds a valid (i.e. unused) file unit.                       *
!*******************************************************************
      implicit none
      logical used
      integer maxuf
      parameter (maxuf=17)

      finduf = 6
      do
        finduf = finduf+1
        if (finduf.gt.maxuf) then
          write(*,*) ''
          write(*,'(a)') 'finduf:ERROR: Unable to find a valid unit'
          stop
        endif
        inquire(finduf,opened=used)
        if (.not.used) return
      enddo
      end


      subroutine getline(buffer)
!*******************************************************************
!*    Reads a new line from the keyboard. The line is returned     *
!*    only if it is not blank or a comment line.                   *
!*                                                                 *
!*    Output:                                                      *
!*      buffer -> line read.                                       *
!*    Comments:                                                    *
!*      -> comment lines start with a '#'                          *
!*******************************************************************
      implicit none
      character*80 buffer
      integer i

      do
        read(*,'(a80)') buffer
        if (buffer(1:1).eq.'#') cycle  ! A comment line
        do i=1,len(buffer)
          if (buffer(i:i).ne.' ') return
        enddo
      enddo
      end


      subroutine iniforce(e0)
!*******************************************************************
!*    Initializes the interaction forcing routines                 *
!*                                                                 *
!*    Input:                                                       *
!*      e0 -> source energy (eV)                                   *
!*******************************************************************
      implicit none
      real*8 e0

      integer ns,nb,nx
      parameter (ns=10000,nb=5000,nx=250)
      integer*4 nbody,mater,kmoth,kdght,ksurf,kflag,kalias,kslast
      common/qtree/nbody,mater(nb),kmoth(nb),kdght(nb,nx),
     1    ksurf(nb,nx),kflag(nb,nx),kalias(ns),kslast
      real*8 force
      common/cforce/force(nb,3,8)
      integer*4 maxmat
      parameter (maxmat=10)
      logical analog,active
      real*8 minwght,maxwght
      common /frc001/ analog(maxmat,3),minwght,maxwght,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION INTERACTION FORCING v.2008-05-15]')
      parameter (eos='[END OF VRIF SECTION]')
      character*80 buffer
      integer i,j,mat,kpar,icol,error
      real*8 forcing,phmfp,hmfp

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'iniforce:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(1x,a3)') buffer
      if (adjustl(buffer(1:3)).eq.'ON') then
        active = .true.
        write(*, '(a)')
     &    '>>>> Interaction Forcing is ON >>>>'
      else if (buffer(1:3).eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Interaction Forcing is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'iniforce:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'iniforce:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      write(*,'(a)') 'Weight window [Wmin,Wmax]:'
      read(*,*) minwght,maxwght
      write(*,'(2(1x,es12.5))') minwght,maxwght

      ! Reset all:
      do icol=1,8
        do kpar=1,3
          do i=1,nb
            force(i,kpar,icol) = 1.0d0
          enddo
        enddo
      enddo
      do kpar=1,3
        do mat=1,maxmat
          analog(mat,kpar) = .true.
        enddo
      enddo

      write(*,'(a)') 'MAT : KPAR : ICOL : forcing'
      read(*,'(a80)') buffer  ! Table header
      do
        read(*,*,iostat=error) mat,kpar,icol,forcing
        if (error.ne.0) then
          write(*,'(a)')
     &     'iniforce:ERROR: unable to read line containing:'
          write(*,'(a)')
     &     '  MAT,KPAR,ICOL,forcing'
          write(*,'(a)')
     &     '  Recall to use a line with MAT=0 to end the list.'
          stop
        endif
        if (mat.eq.0) exit  ! End of list
        write(*,'(3(1x,i3),1x,es12.5)')
     &    mat,kpar,icol,forcing
        if (mat.lt.0) then
          write(*,'(a)') 'iniforce:ERROR: MAT must be positive'
          stop
        endif
        if (kpar.lt.1.or.kpar.gt.3) then
          write(*,'(a)') 'iniforce:ERROR: KPAR must be in [1,3]'
          stop
        endif
        if (icol.lt.0.or.icol.gt.8) then
          write(*,'(a)') 'iniforce:ERROR: ICOL must be in [0,8]'
          stop
        endif
        if (forcing.lt.1.0d0) then
          write(*,'(a)') 'iniforce:ERROR: FORCING must not be < 1'
          stop
        endif
        do i=1,nb  ! Sets forcing for all bodies of the selected material
          if (mater(i).eq.mat) then
            if (icol.ne.0) then
              force(i,kpar,icol) = forcing
            else  ! Force all interactions
              do j=1,8
                force(i,kpar,j) = forcing
              enddo
            endif
          endif
        enddo
        analog(mat,kpar) = .false.
      enddo

      write(*,*) ''
      write(*,'(a,es12.5)')
     &  'INFO: Unforced (analog) hard mean free paths at E(eV) =',e0
      write(*,'(a)') 'MAT : KPAR : ICOL : HMFP(cm)'
      do mat=1,maxmat
        do kpar=1,3
          if (analog(mat,kpar)) cycle  ! Skip unforced interactions
          do icol=1,8
            if (kpar.eq.1.and.(icol.lt.2.or.icol.gt.5)) cycle
            if (kpar.eq.2.and.icol.gt.4) cycle
            if (kpar.eq.3.and.(icol.lt.2.or.icol.gt.6)) cycle
            hmfp = phmfp(e0,kpar,mat,icol)
            write(*,'(i3,1x,i1,1x,i1,1x,1pe12.5)')
     &        mat,kpar,icol,hmfp
          enddo
        enddo
      enddo
      write(*,*) ''
      write(*,'(a)') 'iniforce:WARNING:'//
     &  ' interaction forcing may bias pulse height spectra'

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'iniforce:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)')
     &  '>>>> Interaction forcing initialization finished >>>>'
      end


      logical function isforcing()
!*******************************************************************
!*    Determines whether interaction forcing is to be applied      *
!*******************************************************************
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 maxmat
      parameter (maxmat=10)
      logical analog,active
      real*8 minwght,maxwght
      common /frc001/ analog(maxmat,3),minwght,maxwght,active

      isforcing = .false.
      if (.not.active.or.mat.eq.0) return
      if (analog(mat,kpar).or.wght.lt.minwght.or.wght.gt.maxwght)
     &  return
      isforcing = .true.
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

