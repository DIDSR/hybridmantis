!*******************************************************************
!*                          TALLY                                  *
!*                    PHASE SPACE FILE                             *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Writes to a file the state of all particles that reach        *
!*   a given material, considered as a particle sink.              *
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


      subroutine PSFtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      arg -> current history number (when mode=1)                *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer detmat,psfunit,formatiaea
      real*8 nhist,nhlast,npar,nele,npos
      common /scopsf/ nhist,nhlast,npar,nele,npos,detmat,psfunit,
     &                formatiaea,active
      character*80 fmtstr
      parameter (fmtstr='(i0,8(1x,es12.5),6(1x,i0))')
      integer*4 dn

      if (.not.active) return

      if (mode.eq.4.or.mode.eq.-99) then  ! New material
        if (mat.ne.detmat) return
        dn = nhist-nhlast+0.5d0
        write(psfunit,fmtstr) kpar,e,x,y,z,u,v,w,wght,dn,
     &                        ilb(1),ilb(2),ilb(3),ilb(4),ilb(5)
        nhlast = nhist
        npar = npar+1
        if (kpar.eq.1) then
          nele = nele+1
        else if (kpar.eq.3) then
          npos = npos+1
        endif
        ! mat = 0  ! Force main to stop this particle's simulation
        ! Particle is stopped by setting Eabs=+infty in the PSF detector

      else if (mode.eq.1.or.mode.eq.2) then  ! New history or hist. modified
        nhist = arg

      endif
      end


      subroutine PSFreport(n,cputim)
!*******************************************************************
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Comments:                                                    *
!*      -> 'cputim' should not include initialization procedures;  *
!*         enter 0 or neg. if not available.                       *
!*******************************************************************
      implicit none
      real*8 n,cputim

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      logical active
      integer detmat,psfunit,formatiaea
      real*8 nhist,nhlast,npar,nele,npos
      common /scopsf/ nhist,nhlast,npar,nele,npos,detmat,psfunit,
     &                formatiaea,active
      integer out,finduf,error

      if (.not.active) return

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyPhaseSpaceFile.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'PSFreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PSF]'
      write(out,'(a)')
     &'# No. of electrons/photons/positrons/total written to PSF:'
      write(out,'(f18.0)') nele
      write(out,'(f18.0)') npar-nele-npos
      write(out,'(f18.0)') npos
      write(out,'(f18.0)') npar

      ! Generic report:
      write(out,*)
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
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      close(out)
      end


      subroutine PSFinitally
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      implicit none
      logical active
      integer detmat,psfunit,formatiaea
      real*8 nhist,nhlast,npar,nele,npos
      common /scopsf/ nhist,nhlast,npar,nele,npos,detmat,psfunit,
     &                formatiaea,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PHASE SPACE FILE v.2008-06-01]')
      parameter (eos='[END OF PSF SECTION]')
      character*80 buffer
      integer finduf,error

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PSFinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Phase Space File is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PSFinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PSFinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      read(*,*) formatiaea
      if (formatiaea.eq.0) then
        write(*,'(a)') 'PSF format: standard penEasy in ASCII.'
      else if (formatiaea.eq.1) then
        write(*,'(a)') 'PSFinitally:ERROR: IAEA PSF format requested '//
     &    'but not available.'
        write(*,'(a)') '  Read the instructions provided in the '//
     &    'README file to activate this option.'
        stop
      else
        write(*,'(a)') 'PSFinitally:ERROR: PSF format must be 0 or 1.'
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) detmat
      write(*,'(i3)') detmat
      if (detmat.le.0) then
        write(*,*) 'PSFinitally:ERROR: detection material must be >0'
        stop
      endif
      write(*,'(a)') 'PSF filename:'
      read(*,'(1x,a30)') buffer
      write(*,'(a)') buffer
      psfunit = finduf()
      write(*,'(a)') 'Opening PSF as unit:'
      write(*,'(1x,i0)') psfunit
      open(psfunit,file=buffer)
      write(psfunit,'(a)')
     &  '# [PHASE SPACE FILE FORMAT penEasy v.2008-05-15]'
      write(psfunit,'(a)')
     &  '# KPAR : E : X : Y : Z : U : V : W : WGHT : '//
     &  'DeltaN : ILB(1..5)'

      !*** Init vars:
      nhist  = 0.0d0
      nhlast = 0.0d0
      npar   = 0.0d0
      nele   = 0.0d0
      npos   = 0.0d0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PSFinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PSF tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

