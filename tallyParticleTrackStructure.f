!*******************************************************************
!*                          TALLY                                  *
!*                 PARTICLE TRACK STRUCTURE                        *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Writes to a file the position and value of energy loss events *
!*   so that particle tracks can be visualized.                    *
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


      subroutine PTStally(mode,arg)
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
      integer ptsunit
      integer*4 ntrack,trackmax
      real*8 xlast,ylast,zlast
      common /scopts/ xlast,ylast,zlast,ntrack,trackmax,
     &                ptsunit,active
      real*8 de

      if (.not.active.or.ntrack.gt.trackmax) return

      if (mode.eq.1) then  ! Start of a new history
        ntrack = ntrack+1
        write(ptsunit,*) ''
        write(ptsunit,'(a)') '# New history started'
        return
      endif

      if (mode.lt.0) then
        if (mode.eq.-99) then  ! A new particle starts
          if (mat.eq.0) return ! Particle not in object, no report
          de = 0.0d0
          write(ptsunit,*) ''
          write(ptsunit,*) ''  ! Separate data sets with double blank
        else
          de = arg             ! Energy scored
        endif
        xlast = x
        ylast = y
        zlast = z

      else if (mode.eq.4) then
        if (mat.eq.0) then     ! The particle escaped
          de = 0.0d0
          x = xlast+arg*u      ! Move the particle up to the exit surface
          y = ylast+arg*v
          z = zlast+arg*w
        else
          xlast = x
          ylast = y
          zlast = z
          return  ! No need to report
        endif

      else
        return    ! No need to report

      endif

      write(ptsunit,'(i0,2(1x,i0),6(1x,es9.2),1x,i0)')
     &  kpar,ibody,mat,x,y,z,e,de,wght,ilb(5)
      end


      subroutine PTSreport
!*******************************************************************
!*    Dummy routine, does nothing. PTStally does the job.          *
!*******************************************************************
      implicit none
      end


      subroutine PTSinitally
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      implicit none
      logical active
      integer ptsunit
      integer*4 ntrack,trackmax
      real*8 xlast,ylast,zlast
      common /scopts/ xlast,ylast,zlast,ntrack,trackmax,
     &                ptsunit,active
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PARTICLE TRACK STRUCTURE v.2008-05-15]')
      parameter (eos='[END OF PTS SECTION]')
      character*80 buffer
      integer finduf,error

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PTSinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Particle Track Structure is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PTSinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PTSinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ptsunit = finduf()
      buffer = 'tallyParticleTrackStructure.dat'
      open(ptsunit,file=buffer,iostat=error)
      if (error.ne.0) then
        write(*,'(a)')
     &    'PTSinitally:ERROR: cannot open track data file:'
        write(*,'(a80)') buffer
        stop
      endif
      write(ptsunit,'(a)')
     &  '# kpar : body : mat : x : y : z : e : elost : wght : ilb(5)'
      write(*,'(a)') 'Number of history tracks to display:'
      read(*,*) trackmax
      write(*,'(1x,i0)') trackmax
      if (trackmax.gt.10000) then
        write(*,'(a)')
     &  'PTSinitally:WARNING: No. tracks to report set to max=10000'
        trackmax = 10000
      endif

      ! Init vars:
      ntrack = 0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PTSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PTS tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

