CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2006)                                     C
C  Copyright (c) 2001-2006                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C              SUBROUTINE PACKAGE PENFIELD   (version 2006)
C  *********************************************************************
C
C                         Francesc Salvat, Josep Sempau. January, 2003.
C
C  This subroutine package generates electron/positron trajectories in
C  an external static electromagnetic (EM) field. The material struc-
C  ture is described by means of the geometry package PENGEOM.
C
C  The user must provide the following subroutines (see the examples at
C  the end of this file);
C     SUBROUTINE GETEMF(X,Y,Z,EX,EY,EZ,BX,BY,BZ)
C     SUBROUTINE GETEMP(X,Y,Z,PHI)
C  which deliver the electric field (EX,EY,EZ) in V/cm, the magnetic
C  field (BX,BY,BZ) in gauss and the scalar potential PHI in volts
C  at an arbitrary point with coordinates (X,Y,Z).
C
C  The sequence of calls to generate a particle track is the following,
C  1  CALL START  ! Initializes the simulation of the track.
C
C  2  CALL TPEMF0(ULDV,ULDE,ULEM,DSMAX)  ! Determines the maximum
C       ! step length DSMAX consistent with the adopted delta values.
C  3  CALL JUMP(DSMAX,DS)  ! Determines the path length DS to the next
C       ! interaction event.
C  4  CALL TPEMF1(DS,DSEF,NCROSS)  ! Moves the particle to the
C       ! end of the step and determines its final direction and energy.
C
C  5  CALL KNOCK(DE,ICOL)  ! Simulates next interaction event.
C  6  GO TO 1 OR 2
C  This sequence has to be discontinued when the particle leaves the
C  material system or is absorbed.
C
C  Although the transport of photons is not affected by EM fields, the
C  present routines can also be called to generate photon histories.
C  This is useful to simplify the structure of the main program (all
C  kinds of particle can be tracked by using the same sequence of calls
C  to the simulation routines).
C
C  IMPORTANT NOTE: One of the practical features of PENGEOM is that
C  particles are automatically transported through void regions (i.e.
C  those with MAT=0) following straight trajectories without loosing
C  energy. Of course, in the presence of EM fields this would give
C  erroneous results. The easiest method to avoid this conflict is to
C  fill up the voids with a fake material of very small density.
C
C  *********************************************************************
C                       SUBROUTINE TPEMF1
C  *********************************************************************
      SUBROUTINE TPEMF1(DS,DSEF,NCROSS)
C
C  This subroutine generates electron trajectory segments in static uni-
C  form electromagnetic (EM) fields. It is assumed that subroutine
C  TPEMF0 has been called immediately before entering TPEMF1.
C
C  Input/output (through common /TRACK/):
C     X, Y, Z ... position coordinates,
C     U, V, W ... direction cosines of the direction vector,
C     IBODY ..... body where the point (X,Y,Z) is located,
C     MAT ....... material in body IBODY.
C
C  Input values (arguments):
C     DS ........ path length to travel,
C
C  Output values (arguments):
C     DSEF....... travelled curved path length before leaving the
C                 initial material or completing the jump (less
C                 than DS if the track crosses an interface),
C     NCROSS .... = 0 if the whole step is contained in the
C                     initial material,
C                 .gt.0 if the particle has crossed an interface,
C                       i.e. if it has entered a new material.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=2.99792458D10)  ! Speed of light (m/s)
      PARAMETER (SL2=SL*SL)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/COMFLD/CHARGE,ETOT,BETA2,CONSE,BETA,SPEED,
     1              EX,EY,EZ,BX,BY,BZ,BF2,EDOTV
C
      IF(KPAR.EQ.2) THEN
        CALL STEP(DS,DSEF,NCROSS)
        RETURN
      ENDIF
C
C  ****  Displacement due to inertia and electric force.
C
      CONSS=0.5D0*DS*DS*CONSE
      DRX=DS*U+CONSS*(EX-BETA2*U*EDOTV)
      DRY=DS*V+CONSS*(EY-BETA2*V*EDOTV)
      DRZ=DS*W+CONSS*(EZ-BETA2*W*EDOTV)
      DVX=DS*CONSE*(EX-U*EDOTV)
      DVY=DS*CONSE*(EY-V*EDOTV)
      DVZ=DS*CONSE*(EZ-W*EDOTV)
C
C  ****  Displacement due to magnetic field.
C
      IF(BF2.GT.1.0D-32) THEN
        B=DSQRT(BF2)
        OMEGA=DABS(CHARGE)*SL2*1.0D-8*B/ETOT
        OX=BX/B
        OY=BY/B
        OZ=BZ/B
        IF(CHARGE.GT.0.0D0) THEN
          OX=-OX
          OY=-OY
          OZ=-OZ
        ENDIF
        VDOTO=U*OX+V*OY+W*OZ
        V0PX=U-VDOTO*OX
        V0PY=V-VDOTO*OY
        V0PZ=W-VDOTO*OZ
        V0OX=V0PY*OZ-V0PZ*OY
        V0OY=V0PZ*OX-V0PX*OZ
        V0OZ=V0PX*OY-V0PY*OX
C
        ARG=DS*OMEGA/SPEED
        F1=DCOS(ARG)-1.0D0
        F2=DSIN(ARG)
        FACT=SPEED/OMEGA
        DRX=DRX-DS*V0PX+(F1*V0OX+F2*V0PX)*FACT
        DRY=DRY-DS*V0PY+(F1*V0OY+F2*V0PY)*FACT
        DRZ=DRZ-DS*V0PZ+(F1*V0OZ+F2*V0PZ)*FACT
        DVX=DVX+F1*V0PX-F2*V0OX
        DVY=DVY+F1*V0PY-F2*V0OY
        DVZ=DVZ+F1*V0PZ-F2*V0OZ
      ENDIF
C
C  ****  Combined displacement through the geometrical system.
C
      X0=X
      Y0=Y
      Z0=Z
      U0=U
      V0=V
      W0=W
      DSR=DSQRT(DRX**2+DRY**2+DRZ**2)
      IF(DSR.GT.1.0D-30) THEN
        U=DRX/DSR
        V=DRY/DSR
        W=DRZ/DSR
        CALL STEP(DSR,DSEF,NCROSS)
      ELSE
        WRITE(6,1000)
 1000 FORMAT(1X,'*** Possible inconsistency in subroutine TPEMF1.')
        STOP
      ENDIF
C
C  ****  Direction cosines at the end of the substep.
C
      FACT=DSEF/DSR
      U=U0+FACT*DVX
      V=V0+FACT*DVY
      W=W0+FACT*DVZ
      DSEF=FACT*DS
      FNORM=U*U+V*V+W*W
      IF(DABS(FNORM-1.0D0).GT.1.0D-13) THEN
        FNORM=1.0D0/DSQRT(FNORM)
        U=U*FNORM
        V=V*FNORM
        W=W*FNORM
      ENDIF
C
C  ****  Energy at the end of the substep.
C
      CALL GETEMP(X0,Y0,Z0,PHI0)
      CALL GETEMP(X,Y,Z,PHI)
      E=E+CHARGE*(PHI0-PHI)
C  ****  The kinetic energy cannot be negative!
      IF(E.LT.0.1D0) E=0.1D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE TPEMF0
C  *********************************************************************
      SUBROUTINE TPEMF0(ULDV,ULDE,ULEM,T)
C
C  This subroutine determines the maximum allowed step length of elec-
C  trons/positrons in the electromagnetic field, consistent with the
C  adopted delta values.
C
C  Input values (arguments):
C     ULDV... Upper limit on the amount of deflection over the step
C             due to the em-field (0.02, typical)
C     ULDE... Upper limit on the relative energy variation over the
C             step in the em-field (0.02, typical)
C     ULEM... Upper limit on the amount of change of the em-field
C             over the step (0.02, typical)
C
C  Output values (arguments):
C     T ......... maximum allowed step length.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=2.99792458D10)  ! Speed of light (m/s)
      PARAMETER (TREV=REV+REV)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/COMFLD/CHARGE,ETOT,BETA2,CONSE,BETA,SPEED,
     1              EX,EY,EZ,BX,BY,BZ,BF2,EDOTV
C
      T=1.0D35
      IF(KPAR.EQ.2) RETURN
      CHARGE=DFLOAT(KPAR-2)
C  ****  Ensure normalization of the direction vector.
      DXYZ=U*U+V*V+W*W
      IF(DABS(DXYZ-1.0D0).GT.1.0D-14) THEN
        FNORM=1.0D0/DSQRT(DXYZ)
        U=FNORM*U
        V=FNORM*V
        W=FNORM*W
      ENDIF
C
      ETOT=E+REV
      IF(E.GT.1.0D3) THEN
        BETA2=E*(E+TREV)/ETOT**2
        CONSE=CHARGE/(ETOT*BETA2)
      ELSE
        TECOR=2.0D0*DMAX1(E,0.01D0)
        BETA2=TECOR/REV
        CONSE=CHARGE/TECOR
      ENDIF
      BETA=DSQRT(BETA2)
      SPEED=BETA*SL
      CONSB=CONSE*SPEED*1.0D-8
C
      CALL GETEMF(X,Y,Z,EX,EY,EZ,BX,BY,BZ)
      EF2=EX**2+EY**2+EZ**2
      BF2=BX**2+BY**2+BZ**2
C
C  ****  Step-size truncation based upon energy change
C        by the external field.
C
      EDOTV=EX*U+EY*V+EZ*W
      AEDOTV=DABS(EDOTV)
      IF(AEDOTV.GT.1.0D-30) T=DMIN1(T,DMAX1(ULDE*E,0.1D0)/AEDOTV)
C  ****  If the particle is at rest (approx.) it starts moving in the
C        direction of the electric field.
      IF(E.LT.0.101D0) THEN
        IF(EF2.GT.1.0D-16) THEN
          EDOTV=DSQRT(CHARGE**2*EF2)
          U=CHARGE*EX/EDOTV
          V=CHARGE*EY/EDOTV
          W=CHARGE*EZ/EDOTV
        ELSE
          STOP 'The electron came to rest with no E-field.'
        ENDIF
      ENDIF
C
C  ****  Step-size truncation based upon direction vector change
C        by the external field.
C
      DVX=CONSB*(V*BZ-W*BY)+CONSE*(EX-U*EDOTV)
      DVY=CONSB*(W*BX-U*BZ)+CONSE*(EY-V*EDOTV)
      DVZ=CONSB*(U*BY-V*BX)+CONSE*(EZ-W*EDOTV)
C
      DV2=(DVX*DVX+DVY*DVY+DVZ*DVZ)
      IF(DV2.GT.1.0D-30) T=DMIN1(T,ULDV/SQRT(DV2))
C
C  ****  Step-size truncation based upon keeping the external field
C        approximately a constant.
C
      IF(T.GT.1.0D-8) THEN
        XT=X+T*U
        YT=Y+T*V
        ZT=Z+T*W
        CALL GETEMF(XT,YT,ZT,EXT,EYT,EZT,BXT,BYT,BZT)
C
        VAREF2=(EXT-EX)**2+(EYT-EY)**2+(EZT-EZ)**2
        IF(VAREF2.GT.1.0D-6*EF2) THEN
          IF(EF2.GT.0.01D0) T=DMIN1(T,ULEM*T*DSQRT(EF2/VAREF2))
        ENDIF
C
        VARBF2=(BXT-BX)**2+(BYT-BY)**2+(BZT-BZ)**2
        IF(VARBF2.GT.1.0D-6*BF2) THEN
          IF(BF2.GT.0.01D0) T=DMIN1(T,ULEM*T*DSQRT(BF2/VARBF2))
        ENDIF
      ELSE
        T=1.0D-8
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GETEMF
C  *********************************************************************
      SUBROUTINE GETEMF(XP,YP,ZP,EX,EY,EZ,BX,BY,BZ)
C
C  This routine returns the electric and magnetic field components at
C  the input point (XP,YP,ZP).
C
C  Arguments:
C  (XP,YP,ZP) ...  Position coordinates (cm) -input.
C  (EX,EY,EZ) ...  Components of the E-field, in V/cm -output.
C  (BX,BY,BZ) ...  Components of the B-field, in gauss -output.
C
C  It is assumed that users will provide an equivalent routine, with
C  the same name and arguments, to define their fields.
C
C  For example:
C  For a uniform E-field: (EX,EY,EZ) three constants.
C                         (BX,BY,BZ)=(0,0,0).
C  For a uniform B-field: (EX,EY,EZ)=(0,0,0).
C                         (BX,BY,BZ) three constants.
C  For a point charge at rest at the origin of coordinates:
C                         EX=C*XP/(SQRT(XP**2+YP**2+ZP**2))**3,
C                         EY=C*YP/(SQRT(XP**2+YP**2+ZP**2))**3,
C                         EZ=C*ZP/(SQRT(XP**2+YP**2+ZP**2))**3,
C                         (BX,BY,BZ)=(0,0,0).
C  where C is a constant.
C
C  This example is for an arbitrary uniform EM field, the components of
C  the E-field (in V/cm) and of the B-field (in gauss) are entered
C  through the common block /UFIELD/.
C
C  Notice that 1 gauss = 1.0D-4 tesla and  1 statvolt = 300 volt.
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/UFIELD/EFX,EFY,EFZ,BFX,BFY,BFZ
              
      real*8 Zmax, Zmin                     !!DeBuG!!  EM field: Z extend of the field, passed from main to GETEMF
      COMMON/UFIELD_Z/ Zmin, Zmax      
C
C  ****  The field is set equal to zero in regions where MAT=0, i.e. in
C        the outher vacuum, since otherwise particles that leave the
C        system could return to it under the action of the field.
C
        !!DeBuG!! OLD: IF(MAT.EQ.0) THEN
      IF(MAT.EQ.0 .or. ZP.gt.Zmax .or. ZP.lt.Zmin) THEN     !!DeBuG!!  EM field: Z extend 
        EX=0.0D0
        EY=0.0D0
        EZ=0.0D0
        BX=0.0D0
        BY=0.0D0
        BZ=0.0D0
      ELSE
        EX=EFX
        EY=EFY
        EZ=EFZ
        BX=BFX
        BY=BFY
        BZ=BFZ
      ENDIF  
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GETEMP
C  *********************************************************************
      SUBROUTINE GETEMP(XP,YP,ZP,PHI)
C
C  This routine returns the electric potential at the point (X,Y,Z).
C
C  Arguments:
C  (XP,YP,ZP) ......  Position coordinates (cm) -input.
C  PHI .............  Electric potential, in V -output.
C
C  It is assumed that users will provide an equivalent routine, with
C  the same name and arguments, for their electric field.
C
C  For example:
C  For a constant E-field (EX,EY,EZ) in V/cm,
c               PHI=-(EX*XP+EY*YP+EZ*ZP).
C  For a point charge at the origin of coordinates,
C                          PHI=C/SQRT(XP**2+YP**2+ZP**2),
C  where C is a constant.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/UFIELD/EFX,EFY,EFZ,BFX,BFY,BFZ

      real*8 Zmax, Zmin                     !!DeBuG!!  EM field: Z extend of the field, passed from main to GETEMF
      COMMON/UFIELD_Z/ Zmin, Zmax
            
C
C  ****  The potential is set equal to zero in regions where MAT=0,
C        i.e. in the outher vacuum, since otherwise particles that leave
C        the system could return to it under the action of the field.
C
      !!DeBuG!!  OLD:  IF(MAT.EQ.0) THEN      
      IF(MAT.EQ.0 .or. ZP.gt.Zmax .or. ZP.lt.Zmin) THEN     !!DeBuG!!  EM field: Z extend 
        PHI=0.0D0
      ELSE
        PHI = -EFZ*ZP                                     !!DeBuG!!  EM field: Z extend
         !!DeBuG!! Original code:  PHI=-(EFX*XP+EFY*YP+EFZ*ZP)
      ENDIF  
C
      RETURN
      END
