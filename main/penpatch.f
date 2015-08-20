      SUBROUTINE KNOCKx(DE,ICOL)
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
C  software for any purpose. It is provided "as is" without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*******************************************************************
C*   This routine supersedes the standard KNOCK of PENELOPE.       *
C*   It provides additional information that is needed to tally    *
C*   fluences and dose-fluences. This information is accessed via  *
C*   a common block.                                               *
C*                                                                 *
C*   Notes:                                                        *
C*   -> Parts that have been changed are written in lowercase and  *
C*      marked with 'cx' (except for the this header)              *
C*******************************************************************
C
C     Simulation of random hinges and hard interaction events.
C
C  Output arguments:
C    DE ..... energy deposited by the particle in the material. It is
C             usually equal to the difference between the energies
C             before and after the interaction.
C    ICOL ... kind of interaction suffered by the particle.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI)
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
      PARAMETER (TRUNC=1.01538698D0)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
      COMMON/CHIST/ILBA(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C  ****  Bremsstrahlung emission.
      PARAMETER (NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),DPDFB(MAXMAT,NEGP,NBW),
     2  PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  DW1EL(MAXMAT,NEGP),DW2EL(MAXMAT,NEGP),
     5  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     6  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  DW1PL(MAXMAT,NEGP),DW2PL(MAXMAT,NEGP),
     5  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     6  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
      COMMON/CJUMP1/MODE,KSOFTE,KSOFTI,KDELTA
C
      COMMON/CELSEP/EELMAX(MAXMAT),PELMAX(MAXMAT),
     1              RNDCEd(MAXMAT,NEGP),RNDCPd(MAXMAT,NEGP)
C
      EXTERNAL RAND
C
cx:   *** This common us used to extract info on the restricted Eloss:
      logical called
      real*8 derest
      common /scofl0/ derest,called

cx: (default)
      derest = 0.0d0
cx: This is to test that KNOCKX has been used:
      called = .true.

      IF(KPAR.EQ.2) GO TO 2000
      IF(KPAR.EQ.3) GO TO 3000
C
C  ************  Electrons (KPAR=1).
C
C1000 CONTINUE
      IF(MODE.EQ.1) GO TO 1100
C
C  ****  Artificial soft event (ICOL=1).
C
      ICOL=1
      MODE=1
C
      IF(KSOFTI.EQ.0) THEN
        DE=0.0D0
      ELSE
        EDE0=W1*DST
        VDE0=W2*DST
        FSEDE=MAX(1.0D0-DW1EL(M,KE)*EDE0,0.75D0)
        FSVDE=MAX(1.0D0-DW2EL(M,KE)*EDE0,0.75D0)
        EDE=EDE0*FSEDE
        VDE=VDE0*FSVDE
C  ****  Generation of random values DE with mean EDE and variance VDE.
        SIGMA=SQRT(VDE)
        IF(SIGMA.LT.0.333333333D0*EDE) THEN
C  ****  Truncated Gaussian distribution.
          DE=EDE+RNDG3()*SIGMA
        ELSE
          RU=RAND(1.0D0)
          EDE2=EDE*EDE
          VDE3=3.0D0*VDE
          IF(EDE2.LT.VDE3) THEN
            PNULL=(VDE3-EDE2)/(VDE3+3.0D0*EDE2)
            IF(RU.LT.PNULL) THEN
              DE=0.0D0
            ELSE
              DE=1.5D0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0D0-PNULL)
            ENDIF
          ELSE
            DE=EDE+(2.0D0*RU-1.0D0)*SQRT(VDE3)
          ENDIF
        ENDIF
      ENDIF
C
      E=E-DE
cx:
      derest = de
      IF(E.LT.EABS(1,M)) THEN
        DE=E+DE
        E=0.0D0
        RETURN
      ENDIF
      IF(KSOFTE.EQ.0) RETURN
C
C  ****  Bielajew's randomly alternate hinge.
C
      IF(RAND(2.0D0).GT.0.5D0.AND.DE.GT.1.0D-3) THEN
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        IF(T1E(M,KE+1).GT.-78.3D0) THEN
          T1=EXP(T1E(M,KE)+(T1E(M,KE+1)-T1E(M,KE))*XEK)
          T2=EXP(T2E(M,KE)+(T2E(M,KE+1)-T2E(M,KE))*XEK)
        ELSE
          T1=0.0D0
          T2=0.0D0
        ENDIF
        IF(T1.LT.1.0D-35) RETURN
      ENDIF
C  ****  1st and 2nd moments of the angular distribution.
      EMU1=0.5D0*(1.0D0-EXP(-DST*T1))
      EMU2=EMU1-(1.0D0-EXP(-DST*T2))/6.0D0
C  ****  Sampling from a two-bar histogram with these moments.
      PNUM=2.0D0*EMU1-3.0D0*EMU2
      PDEN=1.0D0-2.0D0*EMU1
      PB=PNUM/PDEN
      PA=PDEN+PB
      RND=RAND(3.0D0)
      IF(RND.LT.PA) THEN
        CDT=1.0D0-2.0D0*PB*(RND/PA)
      ELSE
        CDT=1.0D0-2.0D0*(PB+(1.0D0-PB)*((RND-PA)/(1.0D0-PA)))
      ENDIF
      DF=TWOPI*RAND(4.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ************  Hard event.
C
 1100 CONTINUE
      MODE=0
C  ****  A delta interaction (ICOL=7) occurs when the maximum
C        allowed step length is exceeded.
      IF(KDELTA.EQ.1) THEN
        ICOL=7
        DE=0.0D0
        RETURN
      ENDIF
C  ****  Interaction probabilities.
      STNOW=P(2)+P(3)+P(4)+P(5)
C  ****  Random sampling of the interaction type.
      STS=MAX(STNOW,ST)*RAND(5.0D0)
      SS=P(2)
      IF(SS.GT.STS) GO TO 1200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 1300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 1400
      SS=SS+P(5)
      IF(SS.GT.STS) GO TO 1500
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 1800
C  ****  A delta interaction (ICOL=7) may occur when the total
C        interaction probability per unit path length, ST, is
C        larger than STNOW.
      ICOL=7
      DE=0.0D0
      RETURN
C
C  ****  Hard elastic collision (ICOL=2).
C
 1200 ICOL=2
      IF(E.GE.EELMAX(M)) THEN
        TRNDC=RNDCE(M,KE)+(RNDCE(M,KE+1)-RNDCE(M,KE))*XEK
        TA=EXP(AE(M,KE)+(AE(M,KE+1)-AE(M,KE))*XEK)
        TB=BE(M,KE)+(BE(M,KE+1)-BE(M,KE))*XEK
        CALL EELa(TA,TB,TRNDC,RMU)
      ELSE
        TRNDC=RNDCEd(M,KE)+(RNDCEd(M,KE+1)-RNDCEd(M,KE))*XEK
        CALL EELd(TRNDC,RMU)  ! Uses the ELSEPA database.
      ENDIF
      CDT=1.0D0-(RMU+RMU)
      DF=TWOPI*RAND(6.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      DE=0.0D0
      RETURN
C
C  ****  Hard inelastic collision (ICOL=3).
C
 1300 ICOL=3
      DELTA=DEL(M,KE)+(DEL(M,KE+1)-DEL(M,KE))*XEK
      CALL EINa(E,DELTA,DE,EP,CDT,ES,CDTS,M,IOSC)
cx:
      if (de.lt.eabs(1,m)) derest = de
C  ****  Scattering angles (primary electron).
      DF=TWOPI*RAND(7.0D0)
C  ****  Delta ray.
      IF(ES.GT.EABS(1,M)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  New energy and direction.
      IF(EP.GT.EABS(1,M)) THEN
        E=EP
        CALL DIRECT(CDT,DF,U,V,W)
        RETURN
      ENDIF
      DE=E
      E=0.0D0
      RETURN
C
C  ****  Hard bremsstrahlung emission (ICOL=4).
C
 1400 ICOL=4
      CALL EBRa(E,DE,M)
C  ****  Bremsstrahlung photon.
      IF(DE.GT.EABS(2,M)) THEN
        CALL EBRaA(E,DE,CDTS,M)
        DFS=TWOPI*RAND(8.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(DE,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
cx:
      else
        derest = de
      ENDIF
C  ****  New energy.
      E=E-DE
      IF(E.GT.EABS(1,M)) RETURN
      DE=E+DE
      E=0.0D0
      RETURN
C
C  ****  Ionisation of an inner shell (ICOL=5) -does not affect the
C        primary electron.
C
 1500 ICOL=5
      DE=0.0D0
      CALL ESIa(IZA,ISA)
      IF(IZA.LT.1) RETURN
      ILBA(3)=ICOL
      CALL RELAX(IZA,ISA)
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 1800 ICOL=8
      DE=0.0D0
      CALL EAUX
      RETURN
C
C  ************  Photons (KPAR=2).
C
 2000 CONTINUE
C
      STS=ST*RAND(1.0D0)
      SS=P(1)
      IF(SS.GT.STS) GO TO 2100
      SS=SS+P(2)
      IF(SS.GT.STS) GO TO 2200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 2300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 2400
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 2800
C
C  ****  Rayleigh scattering (ICOL=1).
C
 2100 ICOL=1
      DE=0.0D0
      CALL GRAa(E,CDT,M)
      DF=TWOPI*RAND(2.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ****  Compton scattering (ICOL=2).
C
 2200 ICOL=2
      CALL GCOa(E,DE,EP,CDT,ES,CDTS,M,ISHELL)
      IZA=KZCO(M,ISHELL)
      ISA=KSCO(M,ISHELL)
      DF=TWOPI*RAND(3.0D0)
      IF(IZA.GT.0.AND.ISA.LT.10) THEN
        ILBA(3)=ICOL
        CALL RELAX(IZA,ISA)
      ENDIF
C  ****  Compton electron.
      IF(ES.GT.EABS(1,M)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=IZA*1000000+ISA
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  New direction and energy.
      IF(EP.LT.EABS(2,M)) THEN
        DE=E
        E=0.0D0
      ELSE
        CALL DIRECT(CDT,DF,U,V,W)
        E=EP
      ENDIF
      RETURN
C
C  ****  Photoelectric absorption (ICOL=3).
C
 2300 ICOL=3
      CALL GPHa(ES,IZA,ISA)
C  ****  Delta interaction. Introduced to correct for the use of an
C        upper bound of the photoelectric attenuation coefficient.
      IF(IZA.EQ.0) THEN
        ICOL=7
        DE=0.0D0
        RETURN
      ENDIF
C
      IF(ES.GT.EABS(1,M)) THEN
        CALL SAUTER(ES,CDTS)
        DFS=TWOPI*RAND(4.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=IZA*1000000+ISA
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
      IF(ISA.LT.10) THEN
        ILBA(3)=ICOL
        CALL RELAX(IZA,ISA)
      ENDIF
      DE=E
      E=0.0D0
      RETURN
C
C  ****  Electron-positron pair production (ICOL=4).
C
 2400 ICOL=4
      CALL GPPa(EE,CDTE,EP,CDTP)
      DE=E
C  ****  Electron.
      IF(EE.GT.EABS(1,M)) THEN
        DF=TWOPI*RAND(5.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTE,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(EE,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  Positron.
      IF(EP.GT.EABS(3,M)) THEN
        DF=TWOPI*RAND(6.0D0)
        CALL DIRECT(CDTP,DF,U,V,W)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(EP,X,Y,Z,U,V,W,WGHT,3,ILBA)
C  ****  The positron carries a 'latent' energy of 1022 keV.
        DE=DE-TREV
      ELSE
        CALL PANaR
      ENDIF
      E=0.0D0
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 2800 ICOL=8
      DE=0.0D0
      CALL GAUX
      RETURN
C
C  ************  Positrons (KPAR=3).
C
 3000 CONTINUE
      IF(MODE.EQ.1) GO TO 3100
C
C  ****  Artificial soft event (ICOL=1).
C
      ICOL=1
      MODE=1
C
      IF(KSOFTI.EQ.0) THEN
        DE=0.0D0
      ELSE
        EDE0=W1*DST
        VDE0=W2*DST
        FSEDE=MAX(1.0D0-DW1PL(M,KE)*EDE0,0.75D0)
        FSVDE=MAX(1.0D0-DW2PL(M,KE)*EDE0,0.75D0)
        EDE=EDE0*FSEDE
        VDE=VDE0*FSVDE
C  ****  Generation of random values DE with mean EDE and variance VDE.
        SIGMA=SQRT(VDE)
        IF(SIGMA.LT.0.333333333D0*EDE) THEN
C  ****  Truncated Gaussian distribution.
          DE=EDE+RNDG3()*SIGMA
        ELSE
          RU=RAND(1.0D0)
          EDE2=EDE*EDE
          VDE3=3.0D0*VDE
          IF(EDE2.LT.VDE3) THEN
            PNULL=(VDE3-EDE2)/(VDE3+3.0D0*EDE2)
            IF(RU.LT.PNULL) THEN
              DE=0.0D0
            ELSE
              DE=1.5D0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0D0-PNULL)
            ENDIF
          ELSE
            DE=EDE+(2.0D0*RU-1.0D0)*SQRT(VDE3)
          ENDIF
        ENDIF
      ENDIF
C
      E=E-DE
cx:
      derest = de
C  ****  Annihilation at rest.
      IF(E.LT.EABS(3,M)) THEN
        DE=E+DE+TREV
        E=0.0D0
        CALL PANaR
        RETURN
      ENDIF
      IF(KSOFTE.EQ.0) RETURN
C
C  ****  Bielajew's randomly alternate hinge.
C
      IF(RAND(2.0D0).GT.0.5D0.AND.DE.GT.1.0D-3) THEN
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        IF(T1E(M,KE+1).GT.-78.3D0) THEN
          T1=EXP(T1P(M,KE)+(T1P(M,KE+1)-T1P(M,KE))*XEK)
          T2=EXP(T2P(M,KE)+(T2P(M,KE+1)-T2P(M,KE))*XEK)
        ELSE
          T1=0.0D0
          T2=0.0D0
        ENDIF
        IF(T1.LT.1.0D-35) RETURN
      ENDIF
C  ****  1st and 2nd moments of the angular distribution.
      EMU1=0.5D0*(1.0D0-EXP(-DST*T1))
      EMU2=EMU1-(1.0D0-EXP(-DST*T2))/6.0D0
C  ****  Sampling from a two-bar histogram with these moments.
      PNUM=2.0D0*EMU1-3.0D0*EMU2
      PDEN=1.0D0-2.0D0*EMU1
      PB=PNUM/PDEN
      PA=PDEN+PB
      RND=RAND(3.0D0)
      IF(RND.LT.PA) THEN
        CDT=1.0D0-2.0D0*PB*(RND/PA)
      ELSE
        CDT=1.0D0-2.0D0*(PB+(1.0D0-PB)*((RND-PA)/(1.0D0-PA)))
      ENDIF
      DF=TWOPI*RAND(4.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ************  Hard event.
C
 3100 CONTINUE
      MODE=0
C  ****  A delta interaction (ICOL=7) occurs when the maximum
C        allowed step length is exceeded.
      IF(KDELTA.EQ.1) THEN
        ICOL=7
        DE=0.0D0
        RETURN
      ENDIF
C  ****  Interaction probabilities.
      STNOW=P(2)+P(3)+P(4)+P(5)+P(6)
C  ****  Random sampling of the interaction type.
      STS=MAX(STNOW,ST)*RAND(5.0D0)
      SS=P(2)
      IF(SS.GT.STS) GO TO 3200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 3300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 3400
      SS=SS+P(5)
      IF(SS.GT.STS) GO TO 3500
      SS=SS+P(6)
      IF(SS.GT.STS) GO TO 3600
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 3800
C  ****  A delta interaction (ICOL=7) may occur when the total
C        interaction probability per unit path length, ST, is
C        larger than STNOW.
      ICOL=7
      DE=0.0D0
      RETURN
C
C  ****  Hard elastic collision (ICOL=2).
C
 3200 ICOL=2
      IF(E.GE.PELMAX(M)) THEN
        TRNDC=RNDCP(M,KE)+(RNDCP(M,KE+1)-RNDCP(M,KE))*XEK
        TA=EXP(AP(M,KE)+(AP(M,KE+1)-AP(M,KE))*XEK)
        TB=BP(M,KE)+(BP(M,KE+1)-BP(M,KE))*XEK
        CALL EELa(TA,TB,TRNDC,RMU)
      ELSE
        TRNDC=RNDCPd(M,KE)+(RNDCPd(M,KE+1)-RNDCPd(M,KE))*XEK
        CALL PELd(TRNDC,RMU)  ! Uses the ELSEPA database.
      ENDIF
      CDT=1.0D0-2.0D0*RMU
      DF=TWOPI*RAND(6.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      DE=0.0D0
      RETURN
C
C  ****  Hard inelastic collision (ICOL=3).
C
 3300 ICOL=3
      DELTA=DEL(M,KE)+(DEL(M,KE+1)-DEL(M,KE))*XEK
      CALL PINa(E,DELTA,DE,EP,CDT,ES,CDTS,M,IOSC)
cx:
      if (de.lt.eabs(1,m)) derest = de
C  ****  Scattering angles (primary particle).
      DF=TWOPI*RAND(7.0D0)
C  ****  Delta ray.
      IF(ES.GT.EABS(1,M)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  New energy and direction.
      IF(EP.GT.EABS(3,M)) THEN
        E=EP
        CALL DIRECT(CDT,DF,U,V,W)
        RETURN
      ENDIF
      DE=E+TREV
      E=0.0D0
C  ****  Annihilation at rest.
      CALL PANaR
      RETURN
C
C  ****  Hard bremsstrahlung emission (ICOL=4).
C
 3400 ICOL=4
      CALL EBRa(E,DE,M)
C  ****  Bremsstrahlung photon.
      IF(DE.GT.EABS(2,M)) THEN
        CALL EBRAa(E,DE,CDTS,M)
        DFS=TWOPI*RAND(8.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(DE,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
cx:
      else
        derest = de
      ENDIF
C  ****  New energy.
      E=E-DE
      IF(E.GT.EABS(3,M)) RETURN
      DE=E+DE+TREV
      E=0.0D0
C  ****  Annihilation at rest.
      CALL PANaR
      RETURN
C
C  ****  Ionisation of an inner shell (ICOL=5) -does not affect the
C        primary positron.
C
 3500 ICOL=5
      DE=0.0D0
      CALL PSIa(IZA,ISA)
      IF(IZA.LT.1) RETURN
      ILBA(3)=ICOL
      CALL RELAX(IZA,ISA)
      RETURN
C
C  ****  Positron annihilation in flight (ICOL=6).
C
 3600 ICOL=6
      CALL PANa(E1,CDT1,E2,CDT2)
      DF=TWOPI*RAND(9.0D0)
      IF(E1.GT.EABS(2,M)) THEN
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDT1,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(E1,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      ENDIF
      IF(E2.GT.EABS(2,M)) THEN
        DF=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDT2,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(E2,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      ENDIF
      DE=E+TREV
      E=0.0D0
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 3800 ICOL=8
      DE=0.0D0
      CALL PAUX
      RETURN
      END


      SUBROUTINE KNOCKFx(DEF,ICOL)
C*******************************************************************
C*   This routine supersedes the standard KNOCKF of PENELOPE.      *
C*   The main difference between them is that KNOCKX is called     *
C*   here instead of KNOCK. This change is essential for the track *
C*   length tally to work properly.                                *
C*                                                                 *
C*   Notes:                                                        *
C*   -> Parts that have been changed are written in lowercase and  *
C*      marked with 'cx' (except for the this header)              *
C*******************************************************************
C
C  Modified subroutine 'KNOCK' for interaction forcing.
C
C  Output arguments:
C    DEF .... effective energy deposited by the particle in the mater-
C             ial. Includes a weight correction such that no special
C             action needs to be taken in the main program.
C             --> Use with care; when this quantity is considered as
C             the deposited energy, simulated energy deposition spectra
C             are biased.
C    ICOL ... kind of interaction experienced by the particle.
C
C  This subroutine does not modify the weight of the primary particle.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Secondary stack.
      PARAMETER (NMS=2000)
      COMMON/SECST/ES(NMS),XS(NMS),YS(NMS),ZS(NMS),US(NMS),
     1   VS(NMS),WS(NMS),WGHTS(NMS),KS(NMS),IBODYS(NMS),MS(NMS),
     2   ILBS(5,NMS),NSEC
      COMMON/CERSEC/IERSEC
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
      COMMON/CJUMP1/MODE,KSOFTE,KSOFTI,KDELTA
C  ****  Interaction forcing parameters.
      COMMON/CFORCG/P0(8),IFORC(8)
C
      EXTERNAL RAND

cx:   *** This common us used to extract info on the restricted Eloss:
      logical called
      real*8 derest
      common /scofl0/ derest,called
C
C  ****  Store particle state variables before the interaction.
C
      EA=E
      UA=U
      VA=V
      WA=W
      NSECA=NSEC
C
C  ****  Simulate a new interaction.
C

cx:
      call knockx(de,icol)
c-off CALL KNOCK(DE,ICOL)

      IF(IFORC(ICOL).EQ.0) THEN  ! Unforced interaction.
        DEF=DE
        RETURN
      ENDIF
C  ****  Modify the weights of generated secondary particles, if any.
      WFORCE=P0(ICOL)/P(ICOL)
      DEF=DE*WFORCE

cx:
      derest = derest*wforce

      IF(NSEC.GT.NSECA) THEN
C  ****  Allow writing multiple stack overflow warnings.
        IF(IERSEC.NE.0) IERSEC=0
        DO I=NSECA+1,NSEC
          WGHTS(I)=WGHTS(I)*WFORCE
        ENDDO
      ENDIF
C  ****  And set the final state of the primary particle.
      IF(RAND(1.0D0).GT.WFORCE) THEN
C  ****  The primary particle state is not affected.
        E=EA
        U=UA
        V=VA
        W=WA
        IF(KPAR.EQ.2) MODE=1  ! Same energy as in previous calls.
      ELSE
        IF(KPAR.EQ.2) THEN
          IF(ABS(E-EA).LT.1.0D0) THEN
            MODE=1  ! Same energy as in previous calls.
          ELSE
            MODE=0
          ENDIF
        ENDIF
      ENDIF
C
      RETURN
      END


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
