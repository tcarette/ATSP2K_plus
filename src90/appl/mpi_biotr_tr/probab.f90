!
!     -------------------------------------------------------------
!      P R O B A B
!     -------------------------------------------------------------
!                                                                  *
!     CALCULATE THE TRASITION PROBABILITIES                        *
!                                                                  *
!     DA > 0       1->2     ABSORTION                              *
!     DA < 0       1->2     EMISSION                               *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE PROBAB(ICI, IULS, IULSJ, NPAIR, CONFIGI, CONFIGF, IM, IPRINT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE2=>DOUBLE 
      use debug_C
      use dbg_C
      use INOUT_C
      use ems_C
      use ntrm_C
      use state_C
      use mult_C
      use RYDBERG_C
      use PARAM_C
      use LSJ_C
      use consts_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:24:40  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE trp_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ICI 
      INTEGER , INTENT(IN) :: IULS 
      INTEGER , INTENT(IN) :: IULSJ 
      INTEGER , INTENT(IN) :: NPAIR 
      LOGICAL , INTENT(IN) :: IPRINT 
      CHARACTER , INTENT(INOUT) :: CONFIGI*64 
      CHARACTER , INTENT(INOUT) :: CONFIGF*64 
      CHARACTER  :: IM 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, II, I, JVI 
      REAL(DOUBLE2) :: DA_TMP, AJL, AJR, DA, D, DD, ANGS, ANGSA, SIGMA, TRPT, &
         AKI, AKIV, FE, GV, GL, W, PHAS, FLSJ, G, SLSL, AKILSL, GLSL, SLSV, &
         AKILSV, GLSV 
!-----------------------------------------------
!      LOGICAL REL,VOK,IPRINT
!      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
!      COMMON /DBG/IBUGM
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
!      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
!      COMMON /NTRM/NTERMS
!      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
!     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
!      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
!     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
!     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
!      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
!      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
!     :       (QJVL,JVL(1)),(QJVR,JVR(1))
!      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
!      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,LCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON /LSJ/LL1,LL2,IS1,IS2
!      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
!
    2 FORMAT(/,/,/,32X,'NUMBER OF TERMS IN THE ABOVE SUMMATION =',I10) 
    3 FORMAT(/,'------------------------------------------------------'/,&
         'Pair number ',I3,/,/) 
    4 FORMAT(/,' Initial CSF : ',A,' J = ',F3.1) 
    5 FORMAT(' Final   CSF : ',A,' J = ',F3.1,/) 
    6 FORMAT(' 2*j = ',I5,' lbl = ',I5,' total energy = ',F16.7) 
    7 FORMAT(/,10X,'ENERGY DIFFERENCE OF THE STATES:',9X,1P,D15.7,' CM-1',/,51X&
         ,1P,D15.7,' ANGSTROMS'/,51X,1P,D15.7,' A. U.') 
    8 FORMAT(/,' ',' SL =',1P,E15.7,' TRPT =',1P,E15.7,' AKI =',1P,E15.7) 
    9 FORMAT(/,' ',' SV =',1P,E15.7,' TRPT =',1P,E15.7,' AKI =',1P,E15.7) 
   10 FORMAT(' ',' I = ',I3,' GL = ',1P,D15.7,' J = ',I3,' D1 = ',F3.0) 
  237 FORMAT(/,/,10X,'LENGTH   FORMALISM: ',/,10X,'-------- ----------',/) 
  233 FORMAT(/,10X,'SL                                       = ',1P,D15.7) 
  234 FORMAT(10X,'FINAL OSCILLATOR STRENGTH (GF)           = ',1P,D15.7) 
  134 FORMAT(10X,'TRANSITION PROBABILITY IN EMISSION (Aki) =',1P,D16.7,/,/) 
  235 FORMAT(/,/,10X,'VELOCITY FORMALISM: ',/,10X,'-------- ----------',/) 
  236 FORMAT(/,10X,'SV                                       = ',1P,D15.7) 
  238 FORMAT(10X,'FINAL OSCILLATOR STRENGTH (GF)           = ',1P,D15.7) 
  239 FORMAT(10X,'TRANSITION PROBABILITY IN EMISSION (Aki) =',1P,D16.7,/,/) 
   42 FORMAT(10X,' LS RESULTS (length   formalism):  S   = ',D14.6,/,10X,&
         ' ----------                        Aki = ',D14.6,/,10X,&
         '                                   gf  = ',D14.6) 
   43 FORMAT(10X,'            (velocity formalism):  S   = ',D14.6,/,44X,&
         ' Aki = ',D14.6,/,44X,' gf  = ',D14.6) 
   38 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/,1X,A1,I&
         1,2X,'S = ',1P,D12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5) 
   41 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/,1X,A1,I&
         1,2X,'length:   ','S = ',1P,D12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5) 
   40 FORMAT(4X,'velocity: ',' S = ',1P,D12.5,3X,'GF = ',D12.5,3X,'AKI = ',D&
         12.5) 
   39 FORMAT(9X,1P,D12.5,8X,D12.5,9X,D12.5) 
      WRITE (IWRITE, 2) NTERMS 
      WRITE (ISCW, 2) NTERMS 
      DO K = 1, NPAIR 
         DA_TMP = ET2(IR(K)) - ET1(IL(K)) 
         WRITE (6, *) ' npair = ', NPAIR 
!
! --- print heading for the pair (j,j')
!
         WRITE (IWRITE, 3) K 
         IF (ICI /= 0) THEN 
            WRITE (CONFIGI, '(8A8)') (CFG1(II,IL(K)),II=1,8) 
            WRITE (CONFIGF, '(8A8)') (CFG2(II,IR(K)),II=1,8) 
         ENDIF 
         AJL = DBLE(JVL(K))/TWO 
         AJR = DBLE(JVR(K))/TWO 
         WRITE (IWRITE, 4) CONFIGI(1:50), AJL 
         WRITE (IWRITE, 5) CONFIGF(1:50), AJR 
         WRITE (IWRITE, 6) JVL(K), LBL1(IL(K)), ET1(IL(K)) 
         WRITE (IWRITE, 6) JVR(K), LBL2(IR(K)), ET2(IR(K)) 
!
         DA = ET2(IR(K)) - ET1(IL(K)) 
!GG        ryrat = rydbrg/109737.31534
         SL(K) = SL(K)**2*(JVL(K)+1)*(JVR(K)+1) 
         IF (VOK) THEN 
            IF (LAM == 1) SV(K) = (SV(K)/DA)**2*(JVL(K)+1)*(JVR(K)+1) 
            IF (LAM == 2) SV(K) = (TWO*SV(K)/DA)**2*(JVL(K)+1)*(JVR(K)+1) 
         ENDIF 
         IF (IM == 'M') SL(K) = SL(K)*D4*LAM*(2*LAM - 1) 
         IF (IBUGM /= 0) WRITE (6, *) ' k = ', K, ' sl(k) = ', SL(K), &
            ' sv(k) = ', SV(K) 
         IF (DA == D0) CYCLE  
         D = DABS(DA) 
         DD = D*D2*RYDBRG 
         ANGS = D10**8/DD 
         ANGSA = ANGS 
         IF (ANGS > 2000.D0) THEN 
            SIGMA = (1.D8/ANGS)**2 
            ANGSA = ANGS/(D1 + 8342.13D-8 + 206030./(130.D+8 - SIGMA) + 15997./&
               (38.9D+8 - SIGMA)) 
         ENDIF 
         IF (IPRINT) WRITE (IWRITE, 7) DD, ANGS, D 
!
! --- compute the transition probabilities (sec-1) in emission
!     and the gf(1->2) value
!
         I = 2 
         JVI = JVR(K) 
         IF (DA < D0) THEN 
            I = 1 
            JVI = JVL(K) 
!          da = -da
         ENDIF 
         TRPT = TRP(IM,LAM,DD) 
         AKI = SL(K)*TRPT/(JVI + 1) 
         IF (VOK) AKIV = SV(K)*TRPT/(JVI + 1) 
         IF (IBUGM /= 0) WRITE (IWRITE, 8) SL(K), TRPT, AKI 
         IF (IBUGM/=0 .AND. VOK) WRITE (IWRITE, 9) SV(K), TRPT, AKIV 
         FE = ANGS**2*1.499193D-16 
         IF (IBUGM /= 0) WRITE (IWRITE, 10) I, FE, JVI, D1 
         IF (VOK) GV = FE*AKIV*(JVI + 1)*(-D1)**I 
         GL = FE*AKI*(JVI + 1)*(-D1)**I 
         IF (DABS(GL) <= TOL) CYCLE  
!
         WRITE (IWRITE, 237) 
         WRITE (IWRITE, 233) SL(K) 
         WRITE (IWRITE, 234) GL 
         WRITE (IWRITE, 134) AKI 
!
         IF (VOK) THEN 
            WRITE (IWRITE, 235) 
            WRITE (IWRITE, 236) SV(K) 
            WRITE (IWRITE, 238) GV 
            WRITE (IWRITE, 239) AKIV 
         ENDIF 
!
!  *****  OUTPUT ON tr.lsj FILE
!
         IF (REL) THEN 
!
            WRITE (IULSJ, '(/)') 
            IF (DA_TMP > 0) THEN 
               WRITE (IULSJ, '(I4,F14.8,2X,A)') JVL(K), ET1(IL(K)), CONFIGI(1:&
                  50) 
               WRITE (IULSJ, '(I4,F14.8,2X,A)') JVR(K), ET2(IR(K)), CONFIGF(1:&
                  50) 
            ELSE 
               WRITE (IULSJ, '(I4,F14.8,2X,A)') JVR(K), ET2(IR(K)), CONFIGF(1:&
                  50) 
               WRITE (IULSJ, '(I4,F14.8,2X,A)') JVL(K), ET1(IL(K)), CONFIGI(1:&
                  50) 
            ENDIF 
            WRITE (IULSJ, 38) DD, ANGS, ANGSA, IM, LAM, SL(K), ABS(GL), AKI 
            if (VOK) WRITE (IULSJ, 39) SV(K), ABS(GV), AKIV 
 
            CYCLE  
!
! --- calculate the transition properties in LS
!
         ELSE 
!
            IF (DA_TMP > 0) THEN 
               CALL SIXJ (LL2, JVR(K), IS2, JVL(K), LL1, LAM + LAM, 1, W) 
               PHAS = ONE 
               FLSJ = (JVL(K)+1)*(LL2 + 1)*W*W 
               G = (LL2 + 1)*(IS2 + 1) 
            ELSE 
               CALL SIXJ (LL1, JVL(K), IS1, JVR(K), LL2, LAM + LAM, 1, W) 
               PHAS = -ONE 
               FLSJ = (JVR(K)+1)*(LL1 + 1)*W*W 
               G = (LL1 + 1)*(IS1 + 1) 
            ENDIF 
            SLSL = SL(K)*(IS1 + 1)/((JVR(K)+1)*(JVL(K)+1)*W**2) 
            AKILSL = AKI/FLSJ 
            GLSL = FE*AKILSL*G*PHAS 
            WRITE (IWRITE, 42) SLSL, AKILSL, GLSL 
!
            IF (VOK) THEN 
               SLSV = SV(K)*(IS1 + 1)/((JVR(K)+1)*(JVL(K)+1)*W*W) 
               AKILSV = AKIV/FLSJ 
               GLSV = FE*AKILSV*G*PHAS 
               WRITE (IWRITE, 43) SLSV, AKILSV, GLSV 
            ENDIF 
!
!  *****  OUTPUT ON tr.ls FILE
!
            WRITE (IULS, '(/)') 
            IF (DA_TMP > 0) THEN 
               WRITE (IULS, '(I4,F14.8,2X,A)') JVL(K), ET1(IL(K)), CONFIGI(1:50&
                  ) 
               WRITE (IULS, '(I4,F14.8,2X,A)') JVR(K), ET2(IR(K)), CONFIGF(1:50&
                  ) 
            ELSE 
               WRITE (IULS, '(I4,F14.8,2X,A)') JVR(K), ET2(IR(K)), CONFIGF(1:50&
                  ) 
               WRITE (IULS, '(I4,F14.8,2X,A)') JVL(K), ET1(IL(K)), CONFIGI(1:50&
                  ) 
            ENDIF 
            WRITE (IULS, 41) DD, ANGS, ANGSA, IM, LAM, SLSL, ABS(GLSL), AKILSL 
            IF (VOK) WRITE (IULS, 40) SLSV, ABS(GLSV), AKILSV 
         ENDIF 
 
      END DO 
      RETURN  
      END SUBROUTINE PROBAB 
