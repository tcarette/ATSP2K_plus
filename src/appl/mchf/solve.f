*
*     ------------------------------------------------------------------
*    3-34      S O L V E
*     ------------------------------------------------------------------
*
*       When FIRST is .TRUE., SOLVE computes the potential and exchange
*   function and initializes variables for the i'th radial  equation.
*   The vector P1 is the solution of the radial equation and P2 the
*   variation of the solution with respect to the energy parameter
*   E(I,I).
*
*
      SUBROUTINE SOLVE(I,FIRST)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800,NODD=100)
*
        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        COMMON/LABEL/ EL(NWD),ATOM,TERM
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
        COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
*
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
       COMMON P2(NOD),HQ(NOD),XX(NOD),AC(NODD,NODD),BC(NODD),JV(NODD),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
        LOGICAL FIRST
      DIMENSION ZERO(NOD),P1(NOD)
      EQUIVALENCE (ZERO(1),XX(1)),(PDE(1),P1(1))
      SAVE Zinf,fl,V,b4,CN,C,CD,XY,XP
*
*  *****  IF FIRST IS 'TRUE', CALL POTL AND XCH AND SET UP ARRAYS
*
      IF (.NOT. FIRST) GO TO 17
      CALL POTL(I)
      CALL XCH(I,3)
      ZINF = DMAX1(0.05D0, Z-YR(ND))
      FN = N(I)
      FL = L(I)
      V = YR(1)/R(1)
      B4 = Z*(FL+D4/D3)/((FL+D1)*(FL+D2))
      CN = (D2*Z/FN)**(L(I) +1)
      C = D4*FL +D6
      CD = (FL+D5)**2
      XY = X(1)
      XP = X(2)
      ED = E(I,I)
      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)
      DO 1 J = 3,ND
      X5 = X(J+2)
      X(J) =CH*(-X5+24.D0*(X4+X2) + 194.D0*X3 - X1)/20.D0
      X1 = X2
      X2= X3
      X3 = X4
1     X4 = X5
      X(NO-1) = CH*(X4 + D10*X3 + X2)
      DO 4 J = 1,NO
4     YK(J) = -D2*(Z - YR(J))*R(J) + CD
      X1 =    CH*P(1,I)*(YK(1)+ED*RR(1))
      X2 =    CH*P(2,I)*(YK(2)+ED*RR(2))
      X3 =    CH*P(3,I)*(YK(3)+ED*RR(3))
      X4 =    CH*P(4,I)*(YK(4)+ED*RR(4))
      DO 7 J = 3,ND
      X5 =    CH* P(J+2,I)*(YK(J+2)+ED*RR(J+2))
      X(J) = X(J) - (X5 -D4*(X2 + X4) + D6*X3 +X1)/20.D0
      X1 = X2
      X2 = X3
      X3 = X4
7     X4 = X5
      RL = L(I) + 2.5
      X(2) = R(2)**RL*(X(5)/R(5)**RL - D3*(X(4)/R(4)**RL -
     1     X(3)/R(3)**RL))
*
*  *****  DETERMINE LOWER BOUND ON THE ENERGY PARAMETER
*
      IF (KK .NE. 3) GO TO 80
      DO 11 JJ = 15,ND
      J = NO - JJ
      IF (YK(J) .LT. D0 ) GO TO 63
11    CONTINUE
      WRITE(OUT,12)
12    FORMAT(10X,'POTENTIAL FUNCTION TOO SMALL - 2R*(Z-Y)<(L+.5)**2')
*     STOP
      GO TO 80
63    EM = -YK(J)/RR(J)
      GO TO 81
80    EM = (ZINF/(FN + D5))**2
81    FM = EM
*
*  *****  DETERMINE DIAGONAL ENERGY PARAMETER
*
      F1 = D0
      C11 = D0
      M = MIN0(MAX(I),NO-1)
      DO 5 J = 2,M
      FNUM = P(J+1,I) - P(J,I) - P(J,I) + P(J-1,I)
      FNUM = FNUM - CH*(YK(J+1)*
     1   P(J+1,I) + D10*YK(J)*P(J,I) + YK(J-1)*P(J-1,I))-X(J)
      DEL1 = RR(J+1)*P(J+1,I) + D10*RR(J)*P(J,I) + RR(J-1)*P(J-1,I)
      F1 = F1 +P(J,I)*FNUM
      C11 = C11 + P(J,I)*DEL1
5     CONTINUE
      ED = F1/(C11*CH)
      IF (ED .GT. EM) GO TO 19
*
*  *****  ERROR MESSAGE AND ENERGY ADJUSTMENT FOR AN ENERGY PARAMETER
*  *****  TOO SMALL FOR THE RANGE OF THE FUNCTION
*
      WRITE(OUT,24) ED
24    FORMAT(10X,5HED = ,F10.6,36H; ADJUSTED TO ALLOWED MINIMUM ENERGY )
      ED = EM
      IF ( DABS(FM - E(I,I)) .GT. 1.D-6 .OR. KK .EQ. 3 ) GO TO 19
*
*  ***** RETURN HYDROGENIC FUNCTION
*
      PN = HNORM(N(I),L(I),ZINF)
      DO 65 J = 1,NO
65    PDE(J) = PN*HWF(N(I),L(I),ZINF,R(J))/R2(J)
      AZD = PN*(D2*ZINF/N(I))**(L(I)+1)
      PP = D0
      acc(i) = (acc(i) + 1.d0)/2.d0
      WRITE(OUT,66) EL(I), ZINF
66    FORMAT(//10X, 'RETURN HYDROGENIC FUNCTION FOR ',A3,
     1   ' WITH EFFECTIVE CHARGE ',F10.3)
      RETURN
*
*  *****  CHECK IF UPPER BOUND IS CORRECT
*
19    IF ( D10*ED .LT. EU) GO TO 18
      EU = D10*ED
      FU = EU
18    AZD = AZ(I)
17    DO 26 J=1,NO
      YR(J) = (YK(J) + ED*RR(J))*CH
26    ZERO(J) = D0
*
*  *****  SEARCH FOR THE POINT AT WHICH YR BECOMES POSITIVE
*
      CALL SEARCH(NJ,I)
*
*  *****  COMPUTE STARTING VALUES FROM SERIES EXPANSION
*
      B3 = (V + V + ED - (Z/FN)**2)/C
      DO 6 J = 1,2
      HW  = HWF(N(I),L(I),Z,R(J))/CN
6     HQ(J)   = AZD*(HW + R(J)**(L(I)+3)*B3*(D1-R(J)*B4))/R2(J)
*
*  *****  OBTAIN HOMOGENEOUS SOLUTION
*
      CALL NMRVS(NJ,DELH,MH,HQ,ZERO)
      P1(1) = HQ(1) + XY/C
      P1(2) = HQ(2) + XP/C
*
*  *****  OBTAIN PARTICULAR SOLUTION
*
      CALL NMRVS(NJ,DEL1,M1,P1,X)
*
*  *****  DETERMINE THE ENERGY ADJUSTMENT REQUIRED FOR A SOLUTION WITH
*  *****  GIVEN A0
*
      M = MAX0(M1,MH)
      PNORM = D0
      DO 50 J = 1,M
50    PNORM = PNORM + RR(J)*HQ(J)*P1(J)
      Y1 = P1(NJ-1)
      Y2 = P1(NJ)
      Y3 = P1(NJ+1)
      DELTA = Y2 - Y1 + Y2 - Y3 +YR(NJ-1)*Y1 + D10*YR(NJ)*Y2
     1   + YR(NJ+1)*Y3 + X(NJ)
      DELTAE = HQ(NJ)*DELTA/(H*H*PNORM)
      PP = -DEL1/DELH
*
*  *****  MATCH AT THE JOIN FOR A SOLUTION OF THE DIFFERENTIAL EQUATION
*
      DO 13 J = 1,NO
13    P1(J)   = P1(J) + PP*HQ(J)
*
*  *****  IF  THE EQUATIONS APPEAR TO BE NEARLY
*  ****  SINGULAR, SOLVE THE VARIATIONAL EQUATIONS
*
      IF (KK .NE. 2) RETURN
      X1 = P(1,I)*RR(1)
      X2 = P(2,I)*RR(2)
      P2(1) = X1/C
      P2(2) = X2/C
      DO 8 J = 3,NO
      X3 = P(J,I)*RR(J)
      XX(J-1) = (D10*X2 + X1 + X3)*CH
      X1 = X2
8     X2 = X3
      CALL NMRVS(NJ,DEL2,M2,P2,XX)
      AA = -DEL2/DELH
      M = MAX0(M,M2)
      DO 9 J = 1,NO
9     P2(J) = P2(J) + AA*HQ(J)
      A11 = QUAD(I,M,P2,P2)
      B11 = QUAD(I,M,P1,P2)
      C11 = QUAD(I,M,P1,P1) - D1
      DISC = B11*B11 - A11*C11
      IF ( DISC .LT. D0 ) GO TO 70
      DE1 = -(B11+DSQRT(DISC))/A11
      DE2 = C11/A11/DE1
      IF ( P1(3)+DE1*P2(3) .LT. D0) DE1 = DE2
      GO TO 71
70    DE1 = C11/A11
71    DO 301 J = 1,NO
      P1(J) = P1(J) + DE1*P2(J)
301   CONTINUE
      PP = PP + DE1*AA
      RETURN
      END
