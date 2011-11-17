*
*     ------------------------------------------------------------------
*    3-31      S C A L E
*     ------------------------------------------------------------------
*
*       The current radial functions are scaled according to the
*   procedures of Sec. 7-2  .   Values of AZ and E(I,I), the starting
*   values and the diagonal energy parameters are also scaled.
*
*
      SUBROUTINE SCALE(ZZ)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800)
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
      COMMON RATIO,SR,SC,SS,F0,F1,F2,F3,PNORM,THETA,K
      DIMENSION RS(NOD),PS(NOD)
      EQUIVALENCE (RS(1),YR(1)),(PS(1),X(1))
*
*  *****  SCALE VALUES OF R=RS, P=PS AND ONE-ELECTRON PARAMETERS.
*  *****  GENERATE NEW VALUES OF R, R*R, AND DSQRT(R)
*
      RATIO = Z/ZZ
      SR = DSQRT(RATIO)
      DO 1 J = 1,NO
      R(J) = R(J)*RATIO
      RR(J) = R(J)*R(J)
1     R2(J) = R2(J)*SR
      DO 2 I = 1,NWF
      SC = (ZZ-S(I))/(Z-S(I))
      SS = SC*RATIO
      ED = E(I,I)*SC**2
      CALL EIJSET(I,I,ED)
      DO 3 J = 1,NO
      RS(J) = R(J)/SS
3     PS(J) = P(J,I)*SC
      SC = (ZZ - D5*S(I))/(Z - D5*S(I))
      AZ(I) = AZ(I)*SC**(L(I)+1)*DSQRT(SC)
      K = 3
*
*  *****  INTERPOLATE THE (RS,PS) FUNCTIONS FOR VALUES OF P AT THE SET
*  *****  OF POINTS R
*
      DO 4 J = 1,NO
*
*  *****  SEARCH FOR THE NEAREST ENTRIES IN THE (RS,PS) TABLE
*
5     IF (K .EQ. ND) GO TO 7
      IF (RS(K) .GT. R(J)) GO TO 6
      K = K + 1
      GO TO 5
*
*  *****  INTERPOLATE
*
6     THETA = DLOG(R(J)/RS(K-1))/H
      F0 = PS(K-2)
      F1 = PS(K-1)
      F2 = PS(K)
      F3 = PS(K+1)
      P(J,I) = D5*(F1+F2) + (THETA -D5)*(F2 - F1) +
     1   THETA*(THETA - D1)*(F0 - F1 - F2 + F3)/D4
      GO TO 4
7     P(J,I) = D0
4     CONTINUE
      MAX(I) = NO
*
*  *****NORMALIZE THE INTERPOLATED FUNCTION
*
      PNORM = DSQRT(QUADR(I,I,0))
      DO 10 J = 1,NO
10    P(J,I) = P(J,I)/PNORM
2     CONTINUE
      Z = ZZ
      RETURN
      END
