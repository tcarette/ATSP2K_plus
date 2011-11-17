C
C     ------------------------------------------------------------------
C    3-24      P O T L
C     ------------------------------------------------------------------
C
C       Computes and stores the potential function
C                              2(k-1)
C              YR = SUM  a    Y      (j,j;r)
C                   j,k   ijk
C
        SUBROUTINE POTL(I)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
      POINTER (pkval,kval(1)),(pcoef,coef(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      DO 1 J=1,NO
1       YR(J) = D0
        DO 2 J = 1, NWF
           IF (I.GT.NCLOSD .AND. J.GT.NCLOSD) GO TO 2
           C = SUM(J)
           IF ( I.EQ.J ) C = C - D1
           CALL YKF(J,J,0,REL)
           DO 3 JJ = 1,NO
              YR(JJ) = YR(JJ) + C*YK(JJ)
3          CONTINUE
           IF ( I.EQ.J .AND. L(I) .GT. 0) THEN
              DO 4 K = 2,2*L(I),2
                 CC = -C*CA(L(I),K)
                 CALL YKF(I,I,K,REL)
                 DO 5 JJ = 1,NO
                    YR(JJ) = YR(JJ) + CC*YK(JJ)
5                CONTINUE
4             CONTINUE
           END IF
2       CONTINUE
*
        SUMI = SUM(I)
        IBEGIN = 1
        IEND = INTPTR(1)
        DO 10 J = IBEGIN,IEND
           IE = 0
           call unpacki(1,j,kv,iel1,iel2,iel3,iel4)
           IF (IEL1 .EQ. I) THEN
              IE = IEL2
           ELSE IF (IEL2 .EQ. I) THEN
              IE = IEL1
           END IF
           IF (IE .NE. 0) THEN
              C = COEF(J)/SUMI
              IF (IEL1 .EQ. IEL2) C = 2*C
              CALL YKF(IE,IE,KV,REL)
* 	      print '(10f8.5)', p(1:nod,1)
              DO 12 JJ = 1,NO
                 YR(JJ) = YR(JJ) + C*YK(JJ)
 12           CONTINUE
*	      print *, 'Coef, SUmi', coef(j), sumi
*	    print '(10f8.5)', yr
           END IF
 10     CONTINUE
        END
