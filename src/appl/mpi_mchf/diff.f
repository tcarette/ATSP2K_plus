C
C     ------------------------------------------------------------------
C    3-6       D I F F
C     ------------------------------------------------------------------
C
C
C       Stores LP  in the array YK.  The difference approximation of
C                i
C   Eq. (6-14) is used.
C
C
      SUBROUTINE DIFF(I)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=94,NOD=220,NOFFD=800)
*
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
        include 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
Ctc
*
        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        COMMON/LABEL/ EL(NWD),ATOM,TERM
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
        COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
        COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
*
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :      QMETH,QIEPTR
*
      POINTER(pkval,kval(1)),(pvalue,value(1)),
     :       (pih,ih(1)),(pjh,jh(1)),
     :       (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :       (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
C
C  *****  FORM DD + 2Z/R -L(L+1)/RR|P(I)>
C
      MM = MAX(I) - 3
      FL = L(I)
      TWOZ = Z + Z
      C = (FL+D5)**2
      HH = 180.D0*H*H
      DO 11 K =  4,MM
11    YK(K) = (D2*(P(K+3,I)+P(K-3,I)) - 27.D0*(P(K+2,I)+P(K-2,I)) +
     1   270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I))/HH +
     2   P(K,I)*(TWOZ*R(K) - C)
C
C  *****  BECAUSE OF THE POSSIBILITY OF EXTENSIVE CANCELLATION NEAR THE
C  *****  ORIGIN, SEARCH FOR THE POINT WHERE THE ASYMPTOTIC BEHAVIOUR
C  *****  BEGINS AND SMOOTH THE ORIGIN.
C
      LEXP = L(I) + 2
      Y1 = YK(4)/R2(4)/R(4)**LEXP
      Y2 = YK(5)/R2(5)/R(5)**LEXP
      DO 1 K = 4,100
      KP = K+2
      Y3 = YK(KP)/R2(KP)/R(KP)**LEXP
      IF (Y2 .EQ. D0) GO TO 1
      IF (DABS(Y1/Y2 - D1) .LT..2D0 .OR. DABS(Y3/Y2 - D1) .LT..2D0)
     1       GO TO 2
      Y1 = Y2
      Y2 = Y3
1     CONTINUE
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc    WRITE(OUT,3)  I
       WRITE(80+myid,3)  I
Ctc
3     FORMAT(6X, 'ASYMPTOTIC REGION NOT FOUND FOR FUNCTION NUMBER',I3)
      STOP
C
C  *****  ASYMPTOTIC REGION HAS BEEN FOUND
C
2     KP = K
      KM = KP - 1
      DO 4 K = 1,KM
4     YK(K) = Y1*R2(K)*R(K)**LEXP
      MM = MM + 1
      YK(MM) = (-(P(MM+2,I)+P(MM-2,I)) + D16*(P(MM+1,I)+P(MM-1,I))
     1      -D30*P(MM,I))/(D12*H*H) + P(MM,I)*(TWOZ*R(MM) - C)
      MM = MM + 1
      YK(MM) = (P(MM+1,I) + P(MM-1,I) - D2*P(MM,I))/(H*H) +
     1   P(MM,I)*(TWOZ*R(MM) - C)
      MM = MM+1
      DO 5 K =MM,NO
5     YK(K) = D0
      RETURN
      END
