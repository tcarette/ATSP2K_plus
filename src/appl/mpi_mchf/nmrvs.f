*
*
*     ------------------------------------------------------------------
*    3-20      N M R V S
*     ------------------------------------------------------------------
*
*       Given two starting values, PDE(1) and PDE(2), values of PDE(j),
*   j=3,4,...,NJ+1 are obtained by outward integration of
*               Y" = YR y + F
*   using the discretization  of  Eq.  (6-27 )  with  the  difference
*   correction.  With PDE(NJ) given, the tail procedure is applied to
*   PDE(j),j=NJ+1,  NJ+2,...,MM, where MM is determined automatically
*   and DELTA is the difference between  PDE(NJ+1)  for  outward  and
*   inward integration. (See Eq 6-32, 6-33, and 6-37 for further
*   details.)
*
*
      SUBROUTINE NMRVS(NJ,DELTA,MM,PP,F)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
        include 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
Ctc
*
        PARAMETER (NWD=94,NOD=220,NOFFD=800,NTAIL=150)
*
        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        COMMON/LABEL/ EL(NWD),ATOM,TERM
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :  NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
        COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
      DIMENSION PP(NOD),F(NOD),A(NTAIL),D(NTAIL)
      EQUIVALENCE (G,G3)
*
*  *****  INTEGRATE OUTWARD TO NJ+1
*
      Y1 = PP(1)
      Y2= PP(2)
      G1 = YR(1)
      G2 = YR(2)
      M = NJ + 1
      DO 1 I = 3,M
      G3 = YR(I)
      Y3 = (Y2+Y2-Y1 + (D10*G2*Y2 + G1*Y1) + F(I-1)) / (D1 - G3)
      PP(I) = Y3
      Y1 = Y2
      Y2 = Y3
      G1 = G2
1     G2 = G3
      DELTA = Y3
*
*  *****  APPLY THE TAIL PROCEDURE
*
      K = 1
      PP(M) = -(D1 - G1)*Y1 + F(M)
      A(1) = D1 - G
      D(1) = -(D2 + D10*G)
*
22    RATIO = A(K)/D(K)
      IF (K .GE. (NTAIL)-1 .OR. M .EQ. ND) GO TO 23
      K = K +1
      M = M+1
      G = YR(M)
      A(K) = D1 - G
      D(K) = -(D2 + D10*G) - A(K)*RATIO
      PP(M) = -PP(M-1)*RATIO + F(M)
      IF ((DABS(PP(M))+DABS(PP(M-1))/DABS(D(K))) .GT. TOL
     :  .OR. K .LT. 9) GO TO 22
*
20    CONTINUE
*     CON =DSQRT(EH)*EXP(-DSQRT(DABS(G/CH-.25)/RR(M))*(R(M+1)-R(M)))
*     PP(M) = PP(M)/(D(K) + CON*(D1-  YR(M+1)))
      PP(M) = PP(M)/D(K)
      J = M+1
      DO 2 I= J,NO
2     PP(I) = D0
      DO 3 J = 2,K
      I = M-J+1
      II = K-J+1
3     PP(I) = (PP(I)-A(II+1)*PP(I+1))/D(II)
*
*  *****  SET DELTA = DIFFERENCE OF THE TWO SOLUTIONS AT NJ+1
*  *****         MM = NUMBER OF POINTS IN THE RANGE OF THE SOLUTION
*
      DELTA = DELTA - PP(I)
      MM = M
      RETURN
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc 23    WRITE(OUT,24)
23    WRITE(80+myid,24)
Ctc
24    FORMAT(6X,52HWARNING: FUNCTIONS TRUNCATED BY NMRVS IN TAIL REGION)
      GO TO 20
      END
