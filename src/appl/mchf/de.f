*
*     ------------------------------------------------------------------
*    3-4       D E
*     ------------------------------------------------------------------
*
*       This routine controls the solution of the differenttial equation
*   for the radial function P  .  One of three methods is selected -
*                            I1
*   M1, M2, or M3 -  for solving the equations,  the  initial  choice
*   being determined by an input paramter, METH(I), except when no
*   exchange is present, in which case M2 is selected. (For further
*   information see Sec. 7-4)
*
*        Value of METH(I)     Method
*        ---------------      ------
*        < or =1            M1 with search for an acceptable solution
*             =2            M2 with search for an acceptable solution
*             =3            M3 without any checking
*
*   If M1 fails to find an acceptable solution, the radial  functions
*   are  orthogonalized,  off-diagonal  energy parameters recomputed,
*   and the method tried again.   Should it continue to fail, METH(I)
*   is set to 2.
*
*
      SUBROUTINE DE(I1,ivar)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800,NODD=100)
        INTEGER ivar(NWD)
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
      POINTER(pkval,kval(1)),(pvalue,value(1)),
     :       (pih,ih(1)),(pjh,jh(1)),
     :       (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :       (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
       COMMON P2(NOD),HQ(NOD),XX(NOD),AC(NODD,NODD),BC(NODD),JV(NODD),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
        LOGICAL DIAG
        CHARACTER*2 ASTER(3)
        DATA ASTER/'  ','* ','**'/
        ed1=0
*
      I = I1
      ED2 = E(I,I)
      KK= MAX0(1,METH(I))
      IF (NWF .EQ. 1) KK = 2
      NODE = N(I) - L(I) - 1
*
*  *****  CALL METHD1 TO SOLVE THE DIFFERENTIAL EQUATION
*
      CALL METHD1(I)
      IF ( FAIL ) GO TO 25
*
12    PN = DSQRT(QUAD(I,M,PDE,PDE))
      DO 9 J = 1,M
9     PDE(J) = PDE(J)/PN
      AZZ = AZZ/PN
*
*  *****  CHECK IF METHOD 2 SHOULD BE USED
*
      IF ( KK .NE. 1 ) GO TO 13
      IF (DABS(D1 -ED2/E(I,I)) .GT. 0.01D0  .OR.
     1    DMAX1(DABS(D1 - PN), DABS(D1/PN - D1)) .LT. 0.10D0 )
     2    GO TO 13
       METH(I) = 2
       KK = 2
       GO TO 25
*
*  *****  SET THE ACCELERATING PARAMETER
*
*
13    IF (IPR .NE. I ) GO TO 33
      ED2 = ED2 - E(I,I)
      IF (ED1*ED2 .GT. D0) ACC(I) = .75*ACC(I)
      IF (ED1*ED2 .LT. D0) ACC(I) = (D1 + D3*ACC(I))/D4
33    C = ACC(I)
      CD = D1 - C
*
*   *****  IMPROVE THE ESTIMATES
*
      MAX(I) = M
      DP     = D0
      DO 21 J = 1,M
      DIFF = P(J,I)-PDE(J)
      DP     = DMAX1(DP    ,DABS(DIFF)*R2(J))
21     P(J,I) = PDE(J) + C*DIFF
      IF (M .EQ. NO) GO TO 26
      M = M + 1
      DO 24 J = M,NO
24     P(J,I) = D0
      AZ(I) = CD*AZZ + C*AZ(I)
      MAX(i) = m
      IF (C .NE. 0.d0) then
	PNN = 1.d0/sqrt(quadr(i,i,0))
	DO j = 1,m
	  P(j,i) = PNN*P(j,i)
	END DO
	AZ(i) = PNN*AZ(I)
      END IF
      AZZ = AZ(I)
*
*  *****  CHECK THE ORTHOGONALIZATION
*
26    CONTINUE
      if (omit) go to 68
      nit = nwf-ib+1
      DPNEW = DP/DSQRT(SUM(I))
      IBEGIN = 1
      IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
      IP = IBEGIN
      IJ = 0
50    JI = IJE(IP)
      IF (JI .NE. I ) THEN
        jjv = 0
        do 501 jj = 1,nit
          if(ivar(jj).eq.ji) jjv=jj
501     continue
      IF (jjv.ne.0 .AND. DPM(JI) .GT. DPNEW*DSQRT(SUM(JI))) THEN
*
*               The JI orbital should be orthogonalized
*
         C = QUADR(I,JI,0)
         MM = MAX0(MAX(JI),MAX(I))
         DO 51 J = 1,MM
            P(J,JI) = P(J,JI) - C*P(J,I)
51       CONTINUE
         C2 = DSQRT(QUADR(JI,JI,0))
         DO 52 J = 1,MM
            P(J,JI) = P(J,JI)/C2
52       CONTINUE
         VARIED(JI) = .TRUE.
         MAX(JI) = MM
         AZ(JI) = (AZ(JI) - C*AZ(I))/C2
         WRITE(OUT,63) EL(I),EL(JI),C
      ELSE
*
*              The I'th orbital must be orthogonalized
*
         IJ = IJ + 1
         IF (IJ .GT. (NODD))then
           write(iscw,*) ' TOO MANY ORTHOGONALITY CONDITIONS'
           stop
         end if
         JV(IJ) = JI
      END IF
      END IF
      IP = IP + 1
      IF (IP .LE. IEPTR(I)) GO TO 50
      IF (IJ .NE. 0 ) THEN
         DIAG = .TRUE.
         DO 61 J = 1,IJ
            BC(J) = QUADR(I,JV(J),0)
            AC(J,J) = D1
            DO 62 JJ = J+1, IJ
               IF (E(JV(J),JV(JJ)) .NE. D0 ) THEN
                  AC(J,JJ) = D0
                  AC(JJ,J) = D0
                ELSE
                  AC(J,JJ) = QUADR(JV(J),JV(JJ),0)
                  AC(JJ,J) = AC(J,JJ)
                  DIAG = .FALSE.
               END IF
62          CONTINUE
61       CONTINUE
         IF ( .NOT. DIAG .AND. IJ .GT. 1)
     :              CALL LINEQN(NODD,IJ,AC,BC)
         M = MAX(I)
         DO 65 JJ = 1,IJ
            C = BC(JJ)
            WRITE(OUT,63) EL(JV(JJ)),EL(I),C
63          FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
            M = MAX0(M,MAX(JV(JJ)))
            DO 64 J = 1,M
               P(J,I) = P(J,I) - C*P(J,JV(JJ))
64          CONTINUE
            AZZ = AZZ - C*AZ(JV(JJ))
65       CONTINUE
         PNN = DSQRT(QUADR(I,I,0))
         DO 66 J = 1,M
            P(J,I) = P(J,I)/PNN
66       CONTINUE
         AZZ = AZZ/PNN
      END IF
      M = NO
67    IF (DABS(P(M,I)) .LT. 1.D-15) THEN
         P(M,I) = D0
         M = M-1
         GO TO 67
      END IF
68    MAX(I) = M
      IF (AZZ .GT. D0) AZ(I) = DMAX1(AZZ,D5*AZ(I))
      WRITE(OUT,17) EL(I),E(I,I),AZ(I),PN,ASTER(KK),DP
17    FORMAT(20X,A3,2F15.7,F12.7, A2,1PD10.2)
      DPM(I) = DP
      IF (IPR .EQ. I1) ED1 = ED2
      IF (IPR .NE. I1) ED1 = ED2 - E(I1,I1)
      IPR = I1
      VARIED(I) = .TRUE.
      RETURN
*
*  *****  IF METHD1 FAILED TO FIND AN ACCEPTABLE SOLUTION, ORTHOGONALIZE
*  *****  THE ESTIMATES AND TRY AGAIN
*
25    IF (I .EQ. IB) GO TO 27
      CALL ORTHOG(ivar)
      CALL GRANGE(ivar)
27    CALL METHD1(I)
      IF ( FAIL ) GO TO 23
      GO TO 12
*
*  *****  ERROR RETURN FROM SECOND TRY.  IF M1 WAS USED,SWITCH TO
*         M2 AND TRY ONCE MORE.
*
23    IF ( KK .EQ. 2) RETURN
      KK = 2
      GO TO 27
      END
