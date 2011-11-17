*
*     ------------------------------------------------------------------
*    3-23      O U T P U T
*     ------------------------------------------------------------------
*
*       The radial functions and orthogonality integrals are printed,
*   if PRINT is .TRUE.   The  functions  will  also  be  punched  (or
*   stored) on unit OUF, if OUF .NE. 0.
*
*
      SUBROUTINE OUTPUT(PRINT)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=60,NOD=220,NOFFD=800)
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
*
        LOGICAL PRINT
      DIMENSION POUT(8)
      IF ( .NOT. PRINT ) GO TO 31
*
*  *****  PRINT RADIAL FUNCTIONS, 7 PER PAGE
*
      ML = 1
2     MU = MIN0(ML+7,NWF)
      I = MU - ML + 1
      MX = 0
      DO 1 J = ML,MU
1     MX = MAX0(MX,MAX(J))
      WRITE(PRI,5) ATOM,TERM,(EL(J),J=ML,MU)
5     FORMAT(1H1,9X,19HWAVE FUNCTIONS FOR  ,2A6//10X,1HR,8(10X,A3))
      K= 0
      KK = 0
      DO 6 J = 1,MX
      DO 9 JJ = ML,MU
      IJ = JJ - ML + 1
9     POUT(IJ) = P(J,JJ)*R2(J)
      K = K+1
      IF (K .LE. 10) GO TO 6
      K = 1
      KK = KK+1
      IF (KK .LT. 5) GO TO 21
      KK = 0
      WRITE(PRI,23)
23    FORMAT(1H1//)
      GO TO 6
21    WRITE(PRI,8)
8     FORMAT(1X)
6     WRITE(PRI,10) R(J),(POUT(JJ),JJ=1,I)
10    FORMAT(F13.5,F15.6,7F13.6)
      DO 15 J = ML,MU
      IJ = J - ML + 1
15    POUT(IJ) = DPM(J)
      WRITE(PRI,16) (POUT(J),J=1,I)
16    FORMAT(4X,10HMAX. DIFF. ,F15.7,7F13.7)
      ML = ML+8
      IF (ML .LE. NWF) GO TO 2
31    IF ( NWF .LE. 1) GO TO 30
*
*  *****  PRINT ORTHOGONALITY INTEGRALS
*
*      WRITE(PRI,11) ATOM,TERM
11    FORMAT(////10X,33HORTHOGONALITY INTEGRALS FOR ATOM ,A6,6H TERM ,A6
     1   //20X, 4H(NL),3X,4H(NL),7X,8HINTEGRAL //)
      LM = IB
      ML = MAX0(2,LM)
      DO 12 I = ML,NWF
      JF = I - 1
      DO 13 J = 1,JF
      IF (L(I) .NE. L(J)) GO TO 13
      T = QUADR(I,J,0)
      IF (DABS(T) .GT. 1.D-8)  WRITE(PRI,17) EL(I),EL(J),T
17     FORMAT(21X,A3,4X,A3,F15.8)
13    CONTINUE
12    CONTINUE
30    IF ( OUF .EQ. 0) GO TO 14
*
*  *****  OUTPUT FUNCTIONS ON UNIT OUF FOR FUTURE INPUT
*
*         EKI retained only for compatibility with MCHF format
*
      open (unit=ouf, file='wfn.out', form = 'unformatted') 
      DO 3 I = 1,NWF
      MMX = MAX(I)
	EKI=0
      WRITE (OUF) ATOM,TERM,EL(I),MMX,Z,E(I,I),EKI,AZ(I),
     :            (P(J,I),J=1,MMX)
3     CONTINUE
      close (ouf)
*
14    RETURN
      END
