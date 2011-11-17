
*
*     ------------------------------------------------------------------
*    3-19      M E T H O D
*     ------------------------------------------------------------------
*
*       Uses M1, M2, or M3 to solve the radial equation. If the input
*   data indicated METH(I) = 3, then this  solution  is  returned  to
*   DE.  Otherwise,  the routine searches for an acceptable  solution
*   which  is  both  positive  near  the  origin and has the required
*   number  of nodes.
*
*
      SUBROUTINE METHD1(I)
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
        COMMON/PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
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
     1     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
      LOGICAL V2, FIRST
      DIMENSION P1(220)
      EQUIVALENCE (PDE(1),P1(1))
*
*  *****  'FIRST' MUST BE 'TRUE' THE FIRST TIME SOLVE IS CALLED FOR
*  *****  POTENTIAL AND EXCHANGE TO BE COMPUTED
*  *****  'EU' IS THE UPPER BOUND OF THE ENERGY PARAMETER
*  *****  'EM' IS THE MINIMUM VALUE OF THE ENERGY PARAMETER
*
      FIRST = .TRUE.
      FAIL = .FALSE.
      EM = D0
      EU = ((Z - MIN(D5*S(I),D2*S(I)))/N(I))**2
      FU = EU
      MK = 0
17    CALL SOLVE(I,FIRST)
*
*  *****  IF KK EQUALS 3, OMIT THE NODE CHECKING
*
      IF (KK .EQ. 3) GO TO 51
*
*  *****  COUNT THE NUMBER OF NODES
*
      MN = M
      NC = NODEC(MN)
      IF (TRACE) WRITE(OUT,99) EL(I),NC,MN,NJ,PDE(MN),ED,EU,EM,DELTAE
99    FORMAT(2X,A3,' NC =',I3,' MN =',I3,' NJ =',I3,' PDE(MN) =',
     1   D10.2,' ED =',D10.2,' EU =',D10.2,' EM =',D10.2,
     2   ' DELTAE =',D10.2)
*
*  *****  IF NODE COUNT IS OFF BY NO MORE THAN 1 AND DELTAE IS STILL
*  *****  QUITE LARGE, APPLY THE DELTAE CORRECTION
*
      IF (IABS(NC-NODE) .EQ. 1 .AND. DABS(DELTAE/ED) .GT. 0.02D0)
     1      GO TO 46
*
*  *****  BRANCH ACCORDING TO WHETHER THE NODE COUNT IS TOO SMALL,
*  *****  JUST RIGHT, OR TOO LARGE
*
12    IF (NC - NODE ) 8,9,10
*
*  *****  THE SOLUTION HAS THE CORRECT NUMBER OF NODES
*
9     V2 = DABS(DELTAE)/ED .LT. 1.D-5
      IF (PDE(MN) .LT. D0 .AND. .NOT. V2) GO TO 46
      IF (PDE(MN) .GT. D0) GO TO 51
      DO 52 J = 1,NO
52    PDE(J) = - PDE(J)
      PP = -D2 - PP
51    AZZ = AZD*(D1 + PP)
      CALL EIJSET(I,I,ED)
      RETURN
*
*  *****  THE SOLUTION HAS TOO FEW NODES
*
8     IF (PDE(MN) .LE. D0) GO TO 11
      DEL = D1 - ED/EU
      EU = ED
      IF ( DEL .LT. .05D0) FU = FU*((L(I)+1+NC)/FN)**2.5
       IF (DEL  .GE. .05D0) FU = ED*((L(I)+1+NC)/FN)**2.5
      IF (FU .LT. EM) FU = D5*(EU + EM)
      IF (DABS(FU - ED) .LT. 0.001D0) GO TO 27
      ED = FU
      GO TO 33
*
*  *****  TRY A NEW VALUE OF ED WHICH MUST LIE WITHIN THE UPPER AND
*  *****  LOWER BOUND
*
11    EDP = ED
                    ED = ED*((L(I)+1+NC)/FN)**2.5
      IF (ED .GE. EU ) ED = D5*(EU + EDP)
      IF (ED .LE. EM ) ED = D5*(EM + EDP)
33    MK = MK + 1
      IF ( EU .LE. EM ) WRITE(OUT,30) EM,EU,ED
30    FORMAT(6X,48HWARNING: DIFFICULTY WITH NODE COUNTING PROCEDURE/
     1   6X,42HLOWER BOUND ON ED GREATER THAN UPPER BOUND/
     2   6X,5HEL = ,F10.6,7H  EU = ,F10.6,7H  ED = ,F10.6)
      FIRST = .FALSE.
      IF ( MK .GT. 3*N(I) .OR. EU-EM .LT. FN**(-3)) GO TO 27
      GO TO 17
*
*  *****  THE SOLUTION HAS TOO MANY NODES
*
10    IF (PDE(MN) .LT. D0) GO TO 11
      DEL = D1 - EM/ED
      EM = ED
      IF (DEL .LT. 0.05D0) FM = FM*((L(I)+1+NC)/FN)**2.5
      IF (DEL .GE. 0.05D0) FM = ED*((L(I)+1+NC)/FN)**2.5
      IF (FM .GT. EU) FM = D5*(EU + EM)
      IF (DABS(FM - ED) .LT. 0.001D0) GO TO 27
      ED = FM
       GO TO 33
*
*  *****  ADJUST ENERGY TO LIE BETWEEN UPPER AND LOWER BOUND
*
46    ED = ED - DELTAE
      IF ( ED .GE. EM .AND. ED .LE. EU ) GO TO 33
      EDP = ED
      IF ( NC-NODE .NE. 0 ) ED = (ED+DELTAE)*((L(I)+1+NC)/FN)**2.5
      IF ( ED .GE. EM .AND. ED .LE. EU ) GO TO 33
      ED = EDP + DELTAE + DELTAE
      IF ( ED .GE. EM .AND. ED .LE. EU ) GO TO 33
      ED = ED -DELTAE
      DELTAE = D5*DELTAE
      GO TO 46
*
*  *****  METHOD WAS UNABLE TO FIND AN ACCEPTABLE SOLUTION
*
27    WRITE(OUT,28) KK,EL(I),NC,NJ,ED,EM,EU
28    FORMAT(10X,6HMETHOD,I2,38H UNABLE TO SOLVE EQUATION FOR ELECTRON,
     1   A3/10X,5HNC = ,I3,3X,5HNJ = ,I3,3X,5HED = ,F10.6,3X,5HEL = ,
     2   F10.6,3X,5HEU = ,F10.6)
      FAIL = .TRUE.
      RETURN
      END
