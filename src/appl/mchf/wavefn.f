*
*     ------------------------------------------------------------------
*    3-36      W A V E F N
*     ------------------------------------------------------------------
*
*       This routine initializes radial functions by the procedure
*   indicated by IND(I).
*
*         Value of IND(I)     Method
*         ---------------     ------
*             -1           Functions read from unit IU2
*              0           Screened hydrogenic functions with ZZ=Z-S(I)
*              1           Functions in memory left unchanged
*                                                  0
*   The set of functions are then orthogonalized, Y (i, i;r) and the
*   diagonal energy parameters computed, when necessary.
*
*
      SUBROUTINE WAVEFN
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800, MTERM=20,MEIG=20)
*
        CHARACTER EL*3,ATOM*6,ATM*6,TRM*6,TERM*6
        CHARACTER EL1*3,AT*6,TT*6,TITLE*24
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
      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
*
      COMMON ZZ(NWD),IND(NWD),PN,Z2,FN,M,K,ZT,
     :   ETI,EKI,AZI,PT(NOD),MT,ATM(NWD),TRM(NWD)
*
*
*  *****  GENERATE ARRAYS FOR R,R*R AND DSQRT(R) WITH A CONSTANT MESH
*  *****  SIZE IN THE LOG(Z*R) VARIABLE
*
      DO 1 I=1,NO
      R(I)= EXP(RHO)/Z
      RR(I) = R(I)*R(I)
      R2(I) = DSQRT(R(I))
1     RHO = RHO + H
      RHO = RHO - NO*H
*
*  ***** READ THE WAVEFUNCTIONS
*
      IF (IUF .EQ. 0) GO TO 5
2     READ(IUF,END=5) AT,TT,EL1,MM,ZT,ETI,EKI,AZI,(PT(J),J=1,MM)
      M = MM
      CALL EPTR(EL,EL1,I,*2)
      IF ( I .GT. 0 .AND. IND(I) .EQ. -1) THEN
         ATM(I) = AT
         TRM(I) = TT
         MAX(I) = M
         ZZ(I)  = ZT
         C = D1
         IF ( Z .NE. ZT ) C = Z/ZT
*
*  *****  SCALE RESULTS IF CARDS ARE FOR AN ATOM WITH A DIFFERENT Z
*
         CALL EIJSET(I,I,C*C*ETI)
         AZ(I)  = AZI*C**(L(I)+1)*DSQRT(C)
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
*
*  *****  SET REMAINING VALUES IN THE RANGE = 0.
*
         IF ( M .EQ. NO ) GO TO 12
         M = M +1
         DO 13  J=M,NO
13       P(J,I) = D0
12       IND(I) = -2
      ENDIF
      GO TO 2
*
*  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
*
5     DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
*
*  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
*
7     IF ( IND(I) .EQ. -2 ) GO TO 9
      IND(I) = 0
      WRITE(OUT,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
*
*  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
*  *****  HYDROGENIC APPROXIMATION
*
8     PN = HNORM(N(I),L(I),Z-S(I))
      DO 3 J=1,NO
      P(J,I) = PN*HWF(N(I),L(I),Z-S(I),R(J))/R2(J)
3     CONTINUE
      M = NO
30    IF ( DABS(P(M,I)) .GT. 1.D-15 ) GO TO 31
      P(M,I) = D0
      M = M-1
      GO TO 30
31    MAX(I) = M
*
*  ***** SET THE AZ(I) VALUE
*
      AZ(I) = PN*(D2*(Z - D5*S(I))/N(I))**(L(I) + 1)
      CALL EIJSET(I,I,D0)
9     continue
*
*  *****  ORTHOGONALIZE large n to smaller n
*
      DO 10 i= 1,nwf    
      DO 6 II =1,nwf
        if (e(i,ii) .ne. D0) then
*         .. we have an orthogonality condition
          if (n(i) .gt. n(ii) ) then
*           .. orthogonalize larger n to smaller
            PN = QUADR(I,II,0)
            IF ( DABS(PN) .GT. 1.D-8 ) THEN
              PNN = DSQRT(D1 - PN*PN)
              IF (P(50,I) - PN*P(50,II) .LT. D0) PNN = -PNN
              M = MAX0(MAX(I),MAX(II))
              DO 25 J = 1,M
25              P(J,I) =(P(J,I) - PN*P(J,II))/PNN
            END IF
          end if
        end if
6     CONTINUE
10    CONTINUE
      WRITE(PRI,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     1   4X,5HSIGMA,6X,5HE(NL),4X,6HAZ(NL),4X,9HFUNCTIONS//)
*
*  *****  COMPUTE ONE-ELECTRON ENERGY PARAMETERS IF THEY WERE NOT
*  *****  SPECIFIED ON INPUT.
*
      DO 15 I = 1,NWF
*     IF (E(I,I) .EQ. D0) E(I,I) = HL(EL,I,I,REL) - EKIN(I,I)
      K = IND(I) + 2
      IF ( IND(I) .EQ. -2 ) THEN
           TITLE = ' SCALED '//ATM(I)//TRM(I)
        ELSE IF (IND(I) .EQ. 0) THEN
           TITLE = ' SCREENED HYDROGENIC'
        ELSE
           TITLE = ' UNCHANGED'
      END IF
17    WRITE(PRI,19) EL(I),S(I),E(I,I),AZ(I),TITLE
19    FORMAT(9X,A3,F9.2,F11.3,F10.3,3X,A24)
15    CONTINUE
      IF ( IUF .NE. 0) REWIND(UNIT=IUF)
      RETURN
      END
