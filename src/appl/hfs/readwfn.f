
*-----------------------------------------------------------------------
*     R E A D W F N
*-----------------------------------------------------------------------
*     This subroutine reads the radial wavefunctions and sorts them
*     according to the order in IAJCMP.
*
      SUBROUTINE READWFN 
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=70)
      INTEGER SMAX(NWD)
      DIMENSION SAZ(NWD),SP(NOD,NWD)
      CHARACTER*3 EL(NWD),SEL(NWD),AT*6,TT*6
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /PARAT/D0,D1,D2,D3,D4,D5,D6,D10,H,H1,NO,ND
      COMMON/ADATA/P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,AZ(NWD),L(NWD),MAX(NWD)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISC(6),ISCW

*     Read radial wavefunctions 

      I=1
5     READ(IREADW,END=999) AT,TT,SEL(I),MR,Z,ETI,EKI,SAZ(I),
     :   (SP(JJ,I),JJ=1,MR)

*     Check if the orbital is one of those required in the present run.

      READ(SEL(I),'(A3)')ITEMP
      DO 10 J=1,MAXORB
         IF (IAJCMP(J).EQ.ITEMP) THEN
            DO 15 JJJ = MR,NO
15          SP(JJJ,I) = D0
            SMAX(I) = MR
Cww            IF (I.GT.(60)) STOP ' Too many electrons: max=(30)'
            IF (I.GT.(NWD)) THEN
               WRITE(IWRITE,*) ' Too many electrons: MAX=',NWD
               STOP
            END IF
            I = I+1
            GO TO 5   
         ENDIF
10    CONTINUE
      GO TO 5 
999   CONTINUE
      ZED=DBLE(Z)

*     Sort the orbitals according to the order in IAJCMP

      DO 20 II=1,MAXORB
         READ(SEL(II),'(A3)')ITEMP
         DO 25 I= 1,MAXORB
            IF (IAJCMP(I).EQ.ITEMP) THEN
               EL(I)=SEL(II)
               AZ(I)=SAZ(II)
               MAX(I)=SMAX(II)
               DO 30 JJ=1,NO
                  P(JJ,I)=SP(JJ,II)
30             CONTINUE
            ENDIF
25       CONTINUE
20    CONTINUE
      DO 35 I=1,MAXORB
         IF (EL(I).EQ.'   ') STOP 'Radial orbital missing.'
35    CONTINUE

      DO 40 I=1,MAXORB
         IF (EL(I)(1:1) .NE. ' ') THEN
            L(I) = LVAL(EL(I)(2:2))
         ELSE
            L(I) = LVAL(EL(I)(3:3))
         END IF
40    CONTINUE
      RHO=-4.D0
      DO 45 J=1,NO
         R(J)=DEXP(RHO)/ZED
         RR(J)=R(J)*R(J)
         R2(J)=DSQRT(R(J))
         RHO=RHO+H
45    CONTINUE
      RETURN
      END
