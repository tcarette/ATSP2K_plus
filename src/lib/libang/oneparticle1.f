*     ..........................................................   :
*                                                                  : 
*          Block                                                   : 
*                    O N E   -   T W O                             :
*                                                                  : 
*     For Calculation Angular Momentum Coefficients for            :
*     Relativistic operators                                       :
*                                                                  : 
*     Written by G. Gaigalas,                                      : 
*                  Department  of  Computer Science,               : 
*                  Vanderbilt University,  Nashville               : 
*                                                   October 1996   :  
*                                                                  : 
*     ..........................................................   :
*        i)   one - particle operator                              : 
*     ..........................................................   :
*
*     -------------------------------------------------------------
*      O N E P A R T I C L E 1
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                  N'2 = N2        *
*     Written by G. Gaigalas,                                      *
*     Universite Libre de Bruxelles, Belgium         October 1995  *
*
      SUBROUTINE ONEPARTICLE1(K1,K2,IA,XXX)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      DIMENSION K(2)
      EXTERNAL XXX
      C1=ZERO
      K(1)=K1
      K(2)=K2
      DO 1 I=1,3
        IF(I.EQ.1) THEN
          CALL RLSP0(I,IA,IA,0,IAT)
        ELSE
          J=I-1
          CALL RLSP0(I,IA,IA,2*K(J),IAT)
        ENDIF
        IF(IAT.EQ.0) RETURN
    1 CONTINUE
      DO 2 I=2,3
        J=I-1
        CALL RLSP1(I,IA,2*K(J),0,IAT,REC)
        IF(IAT.EQ.0) RETURN
    2 CONTINUE
      CALL HIBFF(IA,IA,IA,IA,1)
      CALL XXX(IK1(3),ID1(3),INUM,A1)
      IF(DABS(A1).LT.EPS) RETURN
      LA=IJFUL(IA)
      CALL W1(IK1,BK1,ID1,BD1,K(1),K(2),HALF,-HALF,W)
      RECLS=1
      DO 3 I=2,3
        J=I-1
        CALL RLSP1(I,IA,2*K(J),1,IAT,REC)
        RECLS=RECLS*REC
    3 CONTINUE
      C1=A1*W*RECLS/DSQRT(DBLE((2*K(1)+1)*(2*K(2)+1)))
      IF(DABS(C1).GT.EPS)CALL SAVENON(INUM,C1,0,0,LA,0,LA,JA,JB,0)
      RETURN
      END
