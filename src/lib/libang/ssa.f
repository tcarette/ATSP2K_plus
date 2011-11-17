*
*     -------------------------------------------------------------
*      S S A
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF   SPIN OTHER ORBIT  INTERACTIONS BETWEEN THE ELECTRONS    *
*                                                                  *
*     (n l L S  n l L S ::                 ::n l L S  n l L S )    *
*       1 1 1 1  2 2 2 2                      3 3 3 3  4 4 4 4     *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
*
      SUBROUTINE SSA(L1,L2,L3,L4,KL1,KL2,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      AA=ZERO
      IF(ITTK(L1,L3,KL1).EQ.0)RETURN
      IF(ITTK(L2,L4,KL2).EQ.0)RETURN
      AA=RME (L1,L3,KL1)
      AA=AA*RME (L2,L4,KL2)
      IS=4*(2*KL2+5)*(KL2+2)*(2*KL2+3)*(KL2+1)*(2*KL2+1)
      AA=AA*DSQRT(DBLE(IS))*TWO*THREE/DSQRT(DBLE(5))
      KK=IHSH+IHSH-1
C      AA=-AA*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3)))
      AA=AA*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3)))
      RETURN
      END
