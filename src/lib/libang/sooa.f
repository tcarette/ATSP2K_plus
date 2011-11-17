*
*     -------------------------------------------------------------
*      S O O A
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
      SUBROUTINE SOOA(L1,L2,L3,L4,KL,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      AA=ZERO
      IF(KL.LT.0)RETURN
      IF(ITTK(L1,L3,KL).EQ.0)RETURN
      IF(ITTK(L2,L4,KL).EQ.0)RETURN
      AA=RME(L1,L3,KL)
      AA=AA*RME(L2,L4,KL)
      KK=IHSH+IHSH-1
      AA=-AA*TWO*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3)))
      RETURN
      END
