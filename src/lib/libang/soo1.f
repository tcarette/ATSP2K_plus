*
*     ..........................................................   :
*        iii)  Spin - other - orbit                                : 
*     ..........................................................   :
*
*     -------------------------------------------------------------
*      S O O 1
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
      SUBROUTINE SOO1(L1,KL1,KS1,KL2,KS2,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      AA=ZERO
      IF(KL2.LT.0)RETURN
      IF(ITTK(L1,L1,KL2).EQ.0)RETURN
      AA=RME (L1,L1,KL2)
      AA=AA*AA
C   s    space
      IF(KS1.LT.KS2)AA=TWO*AA
C   l    space
      IF(KL1.LT.KL2) THEN
        IS=(2*KL1+1)*(2*KL1+3)*(2*L1-KL1)*(KL1+1)*(KL1+2*L1+2)
      ELSEIF(KL1.EQ.KL2) THEN
        IS=KL1*(KL1+1)*(2*KL1+1)
      ELSE
        IS=(2*KL1+1)*(2*KL1-1)*(2*L1-KL1+1)*KL1*(KL1+2*L1+1)
      ENDIF
      AA=AA*DSQRT(DBLE(IS))*TWO
      KK=IHSH+IHSH-1
      AA=-AA*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3)))
      RETURN
      END
