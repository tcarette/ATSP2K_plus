*
*     -------------------------------------------------------------
*      S O O B
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
      SUBROUTINE SOOB(I1,L1,L2,KL1,KL2,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      AA=ZERO
      IF(KL1.LT.0)RETURN
      IF(ITTK(L1,L2,KL2).EQ.0)RETURN
      IF(KL1.LT.KL2) THEN
        AA=ONE
        IS1=(2*KL1+1)
        IS1=IS1*(2*KL1+3)*(L1+L2-KL1)*(KL1+1-L2+L1)*(KL1+1+L2-L1)*
     :             (KL1+L1+L2+2)
        IS2=(KL1+1)
      ELSEIF(KL1.EQ.KL2) THEN
        IF(I1.EQ.1) THEN
          IF(L1.EQ.L2) THEN
            AA=ONE
            IS1=KL1*(KL1+1)*(2*KL1+1)
            IS2=1
          ELSE
            IS1=(2*KL1+1)
            IS2=KL1*(KL1+1)
            IS3=(L2*(L2+1))+(KL1*(KL1+1))-(L1*(L1+1))
            AA=DBLE(IS3)
          ENDIF
        ELSE
          AA=TWO
          IS1=(2*KL1+1)*KL1*(KL1+1)
          IS2=1
        ENDIF
      ELSE
      IF(KL1.LT.1)RETURN
        AA=ONE
        IS1=(2*KL1+1)
        IS1=IS1*(2*KL1-1)*(L1+L2-KL1+1)*(KL1-L2+L1)*(KL1+L2-L1)*
     :             (KL1+L1+L2+1)
        IS2=KL1
      ENDIF
      AA=AA*DSQRT(DBLE(IS1)/DBLE(IS2))
      RETURN
      END
