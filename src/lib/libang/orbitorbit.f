*
*     -------------------------------------------------------------
*      O R B I T O R B I T
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF  ORBIT - ORBIT  INTERACTION                               *
*                                                                  *
*                              2                                   *
*     (n l L S  n l L S ::   -a   Sqrt(2k+1)  (2k-1)/(k+1)         *
*       1 1 1 1  2 2 2 2                                           *
*                                                                  *
*           (k-1) (1)(k)   (k-1) (k)(1)(k)(0)       k-1  k+2       *
*        [[[C   * L  ] * [[C   * L  ] ]   ]        r  / r          *
*              1    1        2    2                 <    >         *
*                                                                  *
*                                     ::n l L S  n l L S )         *
*                                        3 3 3 3  4 4 4 4          *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1997   * 
*
      SUBROUTINE ORBITORBIT(L1,L2,L3,L4,KL,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      AA=ZERO
      IF(ITTK(L1,L3,KL).EQ.0)RETURN
      IF(ITTK(L2,L4,KL).EQ.0)RETURN
      AA=RME (L1,L3,KL)
      IF(DABS(AA).LT.EPS)RETURN
      AA=AA*RME (L2,L4,KL)
      IF(DABS(AA).LT.EPS)RETURN
      K=KL-1
      S=DBLE(2*K+1)*DBLE(L1+L3+K+2)*DBLE(L1+L3-K)*DBLE(L1-L3+K+1)
     :*DBLE(L3-L1+K+1)*DBLE(L2+L4+K+2)*DBLE(L2+L4-K)*DBLE(L2-L4+K+1)
     :*DBLE(L4-L2+K+1)
      AA=AA*TWO*DSQRT(S)/DBLE(K*KL)
      RETURN
      END
