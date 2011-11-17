*
*     ---------------------------------------------------------------
*     S U B L S
*     ---------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE SUBLS(J1,J2,LL,S)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      S=ZERO
      CALL RUMT(J1,LL,LQ,LS,L)
      CALL RUMT(J2,LL,LQS,LSS,L1S)
      QQ=HALF
      N=2
      Q=HALF*DBLE(LQ)
      QS=HALF*DBLE(LQS)
      QM=-HALF*DBLE(2*LL+1-N)
      QMS=-HALF*DBLE(2*LL+1-(N-1))
      CALL C0T5S(QS,QMS,QQ,Q,QM,A4)
      IF(DABS(A4).LT.EPS) RETURN
      A1=DBLE(N*(LQ+1)*(L+1)*(LS+1))
      S=DSQRT(A1)/A4
      IN=-N-1
      IF((IN/2)*2.NE.IN)S=-S
      RETURN
      END
