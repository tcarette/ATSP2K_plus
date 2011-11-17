*
*     -------------------------------------------------------------
*      S I X J 5
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
*                                                                  *
*     | J/2  K/2  L/2 |                                            *
*     | M/2  N/2  1/2 |             [B.M.X. 75]                    *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      SUBROUTINE SIXJ5(J,K,L,M,N,ITIK,SI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      SI=ZERO
      IF(ITIK.NE.0) THEN
C
C     CHESKED TRIANGULAR CONDITIONS
C
        IF(IXJTIK(J,K,L,M,N,1).EQ.0)RETURN
      ENDIF
      I1=(J+K+L)/2
      AS=DBLE(I1)
      A=DBLE(L)
      B=DBLE(K)
      C=DBLE(J)
      AKA=ONE
      IF(MOD(I1,2).NE.0)AKA=-AKA
      IF(K.LT.M) THEN
        IF(J.LT.N)  THEN
C              M > K,  J < N.
          SI=-AKA*DSQRT((AS+TWO)*(AS-A+ONE)/
     *      ((B+ONE)*(B+TWO)*(C+ONE)*(C+TWO)))
        ELSEIF(J.GT.N)  THEN
C              M > K,  J > N.
          SI=AKA*DSQRT((AS-C+ONE)*(AS-B)/
     *      ((B+ONE)*(B+TWO)*C*(C+ONE)))
        ENDIF
      ELSEIF(K.GT.M) THEN
        IF(J.LT.N) THEN
C             M < K,  J < N.
          SI=AKA*DSQRT((AS-C)*(AS-B+ONE)/
     *      (B*(B+ONE)*(C+ONE)*(C+TWO)))
        ELSEIF(J.GT.N) THEN
C             M < K,  J > N.
          SI=AKA*DSQRT((AS+ONE)*(AS-A)/(B*(B+ONE)*C*(C+ONE)))
        ENDIF
      ENDIF
      RETURN
      END
