*
*     -------------------------------------------------------------
*      S I X J 1
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
*                                                                  *
*     | I/2  J/2  K/2 |                                            *
*     | L/2  M/2   1  |               [B.M.X.  75].                *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      SUBROUTINE SIXJ1(I,J,K,L,M,ITIK,SI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      SI=ZERO
      IF(ITIK.NE.0) THEN
C
C     CHESKED TRIANGULAR CONDITIONS
C
        IF(IXJTIK(I,J,K,L,M,2).EQ.0)RETURN
      ENDIF
      IFA=(I+J+K)/2
      AS=DBLE(IFA)
      AKA=ONE
      IF(MOD(IFA,2).NE.0)AKA=-AKA
      A=DBLE(K)
      B=DBLE(J)
      C=DBLE(I)
      IF(I.LT.M) THEN
        IF(J.LT.L) THEN
C              M > I,   L > J.
          SI=AKA*DSQRT((AS+TWO)*(AS+THREE)*(AS-A+ONE)*(AS-A+TWO)/
     *      ((B+ONE)*(B+TWO)*(B+THREE)*(C+ONE)*(C+TWO)*(C+THREE)))
        ELSEIF(J.EQ.L) THEN
C              M > I,  L = J.
          SI=(-AKA)*DSQRT(TWO*(AS+TWO)*(AS-C)*(AS-B+ONE)*(AS-A+ONE)/
     *      (B*(B+ONE)*(B+TWO)*(C+ONE)*(C+TWO)*(C+THREE)))
        ELSE
C              M > I,  L < J.
          SI=AKA*DSQRT((AS-C-ONE)*(AS-C)*(AS-B+ONE)*(AS-B+TWO)/
     *      ((B-ONE)*B*(B+ONE)*(C+ONE)*(C+TWO)*(C+THREE)))
        ENDIF
      ELSEIF(I.EQ.M) THEN
        IF(J.LT.L) THEN
C              M = L,  L > J.
          SI=(-AKA)*DSQRT((AS+TWO)*(AS-C+ONE)*(AS-B)*(AS-A+ONE)
     *      *TWO/((B+ONE)*(B+TWO)*(B+THREE)*C*(C+ONE)*(C+TWO)))
        ELSEIF(J.EQ.L) THEN
C              M = I,  L = J.
          SI=(-AKA)*((B*B+C*C-A*A)*HALF+B+C-A)/
     *      DSQRT(B*(B+ONE)*(B+TWO)*C*(C+ONE)*(C+TWO))
        ELSE
C              M = I,  L < J.
          SI=AKA*DSQRT((AS+ONE)*(AS-C)*(AS-B+ONE)*(AS-A)*TWO/
     *      ((B-ONE)*B*(B+ONE)*C*(C+ONE)*(C+TWO)))
        ENDIF
      ELSE
        IF(J.LT.L) THEN
C              M < I,   L > J.
          SI=AKA*DSQRT((AS-C+ONE)*(AS-C+TWO)*(AS-B-ONE)*(AS-B)/
     *      ((B+ONE)*(B+TWO)*(B+THREE)*(C-ONE)*C*(C+ONE)))
        ELSEIF(J.EQ.L) THEN
C              M < I,   L = J.
          SI=AKA*DSQRT((AS+ONE)*(AS-C+ONE)*(AS-B)*(AS-A)*TWO/
     *      (B*(B+ONE)*(B+TWO)*(C-ONE)*C*(C+ONE)))
        ELSE
C              M < I,   L < J.
          SI=AKA*DSQRT(AS*(AS+ONE)*(AS-A-ONE)*(AS-A)/
     *      ((B-ONE)*B*(B+ONE)*(C-ONE)*C*(C+ONE)))
        ENDIF
      ENDIF
      RETURN
      END
