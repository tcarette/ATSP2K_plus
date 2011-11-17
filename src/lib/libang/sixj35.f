*
*     -------------------------------------------------------------
*      S I X J 3 5
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
*                                                                  *
*     | J/2  K/2  L/2 |                                            *
*     | M/2  N/2  3/2 |             [B.M.X. 75]                    *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
*
      SUBROUTINE SIXJ35(J,K,L,M,N,ITIK,SI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      SI=ZERO
      IF(ITIK.NE.0) THEN
C
C     CHESKED TRIANGULAR CONDITIONS
C
        IF(IXJTIK(J,K,L,M,N,3).EQ.0)RETURN
      ENDIF
      I1=(J+K+L)/2
      AS=DBLE(I1)
      A=DBLE(L)
      B=DBLE(J)
      C=DBLE(K)
      AKA=ONE
      IF(MOD(I1,2).NE.0)AKA=-AKA
      IF(J-N.EQ.3) THEN
* -3
        IF(K-M.EQ.3) THEN
C  I                      -3/2  -3/2
      SI=AKA*DSQRT((AS-ONE)*AS*(AS+ONE)*(AS-A-TWO)*(AS-A-ONE)*(AS-A)/
     :((B-TWO)*(B-ONE)*B*(B+ONE)*(C-TWO)*(C-ONE)*C*(C+ONE)))
	ELSEIF(M-K.EQ.3) THEN
C  IV  P(12)              3/2   -3/2
      SI=AKA*DSQRT((AS-C-TWO)*(AS-C-ONE)*(AS-C)*(AS-B+ONE)*(AS-B+TWO)*
     :(AS-B+THREE)/ 
     :((C+1)*(C+TWO)*(C+THREE)*(C+FOUR)*(B-TWO)*(B-ONE)*B*(B+ONE)))
        ELSEIF(K-M.EQ.1) THEN
C  II  P(12)             -1/2   -3/2
      SI=AKA*DSQRT(THREE*AS*(AS+ONE)*(AS-A-ONE)*(AS-A)*(AS-C)*
     :(AS-B+ONE)/
     :((C-ONE)*C*(C+ONE)*(C+TWO)*(B-TWO)*(B-ONE)*B*(B+ONE)))
        ELSEIF(M-K.EQ.1) THEN
C  III P(12)              1/2   -3/2
      SI=AKA*DSQRT(THREE*(AS+ONE)*(AS-A)*(AS-C-ONE)*(AS-C)*
     :(AS-B+ONE)*(AS-B+TWO)/
     :(C*(C+ONE)*(C+TWO)*(C+THREE)*(B-TWO)*(B-ONE)*B*(B+ONE)))
        ENDIF
      ELSEIF(N-J.EQ.3) THEN
*  3
        IF(K-M.EQ.3) THEN
C  IV                     -3/2   3/2
      SI=AKA*DSQRT((AS-B-TWO)*(AS-B-ONE)*(AS-B)*(AS-C+ONE)*(AS-C+TWO)*
     :(AS-C+THREE)/
     :((B+ONE)*(B+TWO)*(B+THREE)*(B+FOUR)*(C-TWO)*(C-ONE)*C*(C+ONE)))
	ELSEIF(M-K.EQ.3) THEN
C  2       pataisyta               3/2   3/2
      SI=-AKA*DSQRT((AS+TWO)*(AS+THREE)*(AS+FOUR)*(AS-A+ONE)*
     :(AS-A+TWO)*(AS-A+THREE)/
     :((B+ONE)*(B+TWO)*(B+THREE)*(B+FOUR)*(C+ONE)*(C+TWO)*(C+THREE)*
     :(C+FOUR)))
        ELSEIF(K-M.EQ.1) THEN
C  1   P(12)   pataisytas          -1/2    3/2
      SI=-AKA*DSQRT(THREE*(AS+TWO)*(AS-A+ONE)*(AS-C+ONE)*
     : (AS-C+TWO)*(AS-B-ONE)*(AS-B)/
     :((C-ONE)*C*(C+ONE)*(C+TWO)*(B+ONE)*(B+TWO)*(B+THREE)*(B+FOUR)))
        ELSEIF(M-K.EQ.1) THEN
C  3  P(12)     taisyta           1/2    3/2
      SI=AKA*DSQRT(THREE*(AS+TWO)*(AS+THREE)*(AS-A+ONE)*(AS-A+TWO)*
     :(AS-B)*(AS-C+ONE)/
     :(C*(C+ONE)*(C+TWO)*(C+THREE)*(B+ONE)*(B+TWO)*(B+THREE)*(B+FOUR)))
        ENDIF
* -1
      ELSEIF(J-N.EQ.1) THEN
        IF(K-M.EQ.3) THEN
C  II                   -3/2   -1/2
      SI=AKA*DSQRT((THREE*AS*(AS+ONE)*(AS-A-ONE)*(AS-A)*(AS-B)*
     :(AS-C+ONE))/
     :((B-ONE)*B*(B+ONE)*(B+TWO)*(C-TWO)*(C-ONE)*C*(C+ONE)))
        ELSEIF(M-K.EQ.3) THEN
C  1                     3/2   -1/2
      SI=-AKA*DSQRT(THREE*(AS+TWO)*(AS-A+ONE)*(AS-B+ONE)*
     : (AS-B+TWO)*(AS-C-ONE)*(AS-C)/
     :((B-ONE)*B*(B+ONE)*(B+TWO)*(C+ONE)*(C+TWO)*(C+THREE)*(C+FOUR)))
        ELSEIF(K-M.EQ.1) THEN
C  V                    -1/2   -1/2
      SI=AKA*(TWO*(AS-B)*(AS-C)-(AS+TWO)*(AS-A-ONE))*
     :DSQRT((AS+ONE)*(AS-A)/
     :((B-ONE)*B*(B+ONE)*(B+TWO)*(C-ONE)*C*(C+ONE)*(C+TWO)))
        ELSEIF(M-K.EQ.1) THEN
C  VI P(12)              1/2   -1/2
      SI=AKA*((AS-B+TWO)*(AS-C+ONE)-TWO*(AS-A+ONE)*(AS+ONE))*
     :DSQRT((AS-C)*(AS-B+ONE)/
     :(C*(C+ONE)*(C+TWO)*(C+THREE)*(B-ONE)*B*(B+ONE)*(B+TWO)))
        ENDIF
* 1
      ELSEIF(N-J.EQ.1) THEN
        IF(K-M.EQ.3) THEN
C  III                  -3/2    1/2
      SI=AKA*DSQRT(THREE*(AS+ONE)*(AS-A)*(AS-B-ONE)*(AS-B)*(AS-C+ONE)*
     :(AS-C+TWO)/
     :(B*(B+ONE)*(B+TWO)*(B+THREE)*(C-TWO)*(C-ONE)*C*(C+ONE)))
        ELSEIF(M-K.EQ.3) THEN
C  3              pataisyta       3/2    1/2
      SI=AKA*DSQRT(THREE*(AS+TWO)*(AS+THREE)*(AS-A+ONE)*(AS-A+TWO)*
     :(AS-B+ONE)*(AS-C)/
     :(B*(B+ONE)*(B+TWO)*(B+THREE)*(C+ONE)*(C+TWO)*(C+THREE)*(C+FOUR)))
        ELSEIF(K-M.EQ.1) THEN
C  VI                   -1/2    1/2
      SI=AKA*((AS-C+TWO)*(AS-B+ONE)-TWO*(AS-A+ONE)*(AS+ONE))*
     :DSQRT((AS-B)*(AS-C+ONE)/
     :(B*(B+ONE)*(B+TWO)*(B+THREE)*(C-ONE)*C*(C+ONE)*(C+TWO)))
        ELSEIF(M-K.EQ.1) THEN
C  4      pataisyta               1/2    1/2
      SI=-AKA*(TWO*(AS-B)*(AS-C)-(AS+THREE)*(AS-A))*
     :DSQRT((AS+TWO)*(AS-A+ONE)/
     :(B*(B+ONE)*(B+TWO)*(B+THREE)*C*(C+ONE)*(C+TWO)*(C+THREE)))
      ENDIF
      ENDIF
      RETURN
      END
