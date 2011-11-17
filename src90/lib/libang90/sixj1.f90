!
!     -------------------------------------------------------------
!      S I X J 1
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | I/2  J/2  K/2 |                                            *
!     | L/2  M/2   1  |               [B.M.X.  75].                *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                             March 1995   *
!
      SUBROUTINE SIXJ1(I, J, K, L, M, ITIK, SI) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:06:26  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
      INTEGER  :: K 
      INTEGER  :: L 
      INTEGER  :: M 
      INTEGER , INTENT(IN) :: ITIK 
      REAL(DOUBLE) , INTENT(OUT) :: SI 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IFA 
      REAL(DOUBLE) :: AS, AKA, A, B, C 
!-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
!
!     CHESKED TRIANGULAR CONDITIONS
!
         IF (IXJTIK(I,J,K,L,M,2) == 0) RETURN  
      ENDIF 
      IFA = (I + J + K)/2 
      AS = DBLE(IFA) 
      AKA = ONE 
      IF (MOD(IFA,2) /= 0) AKA = -AKA 
      A = DBLE(K) 
      B = DBLE(J) 
      C = DBLE(I) 
      IF (I < M) THEN 
         IF (J < L) THEN 
!              M > I,   L > J.
            SI = AKA*DSQRT((AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A + TWO&
               )/((B + ONE)*(B + TWO)*(B + THREE)*(C + ONE)*(C + TWO)*(C + &
               THREE))) 
         ELSE IF (J == L) THEN 
!              M > I,  L = J.
            SI = (-AKA)*DSQRT(TWO*(AS + TWO)*(AS - C)*(AS - B + ONE)*(AS - A + &
               ONE)/(B*(B + ONE)*(B + TWO)*(C + ONE)*(C + TWO)*(C + THREE))) 
         ELSE 
!              M > I,  L < J.
            SI = AKA*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO&
               )/((B - ONE)*B*(B + ONE)*(C + ONE)*(C + TWO)*(C + THREE))) 
         ENDIF 
      ELSE IF (I == M) THEN 
         IF (J < L) THEN 
!              M = L,  L > J.
            SI = (-AKA)*DSQRT((AS + TWO)*(AS - C + ONE)*(AS - B)*(AS - A + ONE)&
               *TWO/((B + ONE)*(B + TWO)*(B + THREE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (J == L) THEN 
!              M = I,  L = J.
            SI = (-AKA)*((B*B + C*C - A*A)*HALF + B + C - A)/DSQRT(B*(B + ONE)*&
               (B + TWO)*C*(C + ONE)*(C + TWO)) 
         ELSE 
!              M = I,  L < J.
            SI = AKA*DSQRT((AS + ONE)*(AS - C)*(AS - B + ONE)*(AS - A)*TWO/((B&
                - ONE)*B*(B + ONE)*C*(C + ONE)*(C + TWO))) 
         ENDIF 
      ELSE 
         IF (J < L) THEN 
!              M < I,   L > J.
            SI = AKA*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - B - ONE)*(AS - B&
               )/((B + ONE)*(B + TWO)*(B + THREE)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (J == L) THEN 
!              M < I,   L = J.
            SI = AKA*DSQRT((AS + ONE)*(AS - C + ONE)*(AS - B)*(AS - A)*TWO/(B*(&
               B + ONE)*(B + TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE 
!              M < I,   L < J.
            SI = AKA*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B - ONE)*B*(&
               B + ONE)*(C - ONE)*C*(C + ONE))) 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE SIXJ1 
