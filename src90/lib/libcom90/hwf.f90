!
!     ------------------------------------------------------------------
!              H W F
!     ------------------------------------------------------------------
!
!       Returns the value of an unnormalized (nl) hydrogenic function
!   with nuclear charge ZZ and radius r.
!
!
      REAL(KIND(0.0D0)) FUNCTION HWF (N, L, ZZ, R) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: L 
      REAL(DOUBLE) , INTENT(IN) :: ZZ 
      REAL(DOUBLE) , INTENT(IN) :: R 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I 
      REAL(DOUBLE) :: P, A, B, C, X 
!-----------------------------------------------
      K = N - L - 1 
      P = D1 
      A = D1 
      B = K 
      C = N + L 
      X = -D2*ZZ*R/N 
!
!  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0
!
      IF (X >= (-150.D0)) THEN 
         IF (K >= 0) THEN 
            IF (K /= 0) THEN 
               DO I = 1, K 
                  P = D1 + A/B*P/C*X 
                  A = A + D1 
                  B = B - D1 
                  C = C - D1 
               END DO 
            ENDIF 
            HWF = P*DEXP(X/D2)*(-X)**(L + 1) 
            RETURN  
         ENDIF 
         WRITE (6, 7) N, L, ZZ, R 
    7    FORMAT(' FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM'/,' N =',&
            I4,'   L =',I4,'   Z =',F6.1,'   R =',F8.4) 
         STOP  
      ENDIF 
      HWF = D0 
      RETURN  
      END FUNCTION HWF 
