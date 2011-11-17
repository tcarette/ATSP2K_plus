!     ------------------------------------------------------------------
!              H N O R M
!     ------------------------------------------------------------------
!
!       Returns the value of the normalization constant for an (nl)
!   hydrogenic function with nuclear charge ZZ.
!
!
      REAL(KIND(0.0D0)) FUNCTION HNORM (N, L, ZZ) 
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
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, I 
      REAL(DOUBLE) :: A, B, T, D 
!-----------------------------------------------
      M = L + L + 1 
      A = N + L 
      B = M 
      T = A 
      D = B 
      M = M - 1 
      IF (M /= 0) THEN 
         DO I = 1, M 
            A = A - D1 
            B = B - D1 
            T = T*A 
            D = D*B 
         END DO 
      ENDIF 
      HNORM = DSQRT(ZZ*T)/(N*D) 
      RETURN  
      END FUNCTION HNORM 
