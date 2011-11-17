 
!
!-----------------------------------------------------------------------
!     A C U R A T
!-----------------------------------------------------------------------
!
!     Coefficients of Slater integrals are square roots of
!   rational numbers.  To improve the accuaracy, certain commonly
!   occuring coefficients are improved to machine accuracy.
!
      REAL(KIND(0.0D0)) FUNCTION ACURAT (C) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: C 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NUM 
      INTEGER , DIMENSION(11) :: DEN 
      INTEGER :: I 
      REAL(DOUBLE) :: D1, C2, PROD, EPS 
!-----------------------------------------------
      DATA DEN/ 2, 3, 7, 9, 15, 35, 49, 175, 189, 315, 441/  
      DATA D1/ 1.D0/  
!
      C2 = C*C 
      ACURAT = C 
      IF (C < 0.) THEN 
         DO I = 1, 11 
            PROD = DEN(I)*C2 
            NUM = NINT(PROD) 
            EPS = DABS(NUM - PROD)/DEN(I) 
            IF (EPS > 1.D-8) CYCLE  
            IF (EPS == 0.) CYCLE  
            ACURAT = DSQRT((NUM*D1)/DEN(I)) 
            ACURAT = -ACURAT 
            RETURN  
         END DO 
      ELSE 
         DO I = 1, 11 
            PROD = DEN(I)*C2 
            NUM = NINT(PROD) 
            EPS = DABS(NUM - PROD)/DEN(I) 
            IF (EPS > 1.D-8) CYCLE  
            IF (EPS == 0.) CYCLE  
            ACURAT = DSQRT((NUM*D1)/DEN(I)) 
            RETURN  
         END DO 
      ENDIF 
      RETURN  
      END FUNCTION ACURAT 
