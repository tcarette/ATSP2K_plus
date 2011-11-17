!
!     -----------------------------------------------------------------
!           F A C T R L
!     -----------------------------------------------------------------
!
!
      SUBROUTINE FACTRL(NFACT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE FACT_C 
!
!      GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NFACT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: ZERO, ONE, TWO, GAMMA, X 
!-----------------------------------------------
!
      DATA ZERO, ONE, TWO/ 0.D0, 1.D0, 2.D0/  
!
      GAMMA = ONE 
      GAM(1) = ZERO 
      DO I = 1, NFACT - 1 
         GAMMA = I*GAMMA 
         GAM(I+1) = DLOG(GAMMA) 
      END DO 
      DO I = NFACT + 1, 100 
         X = I - 1 
         GAM(I) = GAM(I-1) + DLOG(X) 
      END DO 
      RETURN  
      END SUBROUTINE FACTRL 
