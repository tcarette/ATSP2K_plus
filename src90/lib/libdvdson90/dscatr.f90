 
      SUBROUTINE DSCATR(N, X, INCX, INDEX, INCI, Y, INCY) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER  :: INCX 
      INTEGER  :: INCI 
      INTEGER  :: INCY 
      INTEGER , INTENT(IN) :: INDEX(*) 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(OUT) :: Y(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!*******
!       Assume incx=inci=1=incy
!*******
 
      DO I = 1, N 
         Y(INDEX(I)) = X(I) 
      END DO 
      RETURN  
      END SUBROUTINE DSCATR 
