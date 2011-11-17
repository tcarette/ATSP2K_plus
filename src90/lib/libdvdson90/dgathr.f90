      SUBROUTINE DGATHR(N, X, INCX, INDEX, INCI, Y, INCY) 
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
      INTEGER , INTENT(IN) :: INCY 
      INTEGER , INTENT(IN) :: INDEX(*) 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(OUT) :: Y(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IY, I 
!-----------------------------------------------
 
!*******
!     Assume incx=inci=1
!*******
 
      IY = 1 
      Y(IY:(N-1)*INCY+IY:INCY) = X(INDEX(:N)) 
      RETURN  
      END SUBROUTINE DGATHR 
