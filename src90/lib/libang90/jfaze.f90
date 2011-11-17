!
!     -------------------------------------------------------------
!      J F A Z E
!     -------------------------------------------------------------
!                                                                  *
!     DETERMINATE THE PHASE FACTOR WHICH APPEAR FROM PERMUTATION   *
!     OF OPERATORS OF SECOND QUANTIZATION                          *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      INTEGER FUNCTION JFAZE (I1, I2, I3, I4) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:53:44  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I1 
      INTEGER , INTENT(IN) :: I2 
      INTEGER , INTENT(IN) :: I3 
      INTEGER , INTENT(IN) :: I4 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      JFAZE = 1 
      IF (I1 > I2) JFAZE = -JFAZE 
      IF (I1 > I3) JFAZE = -JFAZE 
      IF (I1 > I4) JFAZE = -JFAZE 
      IF (I2 > I3) JFAZE = -JFAZE 
      IF (I2 > I4) JFAZE = -JFAZE 
      IF (I3 > I4) JFAZE = -JFAZE 
      RETURN  
      END FUNCTION JFAZE 
