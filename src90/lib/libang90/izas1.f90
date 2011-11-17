!
!     -------------------------------------------------------------
!      I Z A S 1
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                             March 1995   *
!
      INTEGER FUNCTION IZAS1 (IB, QB, IK, QK) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:51:45  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IB 
      INTEGER , INTENT(IN) :: IK 
      REAL(DOUBLE) , INTENT(IN) :: QB 
      REAL(DOUBLE) , INTENT(IN) :: QK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQB, IQK 
!-----------------------------------------------
      IZAS1 = 0 
      IQB = TWO*ABS(QB) + TENTH 
      IF (IQB > IB) RETURN  
      IF (MOD(IB + IQB,2) /= 0) RETURN  
      IQK = TWO*ABS(QK) + TENTH 
      IF (IQK > IK) RETURN  
      IF (MOD(IK + IQK,2) /= 0) RETURN  
      IZAS1 = 1 
      RETURN  
      END FUNCTION IZAS1 
