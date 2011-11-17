!
!     -------------------------------------------------------------
!      I T T K
!     -------------------------------------------------------------
!
!                                                                  *
!     CHESKED TRIANGULAR CONDITIONS FOR   I/2, J/2, K/2.           *
!     I+J>=K, I+K>=J, J+K>=I,                                      *
!     I/2+J/2+K/2 - WHOLE NUMBER                                   *
!     ITTK=1 -   IF NOT SATISFY                                    *
!     ITTK=0 -   IN OVER CASES                                     *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      INTEGER FUNCTION ITTK (I, J, K) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:51:45  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: J 
      INTEGER , INTENT(IN) :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      ITTK = 0 
      IF (IABS(I - J) > K) RETURN  
      IF (I + J < K) RETURN  
      IF (MOD(I + J + K,2) /= 0) RETURN  
      ITTK = 1 
      RETURN  
      END FUNCTION ITTK 
