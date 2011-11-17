!
!     -------------------------------------------------------------
!      I X J T I K
!     -------------------------------------------------------------
!
!                                                                  *
!     CHESKED TRIANGULAR CONDITIONS FOR 6j COEFFICIENT             *
!                                                                  *
!     | I/2  J/2  K/2 |            IXJTIK=1 - IF NOT SATISFY       *
!     | L/2  M/2  N/2 |            IXJTIK=0 - IN OVER CASES        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      INTEGER FUNCTION IXJTIK (I, J, K, L, M, N) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:51:45  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
      INTEGER  :: K 
      INTEGER  :: L 
      INTEGER  :: M 
      INTEGER  :: N 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      IXJTIK = 0 
      IF (ITTK(I,J,K) == 0) RETURN  
      IF (ITTK(I,M,N) == 0) RETURN  
      IF (ITTK(L,J,N) == 0) RETURN  
      IF (ITTK(L,M,K) == 0) RETURN  
      IXJTIK = 1 
      RETURN  
      END FUNCTION IXJTIK 
