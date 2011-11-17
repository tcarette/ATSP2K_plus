!
!     -----------------------------------------------------------------
!     N U M T E R F
!     -----------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      INTEGER FUNCTION NUMTERF (I2N, I2S, I2L, N, I2Q) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:32:21  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE jthn_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I2N 
      INTEGER , INTENT(IN) :: I2S 
      INTEGER , INTENT(IN) :: I2L 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: I2Q 
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!...  /MT67/ 
      COMMON /MT67/ M76(238) 
      INTEGER   M76 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISL, IA, J, LASTE, LS, NR 
!-----------------------------------------------
      NUMTERF = 0 
      ISL = I2S*100 + I2L 
      DO IA = 1, 119 
         IF ((N/2)*2 == N) THEN 
            J = 119 + IA 
         ELSE 
            J = IA 
         ENDIF 
         LASTE = M76(J) 
         LS = JTHN(LASTE,1,10000) 
         IF (ISL /= LS) CYCLE  
         NR = JTHN(LASTE,4,100) 
         IF (I2N == NR) GO TO 2 
      END DO 
      STOP  
    2 CONTINUE 
      NUMTERF = J 
      I2Q = JTHN(LASTE,3,100) 
!     WRITE(*,5) I2Q,NUMTERF
!   5 FORMAT(2X,'I2Q=',I6,'      J=',I6)
      RETURN  
      END FUNCTION NUMTERF 
