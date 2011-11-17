 
 
!
!     -------------------------------------------------------------
!      N I N E 1 2
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2      J2/2  J3/2 |                                    *
!     |  L1/2      L2/2  J3/2 |                                    *
!     |  K1/2 + 1  K1/2    1  |                                    *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE NINE12(J1, J2, J3, L1, L2, K1, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:09:06  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J1 
      INTEGER  :: J2 
      INTEGER  :: J3 
      INTEGER  :: L1 
      INTEGER  :: L2 
      INTEGER  :: K1 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K2 
      REAL(DOUBLE) :: S1, SKAI1, S2, SKAI2, S, VARD 
!-----------------------------------------------
      K2 = K1 - 2 
      CALL SIXJ (J1, L1, K1, L2, J2, J3, 0, S1) 
      SKAI1 = DBLE(J2 - L2 + K1)*DBLE(L2 - J2 + K1)*DBLE(J2 + L2 + K2 + 4)*&
         DBLE(J2 + L2 - K2) 
      S1 = S1*DSQRT(SKAI1) 
      CALL SIXJ (J1, L1, K2, L2, J2, J3, 0, S2) 
      SKAI2 = DBLE(J1 - L1 + K1)*DBLE(L1 - J1 + K1)*DBLE(J1 + L1 + K2 + 4)*&
         DBLE(J1 + L1 - K2) 
      S2 = S2*DSQRT(SKAI2) 
      S = S1 + S2 
!      VARD=DSQRT(DBLE(8*(K2+3)*(K2+2)*(K2+1)*J3*(J3+2)*(J3+1)))
      VARD = DSQRT(DBLE(8*(K2 + 3))*DBLE(K2 + 2)*DBLE(K2 + 1)*DBLE(J3)*DBLE(J3&
          + 2)*DBLE(J3 + 1)) 
      A = S/VARD 
      IF (MOD(L1 + J2 + K2 + J3,4) /= 0) A = -A 
      RETURN  
      END SUBROUTINE NINE12 
