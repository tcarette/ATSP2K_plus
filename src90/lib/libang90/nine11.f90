 
 
!
!     -------------------------------------------------------------
!      N I N E 1 1
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  J3/2 |                                        *
!     |  K1/2  K1/2    1  |                                        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE NINE11(J1, J2, J3, L1, L2, K1, A) 
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
      INTEGER :: ISKAI 
      REAL(DOUBLE) :: S, VARD 
!-----------------------------------------------
      A = ZERO 
      ISKAI = DBLE((J1 - L1)*(J1 + L1 + 2) - (J2 - L2)*(J2 + L2 + 2)) 
      IF (ISKAI == 0) RETURN  
      CALL SIXJ (J1, L1, K1, L2, J2, J3, 0, S) 
      VARD = DSQRT(DBLE(4*K1)*DBLE(K1 + 2)*DBLE(K1 + 1)*DBLE(J3)*DBLE(J3 + 2)*&
         DBLE(J3 + 1)) 
      A = S*DBLE(ISKAI)/VARD 
      IF (MOD(L1 + J2 + K1 + J3,4) /= 0) A = -A 
      RETURN  
      END SUBROUTINE NINE11 
