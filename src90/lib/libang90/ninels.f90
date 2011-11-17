 
 
!
!     -------------------------------------------------------------
!      N I N E L S
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  L3/2 |                                        *
!     |  K1/2  K2/2  K3/2 |                                        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius, LITHUANIA                              January 1997 *
!
      SUBROUTINE NINELS(J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, A) 
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
      USE ittk_I 
      USE nine0_I 
      USE nine_I 
      USE nine11_I 
      USE nine12_I 
      USE nine13_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J1 
      INTEGER  :: J2 
      INTEGER  :: J3 
      INTEGER  :: L1 
      INTEGER  :: L2 
      INTEGER  :: L3 
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      INTEGER  :: I 
      INTEGER  :: IN 
      REAL(DOUBLE)  :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      A = ZERO 
      IF (I == 1) THEN 
         IN = 0 
         IF (ITTK(J1,J2,J3) == 0) RETURN  
         IF (ITTK(L1,L2,L3) == 0) RETURN  
         IF (ITTK(K1,K2,K3) == 0) RETURN  
         IF (ITTK(J1,L1,K1) == 0) RETURN  
         IF (ITTK(J2,L2,K2) == 0) RETURN  
         IF (ITTK(J3,L3,K3) == 0) RETURN  
         IN = 1 
         RETURN  
      ENDIF 
      IN = 1 
      IF (J1*J2*J3*L1*L2*L3*K1*K2*K3 == 0) THEN 
         CALL NINE0 (J1, J2, J3, L1, L2, L3, K1, K2, K3, A) 
      ELSE IF (K3 /= 2) THEN 
         CALL NINE (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, A) 
      ELSE IF (J3 == L3) THEN 
         IF (K1 == K2) THEN 
!
!    Case  J3 = L3           and   K1 = K2
!
            CALL NINE11 (J1, J2, J3, L1, L2, K1, A) 
         ELSE IF (K1 - 2 == K2) THEN 
!
!    Case  J3 = L3           and   K1/2 + 1 = K2/2
!
            CALL NINE12 (J1, J2, J3, L1, L2, K1, A) 
         ELSE IF (K1 + 2 == K2) THEN 
!
!    Case  J3 = L3           and   K1/2 = K2/2 + 1
!
            CALL NINE12 (J2, J1, J3, L2, L1, K2, A) 
            IF (MOD(J1 + J2 + J3 + L1 + L2 + L3 + K1 + K2 + K3,4) /= 0) A = -A 
         ELSE 
            CALL NINE (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, A) 
         ENDIF 
      ELSE IF (J3 - 2 == L3) THEN 
         IF (K1 == K2) THEN 
!
!    Case  J3/2 + 1 = L3/2   and   K1 = K2
!
            CALL NINE12 (J1, L1, K1, J2, L2, J3, A) 
         ELSE IF (K1 - 2 == K2) THEN 
!
!    Case  J3/2 + 1 = L3/2   and   K1/2 + 1 = K2/2
!
            CALL NINE13 (J1, J2, J3, L1, L2, K1, A) 
         ELSE IF (K1 + 2 == K2) THEN 
!
!    Case  J3/2 + 1 = L3/2   and   K1/2 - 1 = K2/2
!
            CALL NINE13 (J2, J1, J3, L2, L1, K2, A) 
            IF (MOD(J1 + J2 + J3 + L1 + L2 + L3 + K1 + K2 + K3,4) /= 0) A = -A 
         ELSE 
            CALL NINE (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, A) 
         ENDIF 
      ELSE IF (J3 + 2 == L3) THEN 
         IF (K1 == K2) THEN 
!
!    Case  J3/2 = L3/2 + 1   and   K1 = K2
!
            CALL NINE12 (L1, J1, K1, L2, J2, L3, A) 
            IF (MOD(J1 + J2 + J3 + L1 + L2 + L3 + K1 + K2 + K3,4) /= 0) A = -A 
         ELSE IF (K1 - 2 == K2) THEN 
!
!    Case  J3/2 = L3/2 + 1   and   K1/2 - 1 = K2/2
!
            CALL NINE13 (L1, L2, L3, J1, J2, K1, A) 
            IF (MOD(J1 + J2 + J3 + L1 + L2 + L3 + K1 + K2 + K3,4) /= 0) A = -A 
         ELSE IF (K1 + 2 == K2) THEN 
!
!    Case  J3/2 = L3/2 + 1   and   K1/2 - 1 = K2/2
!
            CALL NINE13 (L2, L1, L3, J2, J1, K2, A) 
         ELSE 
            CALL NINE (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, A) 
         ENDIF 
      ELSE 
         CALL NINE (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, A) 
      ENDIF 
      RETURN  
      END SUBROUTINE NINELS 
