 
 
!
!     -------------------------------------------------------------
!      N I N E 1 3
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2      J2/2  J3/2+1 |                                  *
!     |  L1/2      L2/2  J3/2   |                                  *
!     |  K1/2 + 1  K1/2    1    |                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius, LITHUANIA                              January 1997 *
!
      SUBROUTINE NINE13(J1, J2, J3, L1, L2, K1, A) 
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
      USE nine_I 
      USE sixj_I 
      USE aconst_I 
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
      REAL(DOUBLE)  :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K2, L3, IV, IN, IS 
      REAL(DOUBLE) :: S1, A1, A2, S2, A3, S3, S, VAR2 
!-----------------------------------------------
      A = ZERO 
      K2 = K1 - 2 
      L3 = J3 - 2 
      IV = J3*(J1 - L1)*(J1 + L1 + 2) - K1*(J1 - J2)*(J1 + J2 + 2) + K1*J3*(K2&
          - L3) 
      IF (IV == 0) THEN 
         CALL NINE (J1, J2, J3, L1, L2, L3, K1, K2, 2, 0, IN, A) 
      ELSE 
         IS = L3 - K2 
         IF (IS == 0) THEN 
            S1 = ZERO 
         ELSE 
            CALL SIXJ (J1, L1, K2, L2, J2, L3, 1, S1) 
            CALL ACONST (J1 + 1, L1, K2 + 1, J2, L3 + 1, A1) 
            S1 = S1*A1*DBLE(IS)/TWO 
         ENDIF 
         CALL ACONST (L1 + 1, J1, K2 + 1, L2, L3 + 1, A2) 
         IF (DABS(A2) < EPS) THEN 
            S2 = ZERO 
         ELSE 
            CALL SIXJ (J1, J2, J3, L2, L1, K2, 1, S2) 
            S2 = S2*A2*DBLE(J3)/TWO 
         ENDIF 
         CALL ACONST (J2 + 1, J1, L3 + 1, L2, K2 + 1, A3) 
         IF (DABS(A3) < EPS) THEN 
            S3 = ZERO 
         ELSE 
            CALL SIXJ (J1, J2, L3, L2, L1, K1, 1, S3) 
            S3 = S3*A3*DBLE(K1)/TWO 
         ENDIF 
         S = S1 + S2 - S3 
         IF (DABS(S) < EPS) RETURN  
         VAR2 = DBLE(K2 + 3)*DBLE(K2 + 2)*DBLE(K2 + 1)*DBLE(L3 + 3)*DBLE(L3 + 2&
            )*DBLE(L3 + 1) 
         A = S*DBLE(8)/(DBLE(IV)*DSQRT(VAR2)) 
         IF (MOD(J2 + L3 + L1 + K2,4) /= 0) A = -A 
      ENDIF 
      RETURN  
      END SUBROUTINE NINE13 
