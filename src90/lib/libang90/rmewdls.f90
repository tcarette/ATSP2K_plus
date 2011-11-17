!
!     -----------------------------------------------------------------
!      R M E W D L S
!     -----------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                              June 1995   *
!
      SUBROUTINE RMEWDLS(J1, J2, K1, K2, K3, COEF) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RIBOLS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:13:23  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumt_I 
      USE rmewd1ls_I 
      USE rmewd2ls_I 
      USE rmewd3ls_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J1 
      INTEGER  :: J2 
      INTEGER , INTENT(IN) :: K1 
      INTEGER , INTENT(IN) :: K2 
      INTEGER , INTENT(IN) :: K3 
      REAL(DOUBLE)  :: COEF 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(16) :: IPR 
      INTEGER :: JI1, JI2, LQ, LS, LL, LQS, LSS, LLS, IFAZ, L, J, IFAZ2 
!-----------------------------------------------
      DATA IPR/ 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120/  
      COEF = 0.00 
      JI1 = IMPTLS(J1) 
      JI2 = IMPTLS(J2) 
      IF (JI1 /= JI2) RETURN  
      CALL RUMT (J1, 2, LQ, LS, LL) 
      CALL RUMT (J2, 2, LQS, LSS, LLS) 
      IF (J1 > J2) THEN 
         JI1 = J2 
         JI2 = J1 
         IFAZ = LQ + LS + LL - LQS - LSS - LLS 
      ELSE 
         JI1 = J1 
         JI2 = J2 
         IFAZ = 4 
      ENDIF 
      IF (J1 > 24) THEN 
         JI1 = JI1 - 24 
         JI2 = JI2 - 24 
         L = 0 
      ELSE 
         JI1 = JI1 - 8 
         JI2 = JI2 - 8 
         L = 1 
      ENDIF 
      J = IPR(JI2) + JI1 
      IF (K1==1 .AND. K2==2 .AND. K3==0) THEN 
         CALL RMEWD1LS (J, L, COEF) 
      ELSE IF (K1==0 .AND. K2==2 .AND. K3==1) THEN 
         L = IABS(L - 1) 
         CALL RMEWD1LS (J, L, COEF) 
         IF (L == 0) THEN 
            IFAZ2 = LQ - LQS 
         ELSE 
            IFAZ2 = LS - LSS 
         ENDIF 
         IF (MOD(IFAZ2,4) /= 0) COEF = -COEF 
      ELSE IF (K1==1 .AND. K2==1 .AND. K3==1) THEN 
         IF (L == 0) THEN 
            CALL RMEWD2LS (J, 0, COEF) 
         ELSE 
            CALL RMEWD2LS (J, 0, COEF) 
            IF (MOD(LQ - LQS,4) /= 0) COEF = -COEF 
         ENDIF 
      ELSE IF (K1==0 .AND. K2==3 .AND. K3==0) THEN 
         CALL RMEWD2LS (J, 1, COEF) 
      ELSE IF (K1==1 .AND. K2==3 .AND. K3==1) THEN 
         IF (L == 0) THEN 
            CALL RMEWD2LS (J, 2, COEF) 
         ELSE 
            CALL RMEWD2LS (J, 2, COEF) 
            IF (MOD(LQ - LQS,4) /= 0) COEF = -COEF 
         ENDIF 
      ELSE IF (K1==1 .AND. K2==4 .AND. K3==0) THEN 
         CALL RMEWD3LS (J, L, COEF) 
      ELSE IF (K1==0 .AND. K2==4 .AND. K3==1) THEN 
         L = IABS(L - 1) 
         CALL RMEWD3LS (J, L, COEF) 
         IF (L == 0) THEN 
            IFAZ2 = LQ - LQS 
         ELSE 
            IFAZ2 = LS - LSS 
         ENDIF 
         IF (MOD(IFAZ2,4) /= 0) COEF = -COEF 
      ENDIF 
      IF (MOD(IFAZ,4) /= 0) COEF = -COEF 
      RETURN  
      END SUBROUTINE RMEWDLS 
