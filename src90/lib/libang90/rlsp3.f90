!
!
!     --------------------------------------------------------------
!     R L S P 3
!     --------------------------------------------------------------
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE RLSP3(K, JA1, JA2, JA3, K1, K2, K3, K4, KA, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:24:23  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rlsp31_I 
      USE rlsp32_I 
      USE itrexg_I 
      USE dlsa6_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: JA1 
      INTEGER  :: JA2 
      INTEGER  :: JA3 
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      INTEGER  :: K4 
      INTEGER  :: KA 
      REAL(DOUBLE)  :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT, IPL, IKKL, IGL, J12, I 
      REAL(DOUBLE) :: R, R1, R2 
!-----------------------------------------------
      REC = ZERO 
!
      IF (JA3>JA1 .AND. JA3>JA2) THEN 
         IF (JA1 - JA2 < 0) THEN 
! 1 2 3
            CALL RLSP31 (K, JA1, JA2, JA3, K1, K2, K3, K4, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            CALL RLSP32 (K, JA1, JA2, JA3, K1, K2, K3, K4, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            CALL RLSP31 (K, JA1, JA2, JA3, K1, K2, K3, K4, KA, 1, IAT, REC) 
            CALL RLSP32 (K, JA1, JA2, JA3, K1, K2, K3, K4, KA, 1, IAT, R1) 
            REC = REC*R1 
         ELSE IF (JA1 - JA2 > 0) THEN 
! 2 1 3
            CALL RLSP31 (K, JA2, JA1, JA3, K2, K1, K3, K4, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            CALL RLSP32 (K, JA2, JA1, JA3, K2, K1, K3, K4, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            CALL RLSP31 (K, JA2, JA1, JA3, K2, K1, K3, K4, KA, 1, IAT, REC) 
            CALL RLSP32 (K, JA2, JA1, JA3, K2, K1, K3, K4, KA, 1, IAT, R1) 
            REC = REC*R1 
            IF (MOD(K1 + K2 - K3,4) /= 0) REC = -REC 
         ELSE 
            STOP  
         ENDIF 
      ELSE IF (JA3<JA1 .AND. JA3<JA2) THEN 
         IF (JA1 - JA2 < 0) THEN 
! 3 1 2
            IPL = ITREXG(K1,K4,K2,KA,IKKL) + 1 
            IF (IKKL == 0) RETURN  
            IGL = IPL + IKKL - 1 
            CALL RLSP31 (K, JA3, JA1, JA2, K4, K1, J12, K2, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            DO I = IPL, IGL, 2 
               J12 = I - 1 
               CALL DLSA6 (K, K4, K3, KA, K1, K2, J12, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA3, JA1, JA2, K4, K1, J12, K2, KA, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA3, JA1, JA2, K4, K1, J12, K2, KA, 1, IAT, R1) 
               CALL DLSA6 (K, K4, K3, KA, K1, K2, J12, 1, IAT, R2) 
               IF (MOD(2*K1 + K2 + K4 - J12 - K3,4) /= 0) R1 = -R1 
               REC = REC + R1*R2 
            END DO 
            IF (DABS(REC) < EPS) RETURN  
            CALL RLSP31 (K, JA3, JA1, JA2, K4, K1, J12, K2, KA, 1, IAT, R1) 
            REC = REC*R1 
         ELSE IF (JA1 - JA2 > 0) THEN 
! 3 2 1
            IPL = ITREXG(K2,K4,K1,KA,IKKL) + 1 
            IF (IKKL == 0) RETURN  
            IGL = IPL + IKKL - 1 
            CALL RLSP31 (K, JA3, JA2, JA1, K4, K2, J12, K1, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            DO I = IPL, IGL, 2 
               J12 = I - 1 
               CALL DLSA6 (K, K4, K3, KA, K2, K1, J12, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA3, JA2, JA1, K4, K2, J12, K1, KA, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA3, JA2, JA1, K4, K2, J12, K1, KA, 1, IAT, R1) 
               CALL DLSA6 (K, K4, K3, KA, K2, K1, J12, 1, IAT, R2) 
               IF (MOD(3*K2 + 2*K1 + K4 - J12 - 2*K3,4) /= 0) R1 = -R1 
               REC = REC + R1*R2 
            END DO 
            IF (DABS(REC) < EPS) RETURN  
            CALL RLSP31 (K, JA3, JA2, JA1, K4, K2, J12, K1, KA, 1, IAT, R1) 
            REC = REC*R1 
         ELSE 
            STOP  
         ENDIF 
      ELSE 
         IF (JA1 - JA2 < 0) THEN 
! 1 3 2
            IPL = ITREXG(K1,K4,K2,KA,IKKL) + 1 
            IF (IKKL == 0) RETURN  
            IGL = IPL + IKKL - 1 
            CALL RLSP31 (K, JA1, JA3, JA2, K1, K4, J12, K2, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            DO I = IPL, IGL, 2 
               J12 = I - 1 
               CALL DLSA6 (K, K4, K3, KA, K1, K2, J12, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA1, JA3, JA2, K1, K4, J12, K2, KA, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA1, JA3, JA2, K1, K4, J12, K2, KA, 1, IAT, R1) 
               CALL DLSA6 (K, K4, K3, KA, K1, K2, J12, 1, IAT, R2) 
               REC = REC + R1*R2 
            END DO 
            IF (DABS(REC) < EPS) RETURN  
            CALL RLSP31 (K, JA1, JA3, JA2, K1, K4, J12, K2, KA, 1, IAT, R1) 
            REC = REC*R1 
            IF (MOD(K1 + K2 - K3,4) /= 0) REC = -REC 
         ELSE IF (JA1 - JA2 > 0) THEN 
! 2 3 1
            IPL = ITREXG(K2,K4,K1,KA,IKKL) + 1 
            IF (IKKL == 0) RETURN  
            IGL = IPL + IKKL - 1 
            CALL RLSP31 (K, JA2, JA3, JA1, K2, K4, J12, K1, KA, 0, IAT, R) 
            IF (IAT == 0) RETURN  
            DO I = IPL, IGL, 2 
               J12 = I - 1 
               CALL DLSA6 (K, K4, K3, KA, K2, K1, J12, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA2, JA3, JA1, K2, K4, J12, K1, KA, 0, IAT, R) 
               IF (IAT == 0) CYCLE  
               CALL RLSP32 (K, JA2, JA3, JA1, K2, K4, J12, K1, KA, 1, IAT, R1) 
               CALL DLSA6 (K, K4, K3, KA, K2, K1, J12, 1, IAT, R2) 
               REC = REC + R1*R2 
            END DO 
            IF (DABS(REC) < EPS) RETURN  
            CALL RLSP31 (K, JA2, JA3, JA1, K2, K4, J12, K1, KA, 1, IAT, R1) 
            REC = REC*R1 
         ELSE 
            STOP  
         ENDIF 
      ENDIF 
      REC = REC/DSQRT(DBLE(J1QN1(JA1,K)*J1QN1(JA2,K)*J1QN1(JA3,K))) 
      RETURN  
      END SUBROUTINE RLSP3 
