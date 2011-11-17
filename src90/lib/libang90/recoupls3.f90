!
!     --------------------------------------------------------------
!     R E C O U P L S 3
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE RECOUPLS3(K, JA1, JA2, JA3, K1, K2, K3, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:00:59  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls31_I 
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
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE)  :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IFAZ 
!-----------------------------------------------
      IF (JA3>JA1 .AND. JA3>JA2) THEN 
         IF (JA1 - JA2 < 0) THEN 
            CALL RECOUPLS31 (K, JA1, JA2, JA3, K1, K2, K3, IRE, IAT, REC) 
         ELSE IF (JA1 - JA2 > 0) THEN 
            CALL RECOUPLS31 (K, JA2, JA1, JA3, K2, K1, K3, IRE, IAT, REC) 
            IFAZ = K1 + K2 - K3 
            IF ((IFAZ/4)*4 /= IFAZ) REC = -REC 
         ELSE 
            STOP  
         ENDIF 
      ELSE IF (JA3<JA1 .AND. JA3<JA2) THEN 
         IF (JA1 - JA2 < 0) THEN 
            CALL RECOUPLS31 (K, JA3, JA1, JA2, K3, K1, K2, IRE, IAT, REC) 
            IF ((K3/2)*2 /= K3) REC = -REC 
         ELSE IF (JA1 - JA2 > 0) THEN 
            CALL RECOUPLS31 (K, JA3, JA2, JA1, K3, K2, K1, IRE, IAT, REC) 
            IFAZ = K1 + K2 + K3 
            IF ((IFAZ/4)*4 /= IFAZ) REC = -REC 
         ELSE 
            STOP  
         ENDIF 
      ELSE 
         IF (JA1 - JA2 < 0) THEN 
            CALL RECOUPLS31 (K, JA1, JA3, JA2, K1, K3, K2, IRE, IAT, REC) 
            IFAZ = K1 - K2 - K3 
            IF ((IFAZ/4)*4 /= IFAZ) REC = -REC 
         ELSE IF (JA1 - JA2 > 0) THEN 
            CALL RECOUPLS31 (K, JA2, JA3, JA1, K2, K3, K1, IRE, IAT, REC) 
            IF ((K1/2)*2 /= K1) REC = -REC 
         ELSE 
            STOP  
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE RECOUPLS3 
