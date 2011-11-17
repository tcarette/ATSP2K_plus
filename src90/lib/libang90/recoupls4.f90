!
!     --------------------------------------------------------------
!     R E C O U P L S 4
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE RECOUPLS4(K, JA1, JA2, JA3, JA4, K1, K2, K3, K4, KA, IRE, IAT&
         , REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:00:59  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dlsa3_I 
      USE dlsa2_I 
      USE dlsa4_I 
      USE dlsa1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: JA1 
      INTEGER  :: JA2 
      INTEGER  :: JA3 
      INTEGER  :: JA4 
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      INTEGER  :: K4 
      INTEGER  :: KA 
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA4, IB4, ISKR 
      REAL(DOUBLE) :: S1, S2, S3, S4, S, RE 
!-----------------------------------------------
      S1 = DBLE(J1QN1(JA1,K)) 
      S2 = DBLE(J1QN1(JA2,K)) 
      S3 = DBLE(J1QN1(JA3,K)) 
      S4 = DBLE(J1QN1(JA4,K)) 
      S = S1*S2*S3*S4 
      REC = ONE/DSQRT(S) 
      IA4 = J1QN1(JA4,K) - 1 
      IB4 = J1QN2(JA4,K) - 1 
      REC = REC*DSQRT(DBLE(IA4 + 1))/DSQRT(DBLE((K4 + 1)*(IB4 + 1))) 
!
      ISKR = JA3 - JA2 
      IF (ISKR > 1) THEN 
         IAT = 0 
         CALL DLSA3 (K, JA2, JA3, KA, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
!
      ISKR = JA4 - JA3 
      IF (ISKR > 1) THEN 
         IAT = 0 
         CALL DLSA3 (K, JA3, JA4, K4, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
!
      IAT = 0 
      CALL DLSA2 (K, JA1, JA4, K4, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
!
      IAT = 0 
      CALL DLSA4 (K, JA1, JA2, K1, K2, KA, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
!
      IAT = 0 
      CALL DLSA4 (K, JA2, JA3, KA, K3, K4, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      IF (JA1==1 .AND. JA2==2) RETURN  
!
      IAT = 0 
      CALL DLSA1 (K, JA1, K1, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
!
      ISKR = JA2 - JA1 
      IF (JA1 == 1) ISKR = JA2 - 1 - JA1 
      IF (ISKR <= 1) RETURN  
      IAT = 0 
      CALL DLSA3 (K, JA1, JA2, K1, IRE, IAT, RE) 
      REC = RE*REC 
      RETURN  
      END SUBROUTINE RECOUPLS4 
