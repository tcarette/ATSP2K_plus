!
!     --------------------------------------------------------------
!     R E C O U P L S 3 1
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE RECOUPLS31(K, JA1, JA2, JA3, K1, K2, KA, IRE, IAT, REC) 
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
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: KA 
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA3, IB3, ISKR 
      REAL(DOUBLE) :: S1, S2, S3, RE 
!-----------------------------------------------
      S1 = DBLE(J1QN1(JA1,K)) 
      S2 = DBLE(J1QN1(JA2,K)) 
      S3 = DBLE(J1QN1(JA3,K)) 
      REC = ONE/DSQRT(S1*S2*S3) 
      IA3 = J1QN1(JA3,K) - 1 
      IB3 = J1QN2(JA3,K) - 1 
      REC = REC*DSQRT(DBLE(IA3 + 1))/DSQRT(DBLE((KA + 1)*(IB3 + 1))) 
!
      IAT = 0 
      ISKR = JA3 - JA2 
      IF (ISKR > 1) THEN 
         CALL DLSA3 (K, JA2, JA3, KA, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
      IAT = 0 
      CALL DLSA2 (K, JA1, JA3, KA, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      IAT = 0 
      CALL DLSA4 (K, JA1, JA2, K1, K2, KA, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      IF (JA1==1 .AND. JA2==2) RETURN  
      IAT = 0 
      CALL DLSA1 (K, JA1, K1, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      ISKR = JA2 - JA1 
      IF (JA1 == 1) ISKR = JA2 - 1 - JA1 
      IF (ISKR <= 1) RETURN  
      IAT = 0 
      CALL DLSA3 (K, JA1, JA2, K1, IRE, IAT, RE) 
      REC = RE*REC 
      RETURN  
      END SUBROUTINE RECOUPLS31 
