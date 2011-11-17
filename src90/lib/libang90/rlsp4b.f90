!
!     --------------------------------------------------------------
!     R L S P 4 B
!     --------------------------------------------------------------
!                                                                  *
!
      SUBROUTINE RLSP4B(K, JA3, JA4, K5, K6, KA, IRE, IAT, REC) 
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
      USE dlsa3_I 
      USE dlsa5_I 
      USE dlsa4_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: JA3 
      INTEGER  :: JA4 
      INTEGER  :: K5 
      INTEGER  :: K6 
      INTEGER  :: KA 
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISKR 
      REAL(DOUBLE) :: RE 
!-----------------------------------------------
      REC = ONE/DSQRT(DBLE(J1QN1(JA4,K))) 
!
      ISKR = JA4 - JA3 
      IF (ISKR > 1) THEN 
         IAT = 0 
         CALL DLSA3 (K, JA3, JA4, K5, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
!
      ISKR = IHSH - JA4 
      IF (ISKR > 1) THEN 
         IAT = 0 
         CALL DLSA3 (K, JA4, IHSH, KA, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
!
      IF (JA4 /= IHSH) THEN 
         IAT = 0 
         CALL DLSA5 (K, JA4, KA, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
      IF (JA3 == JA4) RETURN  
!
      IAT = 0 
      CALL DLSA4 (K, JA3, JA4, K5, K6, KA, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      RETURN  
      END SUBROUTINE RLSP4B 
